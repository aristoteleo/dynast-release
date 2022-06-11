from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import pystan

from .. import config, utils
from ..logging import logger


def read_pi(pi_path: str, group_by: Optional[List[str]] = None) -> Dict[Union[str, Tuple[str, ...]], float]:
    """Read pi CSV as a dictionary.

    Args:
        pi_path: path to CSV containing pi values
        group_by: columns that were used to group estimation

    Returns:
        Dictionary with barcodes and genes as keys
    """
    if group_by is None:
        with open(pi_path, 'r') as f:
            # pi is always the last column
            return float(f.readline().strip().split(',')[-1])

    df = pd.read_csv(pi_path, usecols=group_by + ['pi'], dtype={key: 'string' for key in group_by})
    return dict(df.set_index(group_by)['pi'])


# Process initializer.
# When creating a ProcessPoolExecutor, we specify an initializer, which performs
# some expensive calculation ONCE for each process. In our case, serialization of
# the STAN model takes a long time (~500ms), which we do not want to do every
# time the process performs a job, but just once when the process is created.
_model = None


def initializer(model: pystan.StanModel):
    """Multiprocessing initializer.
    https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ThreadPoolExecutor

    This initializer performs a one-time expensive initialization for each
    process.
    """
    global _model
    _model = model


def beta_mean(alpha: float, beta: float) -> float:
    """Calculate the mean of a beta distribution.
    https://en.wikipedia.org/wiki/Beta_distribution

    Args:
        alpha: First parameter of the beta distribution
        beta: Second parameter of the beta distribution

    Returns:
        Mean of the beta distribution
    """
    return alpha / (alpha + beta)


def beta_mode(alpha: float, beta: float) -> float:
    """Calculate the mode of a beta distribution.
    https://en.wikipedia.org/wiki/Beta_distribution

    When the distribution is bimodal (`alpha`, `beta` < 1), this function returns
    `nan`.

    Args:
        alpha: First parameter of the beta distribution
        beta: Second parameter of the beta distribution

    Returns:
        Mode of the beta distribution
    """
    pi = float('nan')

    # We ignore the following two cases:
    # If alpha=1 and beta=1, mode=any value in (0, 1)
    # If alpha<1 and beta<1, mode=0, 1 (bimodal)
    if alpha > 1 and beta > 1:
        pi = (alpha - 1) / (alpha + beta - 2)
    elif alpha > 1 and beta <= 1:
        pi = 1
    elif alpha <= 1 and beta > 1:
        pi = 0

    return pi


def guess_beta_parameters(guess: float, strength: int = 5) -> Tuple[float, float]:
    """Given a `guess` of the mean of a beta distribution, calculate beta
    distribution parameters such that the distribution is skewed by some
    `strength` toward the `guess`.

    Args:
        guess: Guess of the mean of the beta distribution
        strength: Strength of the skew

    Returns
        Beta distribution parameters (alpha, beta)
    """
    alpha = max(strength * guess, 0.1)
    beta = max(strength - alpha, 0.1)
    return alpha, beta


def fit_stan_mcmc(
        values: np.ndarray,
        p_e: float,
        p_c: float,
        guess: float = 0.5,
        model: pystan.StanModel = None,
        n_chains: int = 1,
        n_warmup: int = 1000,
        n_iters: int = 1000,
        seed: Optional[int] = None,
) -> Tuple[float, float, float, float]:
    """Run MCMC to estimate the fraction of labeled RNA.

    Args:
        values: Array of three columns encoding a sparse array in
            (row, column, value) format, zero-indexed, where
            row:    number of conversions
            column: nucleotide content
            value:  number of reads
        p_e: Average mutation rate in unlabeled RNA
        p_c: Average mutation rate in labeled RNA
        guess: Guess for the fraction of labeled RNA
        model: PyStan model to run MCMC with. If not provided, will try to use the
            `_model` global variable
        n_chains: Number of MCMC chains
        n_warmup: Number of warmup iterations
        n_iters: Number of MCMC iterations, excluding any warmups
    seed: random seed used for MCMC

    Returns:
        (guess, alpha, beta, pi)
    """
    model = model or _model

    # Skew beta distribution toward the guess.
    alpha_guess, beta_guess = guess_beta_parameters(guess)

    data = {
        'N': len(values),
        'conversions': values[:, 0],
        'contents': values[:, 1],
        'counts': values[:, 2],
        'p_c': p_c,
        'p_e': p_e,
    }
    init = [{'log_alpha': np.log(alpha_guess), 'log_beta': np.log(beta_guess), 'pi_g': guess}] * n_chains

    with utils.suppress_stdout_stderr():
        fit = model.sampling(
            data=data,
            pars=['alpha', 'beta', 'pi_g'],
            n_jobs=1,
            chains=n_chains,
            warmup=n_warmup,
            iter=n_iters + n_warmup,
            init=init,
            control={'adapt_delta': 0.99},
            seed=seed,
        )
    samples = fit.extract(('alpha', 'beta', 'pi_g'))
    alpha, beta, pi = np.mean(samples['alpha']), np.mean(samples['beta']), np.mean(samples['pi_g'])
    return guess, alpha, beta, pi


def estimate_pi(
        df_aggregates: pd.DataFrame,
        p_e: float,
        p_c: float,
        pi_path: str,
        group_by: Optional[List[str]] = None,
        p_group_by: Optional[List[str]] = None,
        n_threads: int = 8,
        threshold: int = 16,
        seed: Optional[int] = None,
        nasc: bool = False,
        model: Optional[pystan.StanModel] = None,
) -> str:
    """Estimate the fraction of labeled RNA.

    Args:
        df_aggregates: Pandas dataframe containing aggregate values
        p_e: Average mutation rate in unlabeled RNA
        p_c: Average mutation rate in labeled RNA
        pi_path: Path to write pi estimates
        group_by: Columns that were used to group cells
        p_group_by: Columns that p_e/p_c estimation was grouped by
        n_threads: Number of threads
        threshold: Any conversion-content pairs with fewer than this many reads
            will not be processed
        seed: Random seed
        nasc: Flag to change behavior to match NASC-seq pipeline. Specifically,
            the mode of the estimated Beta distribution is used as pi, defaults to `False`
        model: PyStan model to run MCMC with. If not provided, will try to compile the module manually

    Returns:
        Path to pi output
    """
    df_full = df_aggregates[(df_aggregates[['base', 'count']] > 0).all(axis=1)].copy()

    model = model or pystan.StanModel(file=config.MODEL_PATH, model_name=config.MODEL_NAME)

    if p_group_by is not None:
        df_full.set_index(p_group_by, inplace=True)
        df_full['p_e'] = df_full.index.map(p_e)
        df_full['p_c'] = df_full.index.map(p_c)
        df_full.reset_index(inplace=True)
    else:
        df_full['p_e'] = p_e
        df_full['p_c'] = p_c
    df_full.dropna(subset=['p_c'], inplace=True)  # Drop NA values due to p_c
    values = df_full[['conversion', 'base', 'count']].values
    p_es = df_full['p_e'].values
    p_cs = df_full['p_c'].values
    if group_by is None:
        # p_e, p_c must be a single float
        if not isinstance(p_e, float) or not isinstance(p_c, float):
            raise Exception('`p_e` and `p_c` must be floats when `group_by` is not provided')
        guess, alpha, beta, pi = fit_stan_mcmc(
            values,
            p_e,
            p_c,
            guess=0.5,
            seed=seed,
        )
        pi = beta_mode(alpha, beta) if nasc else pi
    else:
        groups = df_full.groupby(group_by, sort=False, observed=True).indices
        pis = {}
        skipped = 0
        failed = 0
        logger.debug(f'Spawning {n_threads} processes')
        # TODO: collapse rows with duplicate conversion and base counts
        # This can be done by constructing a csr_matrix and coverting it back
        with ProcessPoolExecutor(max_workers=n_threads, initializer=initializer, initargs=(model,)) as executor:
            futures = {}
            # Run larger groups first
            for key in sorted(groups.keys(), key=lambda key: len(groups[key]), reverse=True):
                idx = groups[key]
                p_e_unique = np.unique(p_es[idx])
                p_c_unique = np.unique(p_cs[idx])

                if len(p_e_unique) > 1 or len(p_c_unique) > 1:
                    raise Exception(
                        'p_e and p_c for each aggregate group must be a constant, '
                        f'but instead got {p_e_unique} and {p_c_unique}'
                    )

                p_e = p_e_unique[0]
                p_c = p_c_unique[0]
                vals = values[idx]
                count = sum(vals[:, 2])

                if count < threshold:
                    skipped += 1
                    continue

                # Make a naive guess of the fraction of new RNA
                # Clip guess to [0.01, 0.99] because STAN initialization will fail
                # if guess is either 0 or 1.
                guess = min(max(sum(vals[vals[:, 0] > 0][:, 2]) / count, 0.01), 0.99)

                futures[executor.submit(
                    fit_stan_mcmc,
                    vals,
                    p_e,
                    p_c,
                    guess=guess,
                    seed=seed,
                )] = key

            for future in utils.as_completed_with_progress(futures):
                key = futures[future]
                try:
                    guess, alpha, beta, pi = future.result()
                    pis[key] = (guess, alpha, beta, beta_mode(alpha, beta) if nasc else pi)
                except RuntimeError:
                    failed += 1

        if skipped > 0:
            logger.warning(
                f'Estimation skipped for {skipped} cell-genes because they have less than '
                f'{threshold} reads. Use `--cell-gene-threshold` to change.'
            )
        if failed > 0:
            logger.warning(f'Estimation failed {failed} times.')

    with open(pi_path, 'w') as f:
        if group_by is None:
            f.write(f'{guess},{alpha},{beta},{pi}')
        else:
            f.write(f'{",".join(group_by)},guess,alpha,beta,pi\n')
            for key in sorted(pis.keys()):
                guess, alpha, beta, pi = pis[key]
                formatted_key = key if isinstance(key, str) else ",".join(key)
                f.write(f'{formatted_key},{guess},{alpha},{beta},{pi}\n')

    return pi_path
