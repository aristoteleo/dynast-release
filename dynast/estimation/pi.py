from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pystan

from .. import config, utils
from ..logging import logger


def read_pi(pi_path, group_by=None):
    """Read pi CSV as a dictionary.

    :param pi_path: path to CSV containing pi values
    :type pi_path: str
    :param group_by: columns that were used to group estimation, defaults to
        ``None``
    :type group_by: list, optional

    :return: dictionary with barcodes and genes as keys
    :rtype: dictionary
    """
    if group_by is None:
        with open(pi_path, 'r') as f:
            # pi is always the last column
            return float(f.readline().strip().split(',')[-1])

    df = pd.read_csv(pi_path, usecols=group_by + ['pi'])
    return dict(df.set_index(group_by)['pi'])


# Process initializer.
# When creating a ProcessPoolExecutor, we specify an initializer, which performs
# some expensive calculation ONCE for each process. In our case, serialization of
# the STAN model takes a long time (~500ms), which we do not want to do every
# time the process performs a job, but just once when the process is created.
_model = None


def initializer(model):
    """Multiprocessing initializer.
    https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ThreadPoolExecutor

    This initializer performs a one-time expensive initialization for each
    process.
    """
    global _model
    _model = model


def beta_mean(alpha, beta):
    """Calculate the mean of a beta distribution.
    https://en.wikipedia.org/wiki/Beta_distribution

    :param alpha: first parameter of the beta distribution
    :type alpha: float
    :param beta: second parameter of the beta distribution
    :type beta: float

    :return: mean of the beta distribution
    :rtype: float
    """
    return alpha / (alpha + beta)


def beta_mode(alpha, beta):
    """Calculate the mode of a beta distribution.
    https://en.wikipedia.org/wiki/Beta_distribution

    When the distribution is bimodal (`alpha`, `beta` < 1), this function returns
    `nan`.

    :param alpha: first parameter of the beta distribution
    :type alpha: float
    :param beta: second parameter of the beta distribution
    :type beta: float

    :return: mode of the beta distribution
    :rtype: float
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


def guess_beta_parameters(guess, strength=5):
    """Given a `guess` of the mean of a beta distribution, calculate beta
    distribution parameters such that the distribution is skewed by some
    `strength` toward the `guess`.

    :param guess: guess of the mean of the beta distribution
    :type guess: float
    :param strength: strength of the skew, defaults to `5`
    :type strength: int

    :return: beta distribution parameters (alpha, beta)
    :rtype: (float, float)
    """
    alpha = max(strength * guess, 0.1)
    beta = max(strength - alpha, 0.1)
    return alpha, beta


def fit_stan_mcmc(
    values,
    p_e,
    p_c,
    guess=0.5,
    model=None,
    n_chains=1,
    n_warmup=1000,
    n_iters=1000,
    seed=None,
):
    """Run MCMC to estimate the fraction of labeled RNA.

    :param values: array of three columns encoding a sparse array in
                   (row, column, value) format, zero-indexed, where
                       row:    number of conversions
                       column: nucleotide content
                       value:  number of reads
    :type values: numpy.ndarray
    :param p_e: average mutation rate in unlabeled RNA
    :type p_e: float
    :param p_c: average mutation rate in labeled RNA
    :type p_c: float
    :param guess: guess for the fraction of labeled RNA, defaults to `0.5`
    :type guess: float, optional
    :param model: pyStan model to run MCMC with, defaults to `None`
                  if not provided, will try to use the `_model` global variable
    :type model: pystan.StanModel, optional
    :param n_chains: number of MCMC chains, defaults to `1`
    :type n_chains: int, optional
    :param n_warmup: number of warmup iterations, defaults to `1000`
    :type n_warmup: int, optional
    :param n_iters: number of MCMC iterations, excluding any warmups, defaults to `1000`
    :type n_iters: int, optional
    :param seed: random seed used for MCMC, defaults to `None`
    :type seed: int, optional

    :return: (guess, alpha, beta, pi)
    :rtype: (float, float, float, float)
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
    df_aggregates,
    p_e,
    p_c,
    pi_path,
    group_by=None,
    p_group_by=None,
    n_threads=8,
    threshold=16,
    seed=None,
    nasc=False,
    model=None,
):
    """Estimate the fraction of labeled RNA.

    :param df_aggregates: Pandas dataframe containing aggregate values
    :type df_aggregates: pandas.DataFrame
    :param p_e: average mutation rate in unlabeled RNA
    :type p_e: float
    :param p_c: average mutation rate in labeled RNA
    :type p_c: float
    :param pi_path: path to write pi estimates
    :type pi_path: str
    :param group_by: columns that were used to group cells, defaults to
        ``None``
    :type group_by: list, optional
    :param p_group_by: columns that p_e/p_c estimation was grouped by, defaults to `None`
    :type p_group_by: list, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional
    :param threshold: any conversion-content pairs with fewer than this many reads
                      will not be processed, defaults to `16`
    :type threshold: int, optional
    :param seed: random seed, defaults to `None`
    :type seed: int, optional
    :param nasc: flag to change behavior to match NASC-seq pipeline. Specifically,
                 the mode of the estimated Beta distribution is used as pi, defaults to `False`
    :type nasc: bool, optional
    :param model: pyStan model to run MCMC with, defaults to `None`
                  if not provided, will try to compile the module manually
    :type model: pystan.StanModel, optional

    :return: path to pi output
    :rtype: str
    """
    df_full = df_aggregates[(df_aggregates[['base', 'count']] > 0).all(axis=1)]

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
