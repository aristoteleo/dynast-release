import logging
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pystan
from scipy import sparse

from .. import config, utils

logger = logging.getLogger(__name__)


def read_pi(pi_path):
    df = pd.read_csv(pi_path, usecols=['barcode', 'GX', 'pi'])
    return dict(df.set_index(['barcode', 'GX'])['pi'])


# Process initializer.
# When creating a ProcessPoolExecutor, we specify an initializer, which performs
# some expensive calculation ONCE for each process. In our case, serialization of
# the STAN model takes a long time (~500ms), which we do not want to do every
# time the process performs a job, but just once when the process is created.
_model = None


def initializer(model):
    global _model
    _model = model


def beta_mean(alpha, beta):
    return alpha / (alpha + beta)


def beta_mode(alpha, beta):
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
    alpha = max(strength * guess, 0.1)
    beta = max(strength - alpha, 0.1)
    return alpha, beta


def fit_stan_mcmc(
    values,
    p_e,
    p_c,
    guess=0.5,
    model=None,
    pi_func=(lambda alpha, beta: None),
    n_chains=1,
    n_warmup=500,
    n_iters=1000,
    seed=None,
    subset_threshold=8000,
    subset_seed=None,
):
    model = model or _model
    conversions = []
    contents = []
    for k, n, count in values:
        conversions.extend([k] * count)
        contents.extend([n] * count)

    # If we have many conversions, randomly select a subset. This keeps runtimes
    # down while minimally impacting accuracy.
    if len(conversions) > subset_threshold:
        choices = np.random.RandomState(subset_seed).choice(len(conversions), subset_threshold, replace=False)
        conversions = list(np.array(conversions)[choices])
        contents = list(np.array(contents)[choices])

    # Skew beta distribution toward the guess.
    alpha_guess, beta_guess = guess_beta_parameters(guess)
    alpha_guess_log = np.log(alpha_guess)
    beta_guess_log = np.log(beta_guess)

    data = {
        'N': len(conversions),
        'contents': contents,
        'conversions': conversions,
        'p_c': p_c,
        'p_e': p_e,
    }
    init = [{'log_alpha': alpha_guess_log, 'log_beta': beta_guess_log, 'pi_g': guess}] * n_chains

    with utils.suppress_stdout_stderr():
        fit = model.sampling(
            data=data,
            pars=['alpha', 'beta'],
            n_jobs=1,
            chains=n_chains,
            warmup=n_warmup,
            iter=n_iters + n_warmup,
            init=init,
            control={'adapt_delta': 0.99},
            seed=seed,
        )
    samples = fit.extract(('alpha', 'beta'))
    alpha, beta = np.mean(samples['alpha']), np.mean(samples['beta'])
    pi = pi_func(alpha, beta)
    return guess, alpha, beta, pi


def estimate_pi(
    df_aggregates,
    p_e,
    p_c,
    pi_path,
    p_group_by=None,
    value_columns=['TC', 'T', 'count'],
    n_threads=8,
    threshold=16,
    subset_threshold=8000,
):
    df_aggregates = df_aggregates[(df_aggregates[value_columns] > 0).any(axis=1)]

    logger.debug(f'Compiling STAN model from {config.MODEL_PATH}')
    model = pystan.StanModel(file=config.MODEL_PATH, model_name=config.MODEL_NAME)

    if p_group_by is not None:
        df_full = df_aggregates.set_index(p_group_by, drop=True)
        df_full['p_e'] = df_full.index.map(p_e)
        df_full['p_c'] = df_full.index.map(p_c)
        df_full = df_full.reset_index()
    else:
        df_full = df_aggregates
        df_full['p_e'] = p_e
        df_full['p_c'] = p_c
    values = df_full[value_columns].values
    p_es = df_full['p_e'].values
    p_cs = df_full['p_c'].values
    groups = df_full.groupby(['barcode', 'GX']).indices
    pis = {}
    skipped = 0
    logger.debug(f'Spawning {n_threads} processes')
    with ProcessPoolExecutor(max_workers=n_threads, initializer=initializer, initargs=(model,)) as executor:
        futures = {}
        # Run larger groups first
        for key in sorted(groups.keys(), key=lambda key: sum(values[groups[key]][:, 2]), reverse=True):
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
                fit_stan_mcmc, values[idx], p_e, p_c, guess=guess, pi_func=beta_mode, subset_threshold=subset_threshold
            )] = key

        for future in utils.as_completed_with_progress(futures):
            key = futures[future]
            guess, alpha, beta, pi = future.result()

            pis[key] = (guess, alpha, beta, pi)
    logger.warning(
        f'Estimation skipped {skipped} times because there were less reads than '
        f'threshold ({threshold}).'
    )

    with open(pi_path, 'w') as f:
        f.write('barcode,GX,guess,alpha,beta,pi\n')
        for barcode, gx in sorted(pis.keys()):
            guess, alpha, beta, pi = pis[(barcode, gx)]
            f.write(f'{barcode},{gx},{guess},{alpha},{beta},{pi}\n')

    return pi_path


def split_matrix(matrix, pis, barcodes, features):
    pi_matrix = np.full((len(barcodes), len(features)), np.nan)
    barcode_indices = {barcode: i for i, barcode in enumerate(barcodes)}
    feature_indices = {feature: i for i, feature in enumerate(features)}

    for (barcode, gx), pi in pis.items():
        try:
            pi = float(pi)
        except ValueError:
            continue
        pi_matrix[barcode_indices[barcode], feature_indices[gx]] = pi

    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return pi_matrix, matrix * (1 - pi_matrix), matrix * pi_matrix
