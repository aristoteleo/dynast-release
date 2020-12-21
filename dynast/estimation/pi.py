import logging
from concurrent.futures import as_completed, ProcessPoolExecutor

import numpy as np
import pandas as pd
import pystan
from tqdm import tqdm

from .. import config, utils

logger = logging.getLogger(__name__)


def read_pi(pi_path, group_by=None):
    if group_by is None:
        with open(pi_path, 'r') as f:
            return float(f.read().split(',')[-1])

    df = pd.read_csv(pi_path, usecols=group_by + ['pi'])
    return dict(df.set_index(group_by)['pi'])


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


def fit_stan_optimize(
    values,
    p_e,
    p_c,
    guess=0.5,
    model=None,
    pi_func=lambda alpha,
    beta: None,
    n_chains=4,
    n_iters=2000,
):
    model = model or _model
    conversions = []
    contents = []
    for k, n, count in values:
        conversions.extend([k] * count)
        contents.extend([n] * count)
    data = {
        'N': len(conversions),
        'contents': contents,
        'conversions': conversions,
        'p_c': p_c,
        'p_e': p_e,
    }

    # First, try optimizing.
    try:
        with utils.suppress_stdout_stderr():
            init = [{'log_alpha': 0., 'log_beta': 0., 'pi_g': guess}]
            fit = model.optimizing(data=data, init=init, algorithm='BFGS')
        alpha, beta = fit['alpha'], fit['beta']
        pi = pi_func(alpha, beta)
        if pi is not None:
            return alpha, beta, pi, False
    except Exception:
        pass

    # Fallback to MCMC
    with utils.suppress_stdout_stderr():
        init = [{'log_alpha': 0., 'log_beta': 0., 'pi_g': guess}] * n_chains
        fit = model.sampling(
            data=data,
            n_jobs=1,
            iter=n_iters,
            chains=n_chains,
            init=init,
            control={'adapt_delta': 0.99},
        )
    samples = fit.extract(('alpha', 'beta'))
    alpha, beta = np.mean(samples['alpha']), np.mean(samples['beta'])
    pi = pi_func(alpha, beta)
    return guess, alpha, beta, pi, True


def fit_stan_mcmc(
    values,
    p_e,
    p_c,
    guess=0.5,
    model=None,
    pi_func=(lambda alpha, beta: None),
    n_chains=1,
    n_iters=2000,
):
    model = model or _model
    conversions = []
    contents = []
    for k, n, count in values:
        conversions.extend([k] * count)
        contents.extend([n] * count)

    data = {
        'N': len(conversions),
        'contents': contents,
        'conversions': conversions,
        'p_c': p_c,
        'p_e': p_e,
    }
    init = [{'log_alpha': 0.0, 'log_beta': 0.0, 'pi_g': guess}] * n_chains

    with utils.suppress_stdout_stderr():
        fit = model.sampling(
            data=data,
            pars=['alpha', 'beta'],
            n_jobs=1,
            iter=n_iters,
            chains=n_chains,
            init=init,
            control={'adapt_delta': 0.99},
        )
    samples = fit.extract(('alpha', 'beta'))
    alpha, beta = np.mean(samples['alpha']), np.mean(samples['beta'])
    pi = pi_func(alpha, beta)
    return guess, alpha, beta, pi, False


def estimate_pi(
    df_aggregates,
    p_e,
    p_c,
    pi_path,
    filter_dict=None,
    p_group_by=None,
    group_by=None,
    value_columns=['TC', 'T', 'count'],
    n_threads=8,
    nasc=False,
    method='mcmc',
    threshold=16,
):
    pi_func = beta_mode
    if method == 'optimization':
        fit = fit_stan_optimize
    else:
        fit = fit_stan_mcmc

    logger.debug(f'pi estimation will be grouped by {group_by} using columns {value_columns}')
    df_aggregates = df_aggregates[(df_aggregates[value_columns] > 0).any(axis=1)]
    filter_dict = filter_dict or {}
    if filter_dict:
        logger.debug(f'Filtering aggregates by the following keys: {list(filter_dict.keys())}')
        for column, values in filter_dict.items():
            df_aggregates = df_aggregates[df_aggregates[column].isin(values)]

    logger.debug(f'Compiling STAN model from {config.MODEL_PATH}')
    model = pystan.StanModel(file=config.MODEL_PATH, model_name=config.MODEL_NAME)
    if group_by is None:
        if p_group_by is not None:
            raise Exception('Can not group all aggregates when p_e and p_c are not constants')
        guess, alpha, beta, pi, fallback = fit(
            df_aggregates[value_columns].values, p_e, p_c, model=model, pi_func=pi_func
        )
        if method == 'optimization' and fallback:
            logger.warning('STAN optimization failed. Fallback to MCMC.')

        with open(pi_path, 'w') as f:
            f.write('guess,alpha,beta,pi\n')
            f.write(f'{guess},{alpha},{beta},{pi}\n')
        return pi_path

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
    groups = df_full.groupby(group_by).indices
    pis = {}
    fallbacks = 0
    skipped = 0
    logger.debug(f'Spawning {n_threads} processes')
    with ProcessPoolExecutor(max_workers=n_threads, initializer=initializer,
                             initargs=(model,)) as executor, tqdm(total=2 * len(groups), ascii=True) as pbar:
        pbar.set_description('Queueing')
        futures = {}
        for key, idx in list(groups.items()):
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

            if count < 16:
                skipped += 1
                pbar.update(2)
                continue

            # Make a naive guess of the fraction of new RNA
            # Clip guess to [0.01, 0.99] because STAN initialization will fail
            # if guess is either 0 or 1.
            guess = min(max(sum(vals[vals[:, 0] > 0][:, 2]) / count, 0.01), 0.99)

            futures[executor.submit(fit, values[idx], p_e, p_c, guess=guess, pi_func=pi_func)] = key
            pbar.update(1)

        pbar.set_description('Executing')
        for future in as_completed(futures):
            key = futures[future]
            guess, alpha, beta, pi, fallback = future.result()
            if fallback:
                fallbacks += 1

            pis[key] = (guess, alpha, beta, pi)
            pbar.update(1)
    logger.warning(
        f'Estimation skipped {skipped} times because there were less reads than '
        f'threshold ({threshold}).'
    )
    if method == 'optimization':
        logger.warning(f'STAN optimization failed and fallback to MCMC {fallbacks} times.')

    with open(pi_path, 'w') as f:
        f.write(f'{",".join(group_by)},guess,alpha,beta,pi\n')
        for key in sorted(pis.keys()):
            guess, alpha, beta, pi = pis[key]
            f.write(f'{key if isinstance(key, str) else ",".join(key)},{guess},{alpha},{beta},{pi}\n')

    return pi_path


def split_reads(adata, pis, group_by=None):
    barcodes = adata.obs.barcode
    gene_ids = adata.var.gene_id

    pi_matrix = np.full(adata.shape, np.nan)
    barcode_index = {barcode: i for i, barcode in enumerate(barcodes)}
    gene_id_index = {gene_id: i for i, gene_id in enumerate(gene_ids)}

    for key, pi in pis.items():
        try:
            pi = float(pi)
        except ValueError:
            continue

        if group_by == ['barcode', 'GX']:
            barcode, gx = key
            pi_matrix[barcode_index[barcode], gene_id_index[gx]] = pi
        elif group_by == ['barcode']:
            barcode = key
            pi_matrix[barcode_index[barcode]] = pi
        elif group_by == ['GX']:
            gx = key
            pi_matrix[:, gene_id_index[gx]] = pi
        elif group_by is None:
            pi_matrix[:, :] = pi
        else:
            raise Exception(f'Unknown group_by {group_by}')

    adata.layers['pi'] = pi_matrix
    adata.layers['new_estimated'] = adata.X * pi_matrix
    adata.layers['old_estimated'] = adata.X * (1 - pi_matrix)
    return adata
