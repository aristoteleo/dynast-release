import logging
from concurrent.futures import as_completed, ProcessPoolExecutor

import numpy as np
import pandas as pd
import pystan
from numba import njit
from scipy import sparse
from tqdm import tqdm

from . import config, utils

logger = logging.getLogger(__name__)


def read_p_e(p_e_path, group_by=None):
    if group_by is None:
        with open(p_e_path, 'r') as f:
            return float(f.read())

    df = pd.read_csv(p_e_path)
    return dict(df.set_index(group_by)['p_e'])


def estimate_p_e(df_rates, p_e_path, group_by=None):
    logger.debug(f'p_e estimation will be grouped by {group_by} to {p_e_path}')
    if group_by is not None:
        df_rates = df_rates.set_index(group_by)
    p_e = (df_rates['CT'] + df_rates['GA']) / 2
    logger.debug(f'Writing p_e estimates to {p_e_path}')
    if group_by is not None:
        p_e.reset_index().to_csv(p_e_path, header=group_by + ['p_e'], index=False)
    else:
        p_e = p_e[0]
        with open(p_e_path, 'w') as f:
            f.write(str(p_e))

    return p_e, p_e_path


@njit
def binomial_pmf(k, n, p):
    coef = np.prod(np.arange(k)[::-1] + (n - k + 1)) / (np.prod(np.arange(k) + 1))
    return coef * (p**k) * ((1 - p)**(n - k))


# def filter_aggregates_part(values, p_e, threshold=0.01):
#     # k = values[:,0]
#     # n = values[:,1]
#     # count = values[:,2]
#     mask = np.zeros(values.shape[0], dtype=bool)
#     for i, (k, n, count) in enumerate(values):
#         left = np.sum(values[(values[:, 0] > k) & (values[:, 1] == n), 2]) * binomial_pmf(k, n, p_e)
#         right = threshold * count
#         if left > right:
#             mask[i] = True
#     return mask


def filter_aggregates_part(values, p_e, threshold=0.01):
    sp = sparse.csr_matrix((values[:, 2], (values[:, 0], values[:, 1]))).tolil()
    mask = []

    for k in range(sp.shape[0]):
        for n in range(sp.shape[1]):
            expected = np.sum(sp[k + 1:, n]) * binomial_pmf(k, n, p_e)
            if expected > 0.01 * sp[k, n]:
                mask.append((k, n))

    return sp, mask


def filter_aggregates(df_aggregates, p_e, group_by=None, value_columns=['TC', 'T', 'count'], n_threads=8):
    values = df_aggregates[value_columns].values
    if group_by is None:
        return filter_aggregates_part(values, p_e)

    groups = df_aggregates.groupby(group_by).indices
    filtered = {}
    with ProcessPoolExecutor(max_workers=n_threads) as executor, tqdm(total=2 * len(groups), ascii=True) as pbar:
        pbar.set_description('Queueing')
        futures = {}
        for key in groups.keys():
            futures[executor.submit(filter_aggregates_part, values[groups[key]], p_e[key])] = key
            pbar.update(1)

        pbar.set_description('Executing')
        for future in as_completed(futures):
            key = futures[future]
            filtered[key] = future.result()
            pbar.update(1)

    return filtered


def expectation_maximization_part(sp, mask, threshold=10e-8):
    r = 1
    l = 0  # noqa
    p_c = (r + l) / 2
    prev_p_c = p_c
    keys = sorted(mask)
    new_sp = sp.copy().asfptype()

    while r - l >= threshold:  # noqa
        # E step
        for k, n in keys:
            numerator = 0
            denominator = 0

            for kp in range(new_sp.shape[0]):
                if (kp, n) in mask:
                    numerator += binomial_pmf(k, n, p_c) * new_sp[kp, n]
                    denominator += binomial_pmf(kp, n, p_c)
            new_sp[k, n] = numerator / denominator

        # M step
        numerator = 0
        denominator = 0
        for k in range(new_sp.shape[0]):
            for n in range(new_sp.shape[1]):
                numerator += k * new_sp[k, n]
                denominator += n * new_sp[k, n]
        prev_p_c = p_c
        p_c = numerator / denominator if denominator > 0 else 0

        if prev_p_c == p_c:
            break

        if p_c < prev_p_c:
            r = p_c
        else:
            l = p_c  # noqa

    return new_sp, p_c


def expectation_maximization(filtered, group_by=None, n_threads=8):
    if group_by is None:
        sp, mask = filtered
        return expectation_maximization_part(sp, mask)

    em = {}
    with ProcessPoolExecutor(max_workers=n_threads) as executor, tqdm(total=2 * len(filtered), ascii=True) as pbar:
        pbar.set_description('Queueing')
        futures = {}
        for key, (sp, mask) in filtered.items():
            futures[executor.submit(expectation_maximization_part, sp, mask)] = key
            pbar.update(1)

        pbar.set_description('Executing')
        for future in as_completed(futures):
            key = futures[future]
            em[key] = future.result()
            pbar.update(1)

    return em


def read_p_c(p_c_path, group_by=None):
    if group_by is None:
        with open(p_c_path, 'r') as f:
            return float(f.read())

    df = pd.read_csv(p_c_path)
    return dict(df.set_index(group_by)['p_c'])


def estimate_p_c(
    df_aggregates,
    p_e,
    p_c_path,
    aggregate_path,
    group_by=None,
    value_columns=['TC', 'T', 'count'],
    n_threads=8,
):
    logger.debug(f'p_c estimation will be grouped by {group_by} using columns {value_columns}')
    logger.debug('Filtering aggregates')
    filtered = filter_aggregates(
        df_aggregates, p_e, group_by=group_by, value_columns=value_columns, n_threads=n_threads
    )

    logger.debug('Running EM algorithm')
    em = expectation_maximization(filtered, group_by=group_by, n_threads=n_threads)

    logger.debug(f'Writing p_c estimates to {p_c_path} and aggregates to {aggregate_path}')
    with open(p_c_path, 'w') as p_c_out, open(aggregate_path, 'w') as aggregate_out:
        if group_by is None:
            aggregate_out.write(f'{",".join(value_columns)}\n')

            sp, p_c = em
            p_c_out.write(str(p_c))

            for (k, n), count in sp.todok().items():
                aggregate_out.write(f'{k},{n},{count}\n')
        else:
            p_c_out.write(f'{",".join(group_by)},p_c\n')
            aggregate_out.write(f'{",".join(group_by)},{",".join(value_columns)}\n')
            for key in sorted(em.keys()):
                sp, p_c = em[key]
                formatted_key = key if isinstance(key, str) else ",".join(key)
                p_c_out.write(f'{formatted_key},{p_c}\n')

                for (k, n), count in sp.todok().items():
                    aggregate_out.write(f'{formatted_key},{k},{n},{count}\n')
    p_c = {key: pair[1] for key, pair in em.items()} if group_by is not None else em[1]

    return p_c, p_c_path, aggregate_path


def fit_stan(
    values,
    p_e,
    p_c,
    model,
    n_chains=5,
    n_iters=2000,
):
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
    init = [{'log_alpha': 0.0, 'log_beta': 0.0, 'pi_g': 0.5}] * n_chains

    # model = stan.load_model(pickle_path)
    with utils.suppress_stdout_stderr():
        fit = model.sampling(
            data=data, n_jobs=n_chains, iter=n_iters, chains=n_chains, init=init, control={'adapt_delta': 0.95}
        )
    samples = fit.extract(('alpha', 'beta'))
    return np.mean(samples['alpha']), np.mean(samples['beta'])


def fit_stan_group(
    values_group,
    p_e,
    p_c,
    model,
    n_chains=4,
    n_iters=5000,
):
    params = {}
    with utils.suppress_stdout_stderr():
        for key, values in values_group.items():
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
            init = [{'log_alpha': 0.0, 'log_beta': 0.0, 'pi_g': 0.5}] * n_chains

            # model = stan.load_model(pickle_path)
            fit = model.sampling(
                data=data, n_jobs=1, iter=n_iters, chains=n_chains, init=init, control={'adapt_delta': 0.95}
            )
            samples = fit.extract(('alpha', 'beta'))
            params[key] = (np.mean(samples['alpha']), np.mean(samples['beta']))
        return params


def read_pi(pi_path, group_by=None):
    if group_by is None:
        with open(pi_path, 'r') as f:
            return float(f.read().split(',')[-1])

    df = pd.read_csv(pi_path, usecols=group_by + ['pi'])
    return dict(df.set_index(group_by)['pi'])


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
):
    logger.debug(f'pi estimation will be grouped by {group_by} using columns {value_columns}')
    df_aggregates = df_aggregates[(df_aggregates[value_columns] > 0).all(axis=1)]
    filter_dict = filter_dict or {}
    if filter_dict:
        logger.debug(f'Filtering aggregates by the following keys: {list(filter_dict.keys())}')
        for column, values in filter_dict.items():
            df_aggregates = df_aggregates[df_aggregates[column].isin(values)]

    logger.debug('Summing aggregates by group')
    df_aggregates = df_aggregates.groupby((group_by if group_by is not None else []) +
                                          value_columns[:-1]).sum().reset_index()

    logger.debug(f'Loading STAN model from {config.MODEL_PATH}')
    model = pystan.StanModel(file=config.MODEL_PATH, model_name=config.MODEL_NAME)
    if group_by is None:
        if p_group_by is not None:
            raise Exception('Can not group all aggregates when p_e and p_c are not constants')
        alpha, beta = fit_stan(df_aggregates[value_columns].values, p_e, p_c, model)
        pi = None
        # pi = mode of Beta(alpha, beta)
        if alpha > 1 and beta > 1:
            pi = (alpha - 1) / (alpha + beta - 2)
        elif alpha > 1 and beta <= 1:
            pi = 1
        elif alpha <= 1 and beta > 1:
            pi = 0

        with open(pi_path, 'w') as f:
            f.write(f'{alpha},{beta},{pi}')
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
    with ProcessPoolExecutor(max_workers=n_threads) as executor, tqdm(total=2 * len(groups), ascii=True) as pbar:
        pbar.set_description('Queueing')
        futures = {}
        for key, idx in groups.items():
            p_e_unique = np.unique(p_es[idx])
            p_c_unique = np.unique(p_cs[idx])

            if len(p_e_unique) > 1 or len(p_c_unique) > 1:
                raise Exception(
                    'p_e and p_c for each aggregate group must be a constant, '
                    f'but instead got {p_e_unique} and {p_c_unique}'
                )

            p_e = p_e_unique[0]
            p_c = p_c_unique[0]
            futures[executor.submit(fit_stan, values[idx], p_e, p_c, model)] = key
            pbar.update(1)

        pbar.set_description('Executing')
        for future in as_completed(futures):
            key = futures[future]
            alpha, beta = future.result()
            pi = None

            # pi = mode of Beta(alpha, beta)
            if alpha > 1 and beta > 1:
                pi = (alpha - 1) / (alpha + beta - 2)
            elif alpha > 1 and beta <= 1:
                pi = 1
            elif alpha <= 1 and beta > 1:
                pi = 0

            pis[key] = (alpha, beta, pi)
            pbar.update(1)

    with open(pi_path, 'w') as f:
        f.write(f'{",".join(group_by)},alpha,beta,pi\n')
        for key in sorted(pis.keys()):
            alpha, beta, pi = pis[key]
            f.write(f'{key if isinstance(key, str) else ",".join(key)},{alpha},{beta},{pi}\n')

    return pi_path


def split_reads(adata, pis, group_by=None):
    barcodes = adata.obs.barcode
    gene_ids = adata.var.gene_id

    pi_matrix = np.full(adata.shape, np.nan)
    barcode_index = {barcode: i for i, barcode in enumerate(barcodes)}
    gene_id_index = {gene_id: i for i, gene_id in enumerate(gene_ids)}

    for key, pi in pis.items():
        if group_by == ['barcode', 'GX']:
            barcode, gx = key
            pi_matrix[barcode_index[barcode], gene_id_index[gx]] = pi

    adata.layers['pi'] = pi_matrix
    adata.layers['new'] = adata.X * pi_matrix
    adata.layers['old'] = adata.X * (1 - pi_matrix)
    return adata
