import logging
from concurrent.futures import as_completed, ProcessPoolExecutor

import numpy as np
import pandas as pd
from numba import njit
from scipy import sparse
from tqdm import tqdm

logger = logging.getLogger(__name__)


def read_p_c(p_c_path, group_by=None):
    if group_by is None:
        with open(p_c_path, 'r') as f:
            return float(f.read())

    df = pd.read_csv(p_c_path)
    return dict(df.set_index(group_by)['p_c'])


@njit
def binomial_pmf(k, n, p):
    coef = np.prod(np.arange(k)[::-1] + (n - k + 1)) / (np.prod(np.arange(k) + 1))
    return coef * (p**k) * ((1 - p)**(n - k))


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


def expectation_maximization_part(sp, mask, threshold=1e-8):
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
            r = prev_p_c
        else:
            l = prev_p_c  # noqa

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
