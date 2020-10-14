import logging
from concurrent.futures import as_completed, ProcessPoolExecutor

import numpy as np
import pandas as pd
from numba import njit
from scipy import sparse
from tqdm import tqdm

logger = logging.getLogger(__name__)


def estimate_p_e(df_rates, p_e_path):
    p_e = (df_rates['CT'] + df_rates['GA']) / 2

    # Return a dictionary of index-to-pe pairs if the above is a series
    if isinstance(p_e, pd.Series):
        p_e = p_e.to_dict()
        with open(p_e_path, 'w') as f:
            for key in sorted(p_e.keys()):
                f.write(f'{key},{p_e[key]}\n')
    else:
        with open(p_e_path, 'w') as f:
            f.write(p_e)
    return p_e


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


def filter_aggregates(values, groups, p_e, n_threads=8):

    filtered = {}
    with ProcessPoolExecutor(max_workers=n_threads) as executor, tqdm(total=2 * len(groups), ascii=True) as pbar:
        futures = {}
        for key in groups.keys():
            futures[executor.submit(
                filter_aggregates_part, values[groups[key]],
                p_e.get(key, 0) if isinstance(p_e, dict) else p_e
            )] = key
            pbar.update(1)

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


def expectation_maximization(filtered, n_threads=8):
    em = {}
    with ProcessPoolExecutor(max_workers=n_threads) as executor, tqdm(total=2 * len(filtered), ascii=True) as pbar:
        futures = {}
        for key, (sp, mask) in filtered.items():
            futures[executor.submit(expectation_maximization_part, sp, mask)] = key
            pbar.update(1)

        for future in as_completed(futures):
            key = futures[future]
            em[key] = future.result()
            pbar.update(1)

    return em


def estimate_p_c(df_aggregates, p_e, p_c_path, aggregate_path, group_by='barcode', n_threads=8):
    df = df_aggregates.sort_values(group_by).reset_index()
    values = df[['TC', 'T', 'count']].values
    groups = df.groupby(group_by).indices

    logger.debug('Filtering aggregates')
    filtered = filter_aggregates(values, groups, p_e, n_threads=n_threads)

    logger.debug('Running EM algorithm')
    em = expectation_maximization(filtered, n_threads=n_threads)

    with open(p_c_path, 'w') as p_c_out, open(aggregate_path, 'w') as aggregate_out:
        for key in sorted(em.keys()):
            sp, p_c = em[key]
            p_c_out.write(f'{key},{p_c}\n')

            for (k, n), count in sp.todok().items():
                aggregate_out.write(f'{key},{k},{n},{count}\n')

    return p_c_path, aggregate_path


def estimate_pi_g():
    # TODO
    pass
