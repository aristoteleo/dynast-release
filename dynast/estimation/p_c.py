import logging
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from numba import njit
from scipy import sparse

from .. import utils

logger = logging.getLogger(__name__)


def read_p_c(p_c_path, group_by=None):
    """Read p_c CSV as a dictionary, with `group_by` columns as keys.

    :param p_c_path: path to CSV containing p_c values
    :type p_c_path: str
    :param group_by: columns to group by, defaults to `None`
    :type group_by: list, optional

    :return: dictionary with `group_by` columns as keys (tuple if multiple)
    :rtype: dictionary
    """
    if group_by is None:
        with open(p_c_path, 'r') as f:
            return float(f.read())

    df = pd.read_csv(p_c_path)
    return dict(df.set_index(group_by)['p_c'])


@njit
def binomial_pmf(k, n, p):
    """Numbaized binomial PMF function for faster calculation.

    :param k: number of successes
    :type k: int
    :param n: number of trials
    :type n: int
    :param p: probability of success
    :type p: float

    :return: probability of observing `k` successes in `n` trials with probability
             of success `p`
    :rtype: float
    """
    coef = np.prod(np.arange(k)[::-1] + (n - k + 1)) / (np.prod(np.arange(k) + 1))
    return coef * (p**k) * ((1 - p)**(n - k))


def filter_aggregates_part(values, p_e, threshold=0.01):
    """Filter out conversion-content pairs that we estimate less than `threshold`
    fraction of reads is expected to orginate from unlabeled RNA.

    :param values: array of three columns encoding a sparse array in
                   (row, column, value) format, zero-indexed, where
                       row:    number of conversions
                       column: nucleotide content
                       value:  number of reads
    :type values: numpy.ndarray
    :param p_e: background mutation rate of unlabeled RNA
    :type p_e: float
    :param threshold: filter threshold, defaults to `0.01`
    :type threshold: float, optional

    :return: (sparse matrix representing `values`, filter mask)
    :rtype: (scipy.sparse.spmatrix, list)
    """
    sp = sparse.csr_matrix((values[:, 2], (values[:, 0], values[:, 1]))).tolil()
    mask = []

    for k in range(sp.shape[0]):
        for n in range(sp.shape[1]):
            expected = np.sum(sp[k + 1:, n]) * binomial_pmf(k, n, p_e)
            if expected > threshold * sp[k, n]:
                mask.append((k, n))

    return sp, mask


def filter_aggregates(df_aggregates, p_e, group_by=None, n_threads=8):
    """Wrapper around `filter_aggregates_part` function that runs in parallel.

    :param df_aggregates: Pandas dataframe containing aggregate values
    :type df_aggregates: pandas.DataFrame
    :param p_e: background mutation rate of unlabeled RNA
    :type p_e: float
    :param group_by: columns to group by, defaults to `None`
    :type group_by: list, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional

    :return: dictionary of filtered results for each group
    :rtype: dictionary
    """

    values = df_aggregates[['conversion', 'base', 'count']].values
    if group_by is None:
        return filter_aggregates_part(values, p_e)

    groups = df_aggregates.groupby(group_by).indices
    filtered = {}
    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        futures = {}
        for key in groups.keys():
            futures[executor.submit(filter_aggregates_part, values[groups[key]], p_e[key])] = key

        for future in utils.as_completed_with_progress(futures):
            key = futures[future]
            filtered[key] = future.result()

    return filtered


def expectation_maximization_part(sp, mask, threshold=1e-8):
    """Run EM algorithm to estimate average conversion rate in labeled RNA.

    :param sp: sparse matrix of three columns encoding a sparse array in
               (row, column, value) format, zero-indexed, where
                   row:    number of conversions
                   column: nucleotide content
                   value:  number of reads
    :type sp: scipy.sparse.spmatrix
    :param mask: list of tuples containing coordinates of `sp` that were filtered
    :type mask: list
    :param threshold: stop if error is less than this value, defaults to `1e-8`
    :type threshold: float, optional

    :return: (new sparse matrix with same structure as `sp`, estimated conversion rate in labeled RNA)
    :rtype: (scipy.sparse.spmatrix, float)
    """
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
    """Wrapper function around `expectation_maximization_part` that runs in parallel.

    :param filtered: output of `filter_aggregates`, containing filter results
    :type filtered: dictionary
    :param group_by: columns to group by, defaults to `None`
    :type group_by: list, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional

    :return: dictionary of EM results for each group
    :rtype: dictionary
    """
    if group_by is None:
        sp, mask = filtered
        return expectation_maximization_part(sp, mask)

    em = {}
    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        futures = {}
        for key, (sp, mask) in filtered.items():
            futures[executor.submit(expectation_maximization_part, sp, mask)] = key

        for future in utils.as_completed_with_progress(futures):
            key = futures[future]
            em[key] = future.result()

    return em


def estimate_p_c(
    df_aggregates,
    p_e,
    p_c_path,
    group_by=None,
    n_threads=8,
):
    """Estimate the average conversion rate in labeled RNA.

    :param df_aggregates: Pandas dataframe containing aggregate values
    :type df_aggregates: pandas.DataFrame
    :param p_e: background mutation rate of unlabeled RNA
    :type p_e: float
    :param p_c_path: path to output CSV containing p_c estimates
    :type p_c_path: str
    :param group_by: columns to group by, defaults to `None`
    :type group_by: list, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional

    :return: path to output CSV containing p_c estimates
    :rtype: str
    """
    logger.debug(f'p_c estimation will be grouped by {group_by}')
    logger.debug('Filtering aggregates')
    filtered = filter_aggregates(df_aggregates, p_e, group_by=group_by, n_threads=n_threads)

    logger.debug('Running EM algorithm')
    em = expectation_maximization(filtered, group_by=group_by, n_threads=n_threads)

    logger.debug(f'Writing p_c estimates to {p_c_path}')
    with open(p_c_path, 'w') as p_c_out:
        if group_by is None:
            sp, p_c = em
            p_c_out.write(str(p_c))
        else:
            p_c_out.write(f'{",".join(group_by)},p_c\n')
            for key in sorted(em.keys()):
                sp, p_c = em[key]
                formatted_key = key if isinstance(key, str) else ",".join(key)
                p_c_out.write(f'{formatted_key},{p_c}\n')

    return p_c_path
