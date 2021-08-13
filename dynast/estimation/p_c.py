from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from numba import njit
from scipy import sparse

from .. import utils
from ..logging import logger


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

    df = pd.read_csv(p_c_path, dtype={key: 'string' for key in group_by})
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


def expectation_maximization_nasc(values, p_e, threshold=0.01):
    """NASC-seq pipeline variant of the EM algorithm to estimate average
    conversion rate in labeled RNA.

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

    :return: estimated conversion rate
    :rtype: float
    """
    sp = sparse.csr_matrix((values[:, 2], (values[:, 0], values[:, 1]))).tolil()
    mask = []
    for k in range(sp.shape[0]):
        for n in range(sp.shape[1]):
            expected = np.sum(sp[k + 1:, n]) * binomial_pmf(k, n, p_e)
            if expected > threshold * sp[k, n]:
                mask.append((k, n))

    r = 1
    l = 0  # noqa
    p_c = (r + l) / 2
    prev_p_c = p_c
    keys = sorted(mask)
    new_sp = sp.copy().asfptype()

    while r - l >= 1e-8:  # noqa
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

    return p_c


def expectation_maximization(values, p_e, p_c=0.1, threshold=0.01, max_iters=300):
    """Run EM algorithm to estimate average conversion rate in labeled RNA.

    This function runs the following two steps.
    1) Constructs a sparse matrix representation of `values` and filters out certain
       indices that are expected to contain more than `threshold` proportion of unlabeled reads.
    2) Runs an EM algorithm that iteratively updates the filtered out data and
       stimation.
    See https://doi.org/10.1093/bioinformatics/bty256.

    :param values: array of three columns encoding a sparse array in
                   (row, column, value) format, zero-indexed, where
                       row:    number of conversions
                       column: nucleotide content
                       value:  number of reads
    :type values: numpy.ndarray
    :param p_e: background mutation rate of unlabeled RNA
    :type p_e: float
    :param p_c: initial p_c value, defaults to `0.1`
    :type p_c: float, optional
    :param threshold: filter threshold, defaults to `0.01`
    :type threshold: float, optional
    :param max_iters: maximum number of EM iterations, defaults to `300`
    :type max_iters: int, optional

    :return: estimated conversion rate
    :rtype: float
    """
    sp = sparse.csr_matrix((values[:, 2], (values[:, 0], values[:, 1]))).tolil()
    mask = []
    for k in range(sp.shape[0]):
        for n in range(sp.shape[1]):
            expected = np.sum(sp[k:, n]) * binomial_pmf(k, n, p_e)
            if expected > threshold * sp[k, n]:
                mask.append((k, n))

    prev_p_c = p_c
    keys = sorted(mask)
    new_sp = sp.copy().asfptype()
    for _ in range(max_iters):  # noqa
        # E step
        for k, n in keys:
            numerator = 0
            denominator = 0

            for kp in range(new_sp.shape[0]):
                if (kp, n) not in mask:
                    numerator += binomial_pmf(k, n, p_c) * new_sp[kp, n]
                    denominator += binomial_pmf(kp, n, p_c)
            if denominator > 0:
                new_sp[k, n] = numerator / denominator
        # M step
        numerator = 0
        denominator = 0
        for k in range(new_sp.shape[0]):
            for n in range(new_sp.shape[1]):
                numerator += k * new_sp[k, n]
                denominator += n * new_sp[k, n]
        prev_p_c = p_c
        p_c = numerator / denominator if denominator > 0 else p_c
        if prev_p_c == p_c:
            break

    if p_c <= p_e:
        raise ValueError('p_c <= p_e')

    return p_c


def estimate_p_c(df_aggregates, p_e, p_c_path, group_by=None, threshold=1000, n_threads=8, nasc=False):
    """Estimate the average conversion rate in labeled RNA.

    :param df_aggregates: Pandas dataframe containing aggregate values
    :type df_aggregates: pandas.DataFrame
    :param p_e: background mutation rate of unlabeled RNA
    :type p_e: float
    :param p_c_path: path to output CSV containing p_c estimates
    :type p_c_path: str
    :param group_by: columns to group by, defaults to `None`
    :type group_by: list, optional
    :param threshold: read count threshold, defaults to `1000`
    :type threshold: int, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional
    :param nasc: flag to indicate whether to use NASC-seq pipeline variant of
                 the EM algorithm, defaults to `False`
    :type nasc: bool, optional

    :return: path to output CSV containing p_c estimates
    :rtype: str
    """
    em_func = expectation_maximization_nasc if nasc else expectation_maximization
    values = df_aggregates[['conversion', 'base', 'count']].values
    logger.debug('Running EM algorithm')
    if group_by is None:
        p_c = em_func(values, p_e)
    else:
        groups = df_aggregates.groupby(group_by, sort=False, observed=True).indices
        p_cs = {}
        skipped = 0
        failed = 0
        logger.debug(f'Spawning {n_threads} processes')
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
            futures = {}
            for key, idx in groups.items():
                vals = values[idx]
                if sum(vals[:, 2]) < threshold:
                    skipped += 1
                    continue
                futures[executor.submit(em_func, vals, p_e[key])] = key

            for future in utils.as_completed_with_progress(futures):
                key = futures[future]
                try:
                    p_cs[key] = future.result()
                except Exception:
                    failed += 1
        if skipped > 0:
            logger.warning(
                f'Estimation skipped for {skipped} barcodes because they contain less than '
                f'{threshold} reads. Use `--cell-threshold` to change.'
            )
        if failed > 0:
            logger.warning(f'Estimation failed {failed} times.')

    logger.debug(f'Writing p_c estimates to {p_c_path}')
    with open(p_c_path, 'w') as f:
        if group_by is None:
            f.write(str(p_c))
        else:
            f.write(f'{",".join(group_by)},p_c\n')
            for key in sorted(p_cs.keys()):
                p_c = p_cs[key]
                formatted_key = key if isinstance(key, str) else ",".join(key)
                f.write(f'{formatted_key},{p_c}\n')

    return p_c_path
