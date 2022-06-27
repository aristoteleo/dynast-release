from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from numba import njit
from scipy import sparse

from .. import utils
from ..logging import logger


def read_p_c(p_c_path: str,
             group_by: Optional[List[str]] = None) -> Union[float, Dict[str, float], Dict[Tuple[str, ...], float]]:
    """Read p_c CSV as a dictionary, with `group_by` columns as keys.

    Args:
        p_c_path: Path to CSV containing p_c values
        group_by: Columns to group by, defaults to `None`

    Returns:
        Dictionary with `group_by` columns as keys (tuple if multiple)
    """
    if group_by is None:
        with open(p_c_path, 'r') as f:
            return float(f.read())

    df = pd.read_csv(p_c_path, dtype={key: 'category' for key in group_by})
    return dict(df.set_index(group_by)['p_c'])


@njit(cache=True)
def binomial_pmf(k: int, n: int, p: int) -> float:
    """Numbaized binomial PMF function for faster calculation.

    Args:
        k: Number of successes
        n: Number of trials
        p: Probability of success

    Returns:
        Probability of observing `k` successes in `n` trials with probability
            of success `p`
    """
    coef = np.prod(np.arange(k)[::-1] + (n - k + 1)) / (np.prod(np.arange(k) + 1))
    return coef * (p**k) * ((1 - p)**(n - k))


@njit(cache=True)
def expectation_maximization_nasc(values: np.ndarray, p_e: float, threshold: float = 0.01) -> float:
    """NASC-seq pipeline variant of the EM algorithm to estimate average
    conversion rate in labeled RNA.

    Args:
        values: N x C Numpy array where N is the number of conversions, C is the nucleotide
            content, and the value at this position is the number of reads observed
        p_e: Background mutation rate of unlabeled RNA
        threshold: Filter threshold

    Returns:
        Estimated conversion rate
    """
    mask = []
    for k in range(values.shape[0]):
        for n in range(values.shape[1]):
            expected = (values[k + 1:, n]).sum() * binomial_pmf(k, n, p_e)
            if expected > threshold * values[k, n]:
                mask.append((k, n))

    r = 1
    l = 0  # noqa
    p_c = (r + l) / 2
    prev_p_c = p_c
    keys = sorted(mask)
    keys_set = set(mask)
    new_values = values.astype(np.float32)  # This makes a copy

    while r - l >= 1e-8:  # noqa
        # E step
        for k, n in keys:
            numerator = 0
            denominator = 0

            for kp in range(new_values.shape[0]):
                if (kp, n) in keys_set:
                    numerator += binomial_pmf(k, n, p_c) * new_values[kp, n]
                    denominator += binomial_pmf(kp, n, p_c)
            new_values[k, n] = numerator / denominator

        # M step
        numerator = 0
        denominator = 0
        for k in range(new_values.shape[0]):
            for n in range(new_values.shape[1]):
                numerator += k * new_values[k, n]
                denominator += n * new_values[k, n]
        prev_p_c = p_c
        p_c = numerator / denominator if denominator > 0 else 0

        if prev_p_c == p_c:
            break

        if p_c < prev_p_c:
            r = prev_p_c
        else:
            l = prev_p_c  # noqa

    return p_c


@njit(cache=True)
def expectation_maximization(
    values: np.ndarray, p_e: float, p_c: float = 0.1, threshold: float = 0.01, max_iters: int = 300
) -> float:
    """Run EM algorithm to estimate average conversion rate in labeled RNA.

    This function runs the following two steps.
    1) Constructs a sparse matrix representation of `values` and filters out certain
       indices that are expected to contain more than `threshold` proportion of unlabeled reads.
    2) Runs an EM algorithm that iteratively updates the filtered out data and
       stimation.
    See https://doi.org/10.1093/bioinformatics/bty256.

    Args:
        values: array of three columns encoding a sparse array in (row, column, value) format,
            zero-indexed, where
            row:    number of conversions
            column: nucleotide content
            value:  number of reads
        p_e: Background mutation rate of unlabeled RNA
        p_c: Initial p_c value
        threshold: Filter threshold
        max_iters: Maximum number of EM iterations

    Returns:
        Estimated conversion rate
    """
    mask = []
    for k in range(values.shape[0]):
        for n in range(values.shape[1]):
            expected = (values[k:, n]).sum() * binomial_pmf(k, n, p_e)
            if expected > threshold * values[k, n]:
                mask.append((k, n))

    prev_p_c = p_c
    keys = sorted(mask)
    keys_set = set(mask)
    new_values = values.astype(np.float32)  # This makes a copy
    for _ in range(max_iters):  # noqa
        # E step
        for k, n in keys:
            numerator = 0
            denominator = 0

            for kp in range(new_values.shape[0]):
                if (kp, n) not in keys_set:
                    numerator += binomial_pmf(k, n, p_c) * new_values[kp, n]
                    denominator += binomial_pmf(kp, n, p_c)
            if denominator > 0:
                new_values[k, n] = numerator / denominator
        # M step
        numerator = (new_values * np.arange(new_values.shape[0]).reshape(-1, 1)).sum()
        denominator = (new_values * np.arange(new_values.shape[1])).sum()
        prev_p_c = p_c
        p_c = numerator / denominator if denominator > 0 else p_c
        if prev_p_c == p_c:
            break

    if p_c <= p_e:
        raise ValueError('p_c <= p_e')

    return p_c


def estimate_p_c(
    df_aggregates: pd.DataFrame,
    p_e: Union[float, Dict[str, float], Dict[Tuple[str, ...], float]],
    p_c_path: str,
    group_by: Optional[List[str]] = None,
    threshold: int = 1000,
    n_threads: int = 8,
    nasc: bool = False
) -> str:
    """Estimate the average conversion rate in labeled RNA.

    Args:
        df_aggregates: Pandas dataframe containing aggregate values
        p_e: Background mutation rate of unlabeled RNA
        p_c_path: Path to output CSV containing p_c estimates
        group_by: Columns to group by
        threshold: Read count threshold
        n_threads: Number of threads
        nasc: Flag to indicate whether to use NASC-seq pipeline variant of the EM algorithm

    Returns:
        Path to output CSV containing p_c estimates
    """
    em_func = expectation_maximization_nasc if nasc else expectation_maximization
    values = df_aggregates[['conversion', 'base', 'count']].values
    logger.debug('Running EM algorithm')
    if group_by is None:
        # Initialize with sparse matrix because duplicates need to be added.
        values_mtx = sparse.csr_matrix((values[:, 2], (values[:, 0], values[:, 1]))).A
        p_c = em_func(values_mtx, p_e)
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
                # Initialize with sparse matrix because duplicates need to be added.
                vals_mtx = sparse.csr_matrix((vals[:, 2], (vals[:, 0], vals[:, 1]))).A
                futures[executor.submit(em_func, vals_mtx, p_e[key])] = key

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
