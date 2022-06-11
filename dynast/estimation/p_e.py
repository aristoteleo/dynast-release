from typing import Dict, FrozenSet, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .. import utils
from ..logging import logger
from ..preprocessing.conversion import BASE_COLUMNS, CONVERSION_COLUMNS


def read_p_e(p_e_path: str, group_by: Optional[List[str]] = None) -> Dict[Union[str, Tuple[str, ...]], float]:
    """Read p_e CSV as a dictionary, with `group_by` columns as keys.

    Args:
        p_e_path: Path to CSV containing p_e values
        group_by: Columns to group by

    Returns:
        Dictionary with `group_by` columns as keys (tuple if multiple)
    """
    if group_by is None:
        with open(p_e_path, 'r') as f:
            return float(f.read())

    df = pd.read_csv(p_e_path, dtype={key: 'string' for key in group_by})
    return dict(df.set_index(group_by)['p_e'])


def estimate_p_e_control(
    df_counts: pd.DataFrame, p_e_path: str, conversions: FrozenSet[FrozenSet[str]] = frozenset({frozenset({'TC'})})
) -> str:
    """Estimate background mutation rate of unlabeled RNA for a control sample
    by simply calculating the average mutation rate.

    Args:
        df_counts: Pandas dataframe containing number of each conversion and
            nucleotide content of each read
        p_e_path: Path to output CSV containing p_e estimates
        conversions: Conversion(s) in question

    Returns:
        Path to output CSV containing p_e estimates
    """
    flattened = list(utils.flatten_iter(conversions))
    bases = list(set(f[0] for f in flattened))
    p_e = df_counts[flattened].sum().sum() / df_counts[bases].sum().sum()
    with open(p_e_path, 'w') as f:
        f.write(str(p_e))
    return p_e_path


def estimate_p_e(
    df_counts: pd.DataFrame,
    p_e_path: str,
    conversions: FrozenSet[FrozenSet[str]] = frozenset({frozenset({'TC'})}),
    group_by: Optional[List[str]] = None
) -> str:
    """Estimate background mutation rate of unabeled RNA by calculating the
    average mutation rate of all three nucleotides other than `conversion[0]`.

    Args:
        df_counts: Pandas dataframe containing number of each conversion and
            nucleotide content of each read
        p_e_path: Path to output CSV containing p_e estimates
        conversions: Conversion(s) in question, defaults to `frozenset([('TC',)])`
        group_by: Columns to group by, defaults to `None`

    Returns:
        Path to output CSV containing p_e estimates
    """
    flattened = list(utils.flatten_iter(conversions))
    bases = sorted(set(f[0] for f in flattened))
    if group_by is not None:
        df_sum = df_counts.groupby(group_by, sort=False, observed=True).sum(numeric_only=True).astype(np.uint32)
    else:
        df_sum = pd.DataFrame(df_counts.sum(numeric_only=True).astype(np.uint32)).T

    # It's best to use conversions that don't start with a conversion base.
    # For example, if the conversion is TC, don't use conversions starting with a T.
    # However, if multiple conversions are provided, and they span all bases,
    # we have no choice but to use them.
    conversion_columns = [conv for conv in CONVERSION_COLUMNS if conv[0] not in bases]
    if bases == BASE_COLUMNS:
        logger.warning(
            'All four bases have conversions, so background estimation will fall back to '
            'using ALL non-induced conversions. This may lead to an underestimate. '
            'Please consider using a control sample with `--p-e`.'
        )
        conversion_columns = [conv for conv in CONVERSION_COLUMNS if conv not in flattened]

    for conversion in conversion_columns:
        df_sum[conversion] /= df_sum[conversion[0]]
    p_e = df_sum[conversion_columns].mean(axis=1)

    # # Filter for columns not starting with the conversion base.
    # # If conversion='TC', then select columns that don't start with 'T'
    # base_columns = [base for base in BASE_COLUMNS if base not in bases]
    # conversion_columns = [conv for conv in CONVERSION_COLUMNS if conv[0] not in bases]
    # p_e = df_sum[conversion_columns].sum(axis=1) / df_sum[base_columns].sum(axis=1)

    if group_by is not None:
        p_e.reset_index().to_csv(p_e_path, header=group_by + ['p_e'], index=False)
    else:
        p_e = p_e[0]
        with open(p_e_path, 'w') as f:
            f.write(str(p_e))

    return p_e_path


def estimate_p_e_nasc(df_rates: pd.DataFrame, p_e_path: str, group_by: Optional[List[str]] = None) -> str:
    """Estimate background mutation rate of unabeled RNA by calculating the
    average `CT` and `GA` mutation rates. This function imitates the procedure
    implemented in the NASC-seq pipeline (DOI: 10.1038/s41467-019-11028-9).

    Args:
        df_counts: Pandas dataframe containing number of each conversion and
            nucleotide content of each read
        p_e_path: Path to output CSV containing p_e estimates
        group_by: Columns to group by, defaults to `None`

    Returns:
        Path to output CSV containing p_e estimates
    """
    if group_by is not None:
        df_rates = df_rates.set_index(group_by)
    p_e = (df_rates['CT'] + df_rates['GA']) / 2
    if group_by is not None:
        p_e.reset_index().to_csv(p_e_path, header=group_by + ['p_e'], index=False)
    else:
        with open(p_e_path, 'w') as f:
            f.write(str(p_e))

    return p_e_path
