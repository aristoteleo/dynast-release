from typing import FrozenSet, List, Optional

import numpy as np
import pandas as pd

from ..logging import logger
from .conversion import BASE_IDX, COLUMNS


def read_rates(rates_path: str) -> pd.DataFrame:
    """Read mutation rates CSV as a pandas dataframe.

    Args:
        rates_path: Path to rates CSV

    Returns:
        Rates dataframe
    """
    df = pd.read_csv(rates_path, index_col=None)
    if df.shape[0] == 1:
        return df.iloc[0]
    return df


def read_aggregates(aggregates_path: str) -> pd.DataFrame:
    """Read aggregates CSV as a pandas dataframe.

    Args:
        aggregates_path: Path to aggregates CSV

    Returns:
        Aggregates dataframe
    """
    dtypes = {
        'barcode': 'category',
        'GX': 'category',
        'count': np.uint16,
        **{column: np.uint8
           for column in COLUMNS},
    }
    df = pd.read_csv(aggregates_path, dtype=dtypes)
    return df


def merge_aggregates(*dfs: pd.DataFrame) -> pd.DataFrame:
    """Merge multiple aggregate dataframes into one.

    Args:
        dfs: Dataframes to merge

    Returns:
        Merged dataframe
    """
    df = pd.concat(dfs).groupby(['barcode', 'GX', 'conversion', 'base'], sort=False, observed=True).sum().reset_index()
    df['conversion'] = df['conversion'].astype(np.uint8)
    df['base'] = df['base'].astype(np.uint8)
    df['count'] = df['count'].astype(np.uint16)
    return df


def calculate_mutation_rates(df_counts: pd.DataFrame, rates_path: str, group_by: Optional[List[str]] = None) -> str:
    """Calculate mutation rate for each pair of bases.

    Args:
        df_counts: Counts dataframe, with complemented reverse strand bases
        rates_path: Path to write rates CSV
        group_by: Column(s) to group calculations by, defaults to `None`, which
            combines all rows

    Returns:
        Path to rates CSV
    """
    logger.debug(f'Mutation rates will be grouped by {group_by}')

    if group_by is not None:
        df_sum = df_counts.groupby(group_by, sort=False, observed=True).sum(numeric_only=True).astype(np.uint32)
    else:
        df_sum = pd.DataFrame(df_counts.sum(numeric_only=True).astype(np.uint32)).T

    # Compute mutation rate of each base
    dfs = []
    for base in sorted(BASE_IDX.keys()):
        logger.debug(f'Calculating {base} mutation rate')
        dfs.append(df_sum.filter(regex=f'^{base}[A,C,G,T]$').divide(df_sum[base], axis=0))
    df_rates = pd.concat(dfs, axis=1)
    if group_by is not None:
        df_rates = df_rates.reset_index()
    logger.debug(f'Writing mutation rates to {rates_path}')
    df_rates.to_csv(rates_path, index=False)
    return rates_path


def aggregate_counts(
        df_counts: pd.DataFrame, aggregates_path: str, conversions: FrozenSet[str] = frozenset({'TC'})
) -> str:
    """Aggregate conversion counts for each pair of bases.

    Args:
        df_counts: Counts dataframe, with complemented reverse strand bases
        aggregates_path: Path to write aggregate CSV
        conversions: Conversion(s) in question

    Returns:
        Path to aggregate CSV that was written
    """
    bases = list(set(f[0] for f in conversions))
    df_combined = df_counts[['barcode', 'GX']].copy()
    df_combined['conversion'] = df_counts[list(conversions)].sum(axis=1)
    df_combined['base'] = df_counts[bases].sum(axis=1)
    pd.DataFrame(
        df_combined.groupby(
            ['barcode', 'GX', 'conversion', 'base'],
            sort=False,
            observed=True,
        ).size(),
        columns=['count']
    ).reset_index().to_csv(
        aggregates_path, index=False
    )
    return aggregates_path
