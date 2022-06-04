import numpy as np
import pandas as pd

from .. import utils
from ..logging import logger
from .conversion import BASE_IDX, COLUMNS


def read_rates(rates_path):
    """Read mutation rates CSV as a pandas dataframe.

    :param rates_path: path to rates CSV
    :type rates_path: str

    :return: rates dataframe
    :rtype: pandas.DataFrame
    """
    df = pd.read_csv(rates_path, index_col=None)
    if df.shape[0] == 1:
        return df.iloc[0]
    return df


def read_aggregates(aggregates_path):
    """Read aggregates CSV as a pandas dataframe.

    :param aggregates_path: path to aggregates CSV
    :type aggregates_path: str

    :return: aggregates dataframe
    :rtype: pandas.DataFrame
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


def merge_aggregates(*dfs):
    """Merge multiple aggregate dataframes into one.

    :param *dfs: dataframes to merge
    :type *dfs: pandas.DataFrame

    :return: merged dataframe
    :rtype: pandas.DataFrame
    """
    df = pd.concat(dfs).groupby(['barcode', 'GX', 'conversion', 'base'], sort=False, observed=True).sum().reset_index()
    df['conversion'] = df['conversion'].astype(np.uint8)
    df['base'] = df['base'].astype(np.uint8)
    df['count'] = df['count'].astype(np.uint16)
    return df


def calculate_mutation_rates(df_counts, rates_path, group_by=None):
    """Calculate mutation rate for each pair of bases.

    :param df_counts: counts dataframe, with complemented reverse strand bases
    :type df_counts: pandas.DataFrame
    :param rates_path: path to write rates CSV
    :type rates_path: str
    :param group_by: column(s) to group calculations by, defaults to `None`, which
                     combines all rows
    :type group_by: list

    :return: path to rates CSV
    :rtype: str
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


def aggregate_counts(df_counts, aggregates_path, conversions=frozenset([('TC',)])):
    """Aggregate conversion counts for each pair of bases.

    :param df_counts: counts dataframe, with complemented reverse strand bases
    :type df_counts: pandas.DataFrame
    :param aggregates_path: path to write aggregate CSV
    :type aggregates_path: str
    :param conversions: conversion(s) in question, defaults to `frozenset([('TC',)])`
    :type conversions: list, optional

    :return: path to aggregate CSV that was written
    :rtype: str
    """
    flattened = list(utils.flatten_iter(conversions))
    bases = list(set(f[0] for f in flattened))
    df_combined = df_counts[['barcode', 'GX']].copy()
    df_combined['conversion'] = df_counts[flattened].sum(axis=1)
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
