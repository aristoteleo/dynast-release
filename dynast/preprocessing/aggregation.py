import logging
import os

import numpy as np
import pandas as pd
from tqdm import tqdm

from .conversion import BASE_COLUMNS, BASE_IDX, COLUMNS, CONVERSION_COLUMNS

logger = logging.getLogger(__name__)


def read_rates(rates_path):
    return pd.read_csv(rates_path, index_col=None)


def read_aggregates(aggregates_path):
    dtypes = {'barcode': 'string', 'GX': 'string', **{column: np.uint16 for column in COLUMNS}}
    return pd.read_csv(aggregates_path, dtype=dtypes)


def calculate_mutation_rates(df_counts, rates_path, group_by=None):
    logger.debug(f'Mutation rates will be grouped by {group_by}')
    if group_by is not None:
        df_sum = df_counts.groupby(group_by).sum(numeric_only=True).astype(np.uint32)
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


def aggregate_counts(df, df_genes, aggregates_dir):
    df = df.merge(df_genes[['GX', 'strand']], on='GX')
    df_forward = df[df.strand == '+']
    df_reverse = df[df.strand == '-'][COLUMNS]

    # Rename columns of reverse
    df_reverse.columns = list(reversed(CONVERSION_COLUMNS)) + list(reversed(BASE_COLUMNS))

    df_complemented = pd.concat((df_forward, df_reverse[COLUMNS]))

    paths = {}
    for conversion in tqdm(CONVERSION_COLUMNS, ascii=True):
        csv_path = os.path.join(aggregates_dir, f'{conversion}.csv')
        logger.debug(f'Aggregating counts for {conversion} conversion to {csv_path}')
        df_agg = pd.DataFrame(df_complemented.groupby(['barcode', 'GX', conversion, conversion[0]]).size())
        df_agg.columns = ['count']
        df_agg.reset_index().to_csv(csv_path, index=False)
        paths[conversion] = csv_path

    return paths
