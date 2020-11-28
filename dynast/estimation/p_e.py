import logging

import numpy as np
import pandas as pd

from ..preprocessing.conversion import BASE_COLUMNS, CONVERSION_COLUMNS

logger = logging.getLogger(__name__)


def read_p_e(p_e_path, group_by=None):
    if group_by is None:
        with open(p_e_path, 'r') as f:
            return float(f.read())

    df = pd.read_csv(p_e_path)
    return dict(df.set_index(group_by)['p_e'])


def estimate_p_e(df_counts, p_e_path, conversion='TC', group_by=None):
    logger.debug(f'p_e estimation will be grouped by {group_by} to {p_e_path}')
    if group_by is not None:
        df_sum = df_counts.groupby(group_by).sum(numeric_only=True).astype(np.uint32)
    else:
        df_sum = pd.DataFrame(df_counts.sum(numeric_only=True).astype(np.uint32)).T

    # Filter for columns not starting with the conversion base.
    # If conversion='TC', then select columns that don't start with 'T'
    base_columns = [base for base in BASE_COLUMNS if base != conversion[0]]
    conversion_columns = [conv for conv in CONVERSION_COLUMNS if conv[0] != conversion[0]]
    p_e = df_sum[conversion_columns].sum(axis=1) / df_sum[base_columns].sum(axis=1)
    if group_by is not None:
        p_e.reset_index().to_csv(p_e_path, header=group_by + ['p_e'], index=False)
        p_e = dict(p_e)
    else:
        p_e = p_e[0]
        with open(p_e_path, 'w') as f:
            f.write(str(p_e))

    return p_e, p_e_path


# def estimate_p_e(df_rates, p_e_path, conversion='TC', group_by=None):
#     logger.debug(f'p_e estimation will be grouped by {group_by} to {p_e_path}')
#     if group_by is not None:
#         df_rates = df_rates.set_index(group_by)
#     columns = CONVERSION_COLUMNS.copy()
#     columns.remove(conversion)
#     p_e = df_rates[columns[0]]
#     for column in columns[1:]:
#         p_e += df_rates[column]
#     p_e /= len(columns)
#     logger.debug(f'Writing p_e estimates to {p_e_path}')
#     if group_by is not None:
#         p_e.reset_index().to_csv(p_e_path, header=group_by + ['p_e'], index=False)
#     else:
#         p_e = p_e[0]
#         with open(p_e_path, 'w') as f:
#             f.write(str(p_e))
#
#     return p_e, p_e_path
