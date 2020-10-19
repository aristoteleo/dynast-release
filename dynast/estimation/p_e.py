import logging

import pandas as pd

logger = logging.getLogger(__name__)


def read_p_e(p_e_path, group_by=None):
    if group_by is None:
        with open(p_e_path, 'r') as f:
            return float(f.read())

    df = pd.read_csv(p_e_path)
    return dict(df.set_index(group_by)['p_e'])


def estimate_p_e(df_rates, p_e_path, group_by=None):
    logger.debug(f'p_e estimation will be grouped by {group_by} to {p_e_path}')
    if group_by is not None:
        df_rates = df_rates.set_index(group_by)
    p_e = (df_rates['CT'] + df_rates['GA']) / 2
    logger.debug(f'Writing p_e estimates to {p_e_path}')
    if group_by is not None:
        p_e.reset_index().to_csv(p_e_path, header=group_by + ['p_e'], index=False)
    else:
        p_e = p_e[0]
        with open(p_e_path, 'w') as f:
            f.write(str(p_e))

    return p_e, p_e_path
