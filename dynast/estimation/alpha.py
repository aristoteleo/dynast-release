from typing import Dict, FrozenSet, List, Optional, Tuple, Union

import numpy as np
import pandas as pd


def read_alpha(alpha_path: str,
               group_by: Optional[List[str]] = None) -> Union[float, Dict[str, float], Dict[Tuple[str, ...], float]]:
    """Read alpha CSV as a dictionary, with `group_by` columns as keys.

    Args:
        alpha_path: Path to CSV containing alpha values
        group_by: Columns to group by, defaults to `None`

    Returns:
        Dictionary with `group_by` columns as keys (tuple if multiple)
    """
    if group_by is None:
        with open(alpha_path, 'r') as f:
            return float(f.read())

    df = pd.read_csv(alpha_path, dtype={key: 'category' for key in group_by})
    return dict(df.set_index(group_by)['alpha'])


def estimate_alpha(
    df_counts: pd.DataFrame,
    pi_c: Union[float, Dict[str, float], Dict[Tuple[str, ...], float]],
    alpha_path: str,
    conversions: FrozenSet[str] = frozenset({'TC'}),
    group_by: Optional[List[str]] = None,
    pi_c_group_by: Optional[List[str]] = None,
) -> str:
    """Estimate the detection rate alpha.

    Args:
        df_counts: Pandas dataframe containing conversion counts
        pi_c: Labeled mutation rate
        alpha_path: Path to output CSV containing alpha estimates
        conversions: Conversions to consider
        group_by: Columns to group by
        pi_c_group_by: Columns that were used to group when calculating pi_c

    Returns:
        Path to output CSV containing alpha estimates
    """
    columns = ['barcode'] + list(conversions)
    if group_by is not None:
        columns += group_by
    if pi_c_group_by is not None:
        columns += pi_c_group_by
    df_full = df_counts[list(set(columns))].copy()

    if pi_c_group_by is not None:
        df_full.set_index(pi_c_group_by, inplace=True)
        df_full['pi_c'] = df_full.index.map(pi_c)
        df_full.reset_index(inplace=True)
    else:
        df_full['pi_c'] = pi_c
    df_full.dropna(subset=['pi_c'], inplace=True)  # Drop NA values due to pi_c
    pi_cs = df_full['pi_c'].values
    if group_by is None:
        if not isinstance(pi_c, float):
            raise Exception('`pi_c` and `p_c` must be a float when `group_by` is not provided')
        total = df_full.shape[0]
        new = (df_full[list(conversions)] > 0).any(axis=1).sum()
        ntr = new / total
        alpha = pi_c / ntr
    else:
        groupby = df_full.groupby(group_by, sort=False, observed=True)
        groups = groupby.indices
        total = groupby.size()
        new = df_full[(df_full[list(conversions)] > 0).any(axis=1)].groupby(group_by, observed=True, sort=False).size()
        ntr = new.reindex(total.index, fill_value=0.) / total
        alphas = {}
        for key, idx in groups.items():
            pi_c_unique = np.unique(pi_cs[idx])
            if len(pi_c_unique) > 1:
                raise Exception(f'`pi_c` for each aggregate group must be a constant, but instead got {pi_c_unique}.')

            alphas[key] = pi_c_unique[0] / ntr[key]

    with open(alpha_path, 'w') as f:
        if group_by is None:
            f.write(str(alpha))
        else:
            f.write(f'{",".join(group_by)},alpha\n')
            for key in sorted(alphas.keys()):
                alpha = alphas[key]
                formatted_key = key if isinstance(key, str) else ",".join(key)
                f.write(f'{formatted_key},{alpha}\n')
    return alpha_path
