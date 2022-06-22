from typing import Dict, FrozenSet, List, Optional, Tuple, Union

import pandas as pd


def read_alpha(alpha_path: str, group_by: Optional[List[str]] = None) -> Union[float, Dict[str, float], Dict[Tuple[str, ...], float]]:
    """Read alpha CSV as a dictionary, with `group_by` columns as keys.

    Args:
        alpha_path: Path to CSV containing alpha values
        group_by: Columns to group by, defaults to `None`

    Returns:
        Dictionary with `group_by` columns as keys (tuple if multiple)
    """
    if group_by is None:
        with open(p_c_path, 'r') as f:
            return float(f.read())

    df = pd.read_csv(p_c_path, dtype={key: 'string' for key in group_by})
    return dict(df.set_index(group_by)['p_c'])

def estimate_alpha(
    df_counts: pd.DataFrame,
    pi_c: Union[float, Dict[str, float], Dict[Tuple[str, ...], float]],
    alpha_path: str,
    conversions: FrozenSet[str] = frozenset({'TC'}),
    pi_group_by: Optional[List[str]] = None,
    group_by: Optional[List[str]] = None,
) -> str:
    """Estimate the detection rate alpha.

    Args:
        df_counts: Pandas dataframe containing conversion counts
        pi_c: Labeled mutation rate
        alpha: Path to output CSV containing alpha estimates
        conversions: Conversions to consider
        pi_group_by: Columns that were used to group when calculating pi_c
        group_by: Columns to group by

    Returns:
        Path to output CSV containing alpha estimates
    """
    if group_by is not None:
        total = df_counts.groupby(group_by, observed=True, sort=False).size()
        new = df_counts[df_counts[list(conversions)] > 0].groupby(group_by, observed=True, sort=False).size()
    else:
        total = df_counts.shape[0]
        new = (df_counts[list(conversions)] > 0).sum()
    ntr = new / total
    if isinstance(pi_c, dict):
        pi_c = pd.Series(pi_c)
    alpha = pi_c / ntr

    with open(alpha_path, 'w') as f:
        if group_by is None:
            f.write(str(alpha))
        else:
            f.write(f'{",".join(group_by)},alpha\n')
            for key in sorted(alpha.keys()):
                a = alpha[key]
                formatted_key = key if isinstance(key, str) else ",".join(key)
                f.write(f'{formatted_key},{a}\n')



