from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm


def estimate_p_e(df_rates, p_e_path):
    p_e = (df_rates['CT'] + df_rates['GA']) / 2

    # Return a dictionary of index-to-pe pairs if the above is a series
    if isinstance(p_e, pd.Series):
        p_e = p_e.to_dict()
        with open(p_e_path, 'w') as f:
            for key in sorted(p_e.keys()):
                f.write(f'{key},{p_e[key]}\n')
    else:
        with open(p_e_path, 'w') as f:
            f.write(p_e)
    return p_e


def filter_aggregates_part(values, p_e, threshold=0.01):
    # k = values[:,0]
    # n = values[:,1]
    # count = values[:,2]
    mask = np.zeros(values.shape[0], dtype=bool)
    for i, (k, n, count) in enumerate(values):
        left = np.sum(values[(values[:, 0] > k) & (values[:, 1] == n), 2]) * stats.binom.pmf(k, n, p_e)
        right = threshold * count
        if left > right:
            mask[i] = True
    return mask


def filter_aggregates(df_aggregates, p_e, group_by='barcode', n_threads=8):
    df = df_aggregates.sort_values(group_by).reset_index()
    values = df[['TC', 'T', 'count']].values
    groups = df.groupby(group_by).indices

    pbar = tqdm(total=2 * len(groups), ascii=True)

    mask = None
    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        futures = []
        for key in sorted(groups.keys()):
            futures.append(
                executor.submit(
                    filter_aggregates_part, values[groups[key]],
                    p_e.get(key, 0) if isinstance(p_e, dict) else p_e
                )
            )
            pbar.update(1)

        masks = []
        for future in futures:
            masks.append(future.result())

        mask = np.concatenate(masks)

    return df[mask]


def estimate_p_c(df_aggregates, p_e, group_by='barcode', n_threads=8):
    df = filter_aggregates(df_aggregates, p_e, n_threads=n_threads)  # noqa

    # TODO: EM


def estimate_pi_g():
    # TODO
    pass
