import random
import tempfile
from collections import Counter
from concurrent.futures import as_completed, ProcessPoolExecutor
from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pystan
from tqdm import tqdm

from .. import config, estimation
from ..preprocessing import aggregation, conversion


def generate_sequence(k, seed=None):
    """Generate a random genome sequence of length `k`.

    :param k: length of the sequence
    :type k: int
    :param seed: random seed, defaults to `None`
    :type seed: int, optional

    :return: a random sequence
    :rtype: str
    """
    random.seed(seed)
    return ''.join(random.choices(conversion.BASE_COLUMNS, k=k))


def simulate_reads(sequence, p_e, p_c, pi, l=100, n=100, seed=None):  # noqa
    """Simulate `n` reads of length `l` from a sequence.

    :param sequence: sequence to generate the reads from
    :type sequence: str
    :param p_e: background specific mutation rate. This is the rate a specific
                base mutates to another specific base (i.e. T>C, A>G, ...)
    :type p_e: float
    :param p_c: T>C mutation rate in labeled RNA
    :type p_c: float
    :param pi: fraction of labeled RNA
    :type pi: float
    :param l: length of each read, defaults to `100`
    :type l: int, optional
    :param n: number of reads to simulate, defaults to `100`
    :type n: int, optional
    :param seed: random seed, defaults to `None`
    :type seed: int, optional

    :return: a dataframe with each read as a row and the number of conversions and
             base content as the columns
    :rtype: pandas.DataFrame
    """
    generator = np.random.RandomState(seed)

    n_new = int(n * pi)
    n_old = n - n_new
    contents = []
    convs = []

    # Generate new sequences
    for _ in range(n_new):
        i = generator.randint(0, len(sequence) - l)  # noqa
        subsequence = sequence[i:i + l]  # noqa

        # Nucleotide content
        content = Counter(subsequence)

        # Create read with mutations
        conv = {conversion: 0 for conversion in conversion.CONVERSION_COLUMNS}
        for base in subsequence:
            if base == 'T' and generator.random() < p_c:
                conv['TC'] += 1
            else:
                other_bases = [b for b in conversion.BASE_COLUMNS if b != base]
                for other_base in random.sample(other_bases, k=len(other_bases)):
                    bases = f'{base}{other_base}'
                    if generator.random() < p_e and bases != 'TC':
                        conv[bases] += 1
                        break

        contents.append(dict(content))
        convs.append(conv)

    # Generate old sequences
    for _ in range(n_old):
        i = generator.randint(0, len(sequence) - l)  # noqa
        subsequence = sequence[i:i + l]  # noqa

        # Nucleotide content
        content = Counter(subsequence)

        # Create read with mutations
        conv = {conversion: 0 for conversion in conversion.CONVERSION_COLUMNS}
        for base in subsequence:
            other_bases = [b for b in conversion.BASE_COLUMNS if b != base]
            for other_base in random.sample(other_bases, k=len(other_bases)):
                if generator.random() < p_e:
                    conv[f'{base}{other_base}'] += 1
                    break
        contents.append(dict(content))
        convs.append(conv)

    df_contents = pd.DataFrame(contents)
    df_conversions = pd.DataFrame(convs)
    df_counts = pd.concat((df_contents, df_conversions),
                          axis=1)[conversion.COLUMNS].iloc[generator.choice(np.arange(n), size=n,
                                                                            replace=False)].reset_index(drop=True)

    return df_counts


__model = None
_pi_model = None


def initializer(model):
    global _model
    _model = model


def estimate(
    df_counts,
    p_e,
    p_c,
    pi,
    estimate_p_e=False,
    estimate_p_c=False,
    estimate_pi=True,
    model=None,
    nasc=False,
):
    """
    """
    # p_e
    if estimate_p_e:
        if nasc:
            with tempfile.NamedTemporaryFile() as tf:
                rates_path = aggregation.calculate_mutation_rates(df_counts, tf.name)
                df_rates = aggregation.read_rates(rates_path)
            with tempfile.NamedTemporaryFile() as tf:
                p_e_path = estimation.estimate_p_e_nasc(df_rates, tf.name)
                p_e_estimate = estimation.read_p_e(p_e_path)

        else:
            with tempfile.NamedTemporaryFile() as tf:
                p_e_path = estimation.estimate_p_e(df_counts, tf.name)
                p_e_estimate = estimation.read_p_e(p_e_path)
    else:
        p_e_estimate = p_e

    # p_c
    if estimate_p_c:
        df_aggregates = pd.DataFrame(df_counts.groupby(['TC', 'T'], sort=False, observed=True).size())
        df_aggregates.columns = ['count']
        df_aggregates.reset_index(inplace=True)
        df_aggregates.rename(columns={'TC': 'conversion', 'T': 'base'}, inplace=True)
        with tempfile.NamedTemporaryFile() as tf:
            p_c_path = estimation.estimate_p_c(df_aggregates, p_e_estimate, tf.name, nasc=nasc)
            p_c_estimate = estimation.read_p_c(p_c_path)
    else:
        p_c_estimate = p_c

    # pi
    if estimate_pi:
        df_aggregates = pd.DataFrame(
            df_counts[df_counts['GX'] == 'gene_0'].groupby(['TC', 'T'], sort=False, observed=True).size()
        )
        df_aggregates.columns = ['count']
        df_aggregates.reset_index(inplace=True)
        df_aggregates = df_aggregates[(df_aggregates[['T', 'count']] > 0).all(axis=1)]
        vals = df_aggregates.values
        guess = min(max((sum(vals[vals[:, 0] > 0][:, 2]) / sum(vals[:, 2])), 0.01), 0.99)
        guess, alpha, beta, pi_estimate = estimation.pi.fit_stan_mcmc(
            vals,
            p_e_estimate,
            p_c_estimate,
            guess=guess,
            model=model,
        )
        if nasc:
            pi_estimate = estimation.pi.beta_mode(alpha, beta)
    else:
        guess, alpha, beta, pi_estimate = pi, None, None, pi

    return p_e_estimate, p_c_estimate, guess, alpha, beta, pi_estimate


def _simulate(
    p_e,
    p_c,
    pi,
    sequence=None,
    k=10000,
    l=100,  # noqa
    n=100,
    estimate_p_e=False,
    estimate_p_c=False,
    estimate_pi=True,
    seed=None,
    model=None,
    nasc=False,
):
    model = model or _model

    pis = pi
    ns = n
    if isinstance(pi, list) and not isinstance(n, list):
        pis = pi
        ns = [n] * len(pis)
    elif not isinstance(pi, list) and isinstance(n, list):
        ns = n
        pis = [pi] * len(ns)
    elif not isinstance(pi, list) and not isinstance(n, list):
        ns = [n]
        pis = [pi]

    assert len(pis) == len(ns)

    dfs = []
    for i, (pi, n) in enumerate(zip(pis, ns)):
        sequence = sequence or generate_sequence(k, seed=seed)

        df_counts = simulate_reads(sequence, p_e, p_c, pi, l=l, n=n, seed=seed)
        df_counts['GX'] = f'gene_{i}'
        dfs.append(df_counts)
    df_counts = pd.concat(dfs, ignore_index=True)

    return estimate(
        df_counts,
        p_e,
        p_c,
        pis[0],
        estimate_p_e=estimate_p_e,
        estimate_p_c=estimate_p_c,
        estimate_pi=estimate_pi,
        model=model,
        nasc=nasc,
    )


def simulate(
    p_e,
    p_c,
    pi,
    sequence=None,
    k=10000,
    l=100,  # noqa
    n=100,
    n_runs=16,
    n_threads=8,
    estimate_p_e=False,
    estimate_p_c=False,
    estimate_pi=True,
    model=None,
    nasc=False,
):
    model = model or pystan.StanModel(file=config.MODEL_PATH, model_name=config.MODEL_NAME)

    p_es = []
    p_cs = []
    guesses = []
    alphas = []
    betas = []
    pis = []
    with ProcessPoolExecutor(max_workers=n_threads, initializer=initializer, initargs=(model,)) as executor:
        futures = [
            executor.submit(
                _simulate,
                p_e,
                p_c,
                pi,
                sequence=sequence,
                k=k,
                l=l,
                n=n,
                estimate_p_e=estimate_p_e,
                estimate_p_c=estimate_p_c,
                estimate_pi=estimate_pi,
                nasc=nasc,
            ) for _ in range(n_runs)
        ]

        for future in as_completed(futures):
            p_e_estimate, p_c_estimate, guess, alpha_estimate, beta_estimate, pi_estimate = future.result()
            p_es.append(p_e_estimate)
            p_cs.append(p_c_estimate)
            guesses.append(guess)
            alphas.append(alpha_estimate)
            betas.append(beta_estimate)
            pis.append(pi_estimate)

    return p_es, p_cs, guesses, alphas, betas, pis


def simulate_batch(
    p_e,
    p_c,
    pi,
    l,  # noqa
    n,
    estimate_p_e,
    estimate_p_c,
    estimate_pi,
    n_runs,
    n_threads,
    model,
    nasc=False
):
    """Helper function to run simulations in batches.
    """
    p_es = p_e
    p_cs = p_c
    pis = pi
    if not isinstance(p_e, list):
        p_es = [p_e]
    if not isinstance(p_c, list):
        p_cs = [p_c]
    if not isinstance(pi, list):
        pis = [pi]

    dfs = []
    for p_e, p_c, pi in tqdm(list(product(p_es, p_cs, pis))):
        p_e_estimates, p_c_estimates, guesses, alphas, betas, pi_estimates = simulate(
            p_e,
            p_c,
            pi,
            l=l,
            n=n,
            estimate_p_e=estimate_p_e,
            estimate_p_c=estimate_p_c,
            estimate_pi=estimate_pi,
            n_runs=n_runs,
            n_threads=n_threads,
            model=model,
            nasc=nasc,
        )

        dfs.append(
            pd.DataFrame({
                'p_e': p_e,
                'p_c': p_c,
                'pi': pi[0] if isinstance(pi, list) else pi,
                'p_e_estimate': p_e_estimates,
                'p_c_estimate': p_c_estimates,
                'guess': guesses,
                'alpha_estimate': alphas,
                'beta_estimate': betas,
                'pi_estimate': pi_estimates
            })
        )
    return pd.concat(dfs, ignore_index=True)


def plot_estimations(
    X, Y, n_runs, means, truth, ax=None, box=True, tick_decimals=1, title=None, xlabel=None, ylabel=None
):
    if ax is not None:
        _ax = ax
    else:
        fig, _ax = plt.subplots(figsize=(5, 5), tight_layout=True)

    if box:
        X_range = max(X) - min(X)
        _ax.boxplot(
            list(np.array(Y).reshape(-1, n_runs)),
            positions=np.sort(np.unique(X)),
            zorder=-1,
            widths=X_range * 0.05,
            medianprops=dict(c='gray', linewidth=1.5),
            boxprops=dict(facecolor='lightgray', color='gray', linewidth=1.5),
            whiskerprops=dict(c='gray', linewidth=1.5),
            capprops=dict(c='gray', linewidth=1.5),
            patch_artist=True,
            showfliers=False,
        )
    _ax.scatter(X, Y, s=3, label=f'n={n_runs}')
    _ax.scatter(means.index, means.values, s=15, label='mean')
    try:
        iter(truth)
        _ax.plot(truth, truth, c='red', linewidth=1, label='truth')
    except:  # noqa
        _ax.plot([min(X), max(X)], [truth, truth], c='red', linewidth=1, label='truth')

    if box:
        _ax.set_xlim(left=min(X) - X_range * 0.1, right=max(X) + X_range * 0.1)
    xticks = np.sort(np.unique(X))
    _ax.set_xticks(xticks)
    _ax.set_xticklabels([f'{round(x, tick_decimals)}' for x in xticks])
    _ax.legend()
    if title:
        _ax.set_title(title)
    if xlabel:
        _ax.set_xlabel(xlabel)
    if ylabel:
        _ax.set_ylabel(ylabel)
    if ax is None:
        fig.show()
    return _ax
