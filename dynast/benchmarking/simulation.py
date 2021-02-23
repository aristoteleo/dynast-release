import random
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pystan

from .. import config, estimation, utils
from ..preprocessing import conversion


def generate_sequence(k, seed=None):
    random.seed(seed)
    return ''.join(random.choices(conversion.BASE_COLUMNS, k=k))


def simulate_reads(sequence, p_e, p_c, pi, l=100, n=100, m=0, seed=None):  # noqa
    generator = np.random.RandomState(seed)

    n_reads = n + m
    n_new = int(n_reads * pi)
    n_old = n_reads - n_new
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
            elif generator.random() < p_e:
                bases = f'{base}{random.choice([b for b in conversion.BASE_COLUMNS if b != base])}'
                if bases != 'TC':
                    conv[bases] += 1
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
            if generator.random() < p_e:
                conv[f'{base}{random.choice([b for b in conversion.BASE_COLUMNS if b != base])}'] += 1
        contents.append(dict(content))
        convs.append(conv)

    df_contents = pd.DataFrame(contents)
    df_conversions = pd.DataFrame(convs)
    df_counts = pd.concat((df_contents, df_conversions), axis=1)[conversion.COLUMNS].iloc[
        generator.choice(np.arange(n_reads), size=n_reads, replace=False)].reset_index()
    df_counts['subset'] = [True] * n + [False] * m

    return df_counts


_model = None


def initializer(model):
    global _model
    _model = model


def estimate(
    df_counts,
    p_e,
    p_c,
    estimate_p_e=False,
    estimate_p_c=False,
    estimate_pi=True,
    model=None,
):
    model = model or _model

    # p_e
    if estimate_p_e:
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
            p_c_path = estimation.estimate_p_c(df_aggregates, p_e_estimate, tf.name)
            p_c_estimate = estimation.read_p_c(p_c_path)
    else:
        p_c_estimate = p_c

    # pi
    if estimate_pi:
        df_aggregates = pd.DataFrame(
            df_counts[df_counts['subset']].groupby(['TC', 'T'], sort=False, observed=True).size()
        )
        df_aggregates.columns = ['count']
        df_aggregates.reset_index(inplace=True)
        df_aggregates = df_aggregates[(df_aggregates[['T', 'count']] > 0).all(axis=1)]
        vals = df_aggregates.values
        guess = min(max((sum(vals[vals[:, 0] > 0][:, 2]) / sum(vals[:, 2])), 0.01), 0.99)
        p_c_estimate, guess, alpha, beta, pi = estimation.pi.fit_stan_mcmc(
            vals,
            p_e_estimate,
            p_c_estimate,
            guess=guess,
            model=model,
            pi_func=estimation.pi.beta_mean,
        )
    else:
        guess, alpha, beta, pi = None, None, None, None

    return p_e_estimate, p_c_estimate, guess, alpha, beta, pi


def _simulate(
    p_e,
    p_c,
    pi,
    sequence=None,
    k=10000,
    l=100,  # noqa
    n=100,
    m=0,
    estimate_p_e=False,
    estimate_p_c=False,
    estimate_pi=True,
    model=None,
    seed=None,
):
    model = model or _model
    sequence = sequence or generate_sequence(k, seed=seed)
    df_counts = simulate_reads(sequence, p_e, p_c, pi, l=l, n=n, m=m, seed=seed)

    return estimate(
        df_counts,
        p_e,
        p_c,
        estimate_p_e=estimate_p_e,
        estimate_p_c=estimate_p_c,
        estimate_pi=estimate_pi,
        model=model,
    )


def simulate(
    p_e,
    p_c,
    pi,
    sequence=None,
    k=10000,
    l=100,  # noqa
    n=100,
    m=0,
    n_runs=16,
    n_threads=8,
    estimate_p_e=False,
    estimate_p_c=False,
    estimate_pi=True,
    model=None,
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
                m=m,
                estimate_p_e=estimate_p_e,
                estimate_p_c=estimate_p_c,
                estimate_pi=estimate_pi,
            ) for _ in range(n_runs)
        ]

        for future in utils.as_completed_with_progress(futures):
            p_e_estimate, p_c_estimate, guess, alpha_estimate, beta_estimate, pi_estimate = future.result()
            p_es.append(p_e_estimate)
            p_cs.append(p_c_estimate)
            guesses.append(guess)
            alphas.append(alpha_estimate)
            betas.append(beta_estimate)
            pis.append(pi_estimate)

    return p_es, p_cs, guesses, alphas, betas, pis
