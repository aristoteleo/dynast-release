import random
import tempfile
from collections import Counter
from concurrent.futures import as_completed, ProcessPoolExecutor

import numpy as np
import pandas as pd
import pystan

from .. import config, estimation
from ..preprocessing import conversion


def generate_sequence(k, seed=None):
    random.seed(seed)
    return ''.join(random.choices(conversion.BASE_COLUMNS, k=k))


def simulate_reads(sequence, p_e, p_c, pi, l=100, n=100, seed=None):  # noqa
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
                if generator.random() < p_e:
                    conv[f'{base}{random.choice([b for b in conversion.BASE_COLUMNS if b != base])}'] += 1
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
    df_counts = pd.concat((df_contents, df_conversions), axis=1)[conversion.COLUMNS]

    return df_counts


_model = None


def initializer(model):
    global _model
    _model = model


def estimate(
    df_counts,
    p_e,
    p_c,
    fit,
    pi_func,
    estimate_p_e=False,
    estimate_p_c=False,
    model=None,
):
    model = model or _model

    # p_e
    with tempfile.NamedTemporaryFile() as tf:
        p_e_estimate = estimation.estimate_p_e(df_counts, tf.name)[0]
    if not estimate_p_e:
        p_e_estimate = p_e

    # p_c
    df_aggregates = pd.DataFrame(df_counts.groupby(['TC', 'T']).size())
    df_aggregates.columns = ['count']
    df_aggregates.reset_index(inplace=True)
    with tempfile.NamedTemporaryFile() as tf1, tempfile.NamedTemporaryFile() as tf2:
        p_c_estimate = estimation.estimate_p_c(df_aggregates, p_e_estimate, tf1.name, tf2.name)[0]
    if not estimate_p_c:
        p_c_estimate = p_c

    # pi
    alpha, beta, pi, _ = fit(
        df_aggregates.values,
        p_e_estimate,
        p_c_estimate,
        model=model,
        pi_func=pi_func,
    )

    return p_e_estimate, p_c_estimate, alpha, beta, pi


def _simulate(
    p_e,
    p_c,
    pi,
    fit,
    pi_func,
    sequence=None,
    k=10000,
    l=100,  # noqa
    n=100,
    estimate_p_e=False,
    estimate_p_c=False,
    model=None,
    seed=None,
):
    sequence = sequence or generate_sequence(k, seed=seed)
    df_counts = simulate_reads(sequence, p_e, p_c, pi, l=l, n=n, seed=seed)
    return estimate(
        df_counts,
        p_e,
        p_c,
        fit,
        pi_func,
        estimate_p_e=False,
        estimate_p_c=False,
        model=model,
    )


def simulate(
    p_e,
    p_c,
    pi,
    fit,
    pi_func,
    sequence=None,
    k=10000,
    l=100,  # noqa
    n=100,
    n_runs=16,
    n_threads=8,
    estimate_p_e=False,
    estimate_p_c=False,
    model=None,
):
    model = model or pystan.StanModel(file=config.MODEL_PATH, model_name=config.MODEL_NAME)

    p_es = []
    p_cs = []
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
                fit,
                pi_func,
                sequence=sequence,
                k=k,
                l=l,
                n=n,
                estimate_p_e=estimate_p_e,
                estimate_p_c=estimate_p_c
            ) for _ in range(n_runs)
        ]

        for future in as_completed(futures):
            p_e_estimate, p_c_estimate, alpha_estimate, beta_estimate, pi_estimate = future.result()
            p_es.append(p_e_estimate)
            p_cs.append(p_c_estimate)
            alphas.append(alpha_estimate)
            betas.append(beta_estimate)
            pis.append(pi_estimate)

    return p_es, p_cs, alphas, betas, pis
