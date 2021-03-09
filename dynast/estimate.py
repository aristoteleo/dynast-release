import datetime as dt
import os

import pystan

from . import config, constants, estimation, preprocessing, utils
from .logging import logger
from .stats import Stats


@logger.namespaced('estimate')
def estimate(
    count_dir,
    out_dir,
    reads='complete',
    groups=None,
    cell_threshold=1000,
    cell_gene_threshold=16,
    control_p_e=None,
    control=False,
    n_threads=8,
    temp_dir=None,
    nasc=False,
    seed=None,
):
    stats = Stats()
    stats.start()
    stats_path = os.path.join(
        out_dir, f'{constants.STATS_PREFIX}_{dt.datetime.strftime(stats.start_time, "%Y%m%d_%H%M%S_%f")}.json'
    )
    os.makedirs(out_dir, exist_ok=True)

    conversions = utils.read_pickle(os.path.join(count_dir, constants.CONVS_FILENAME))
    logger.info(f'Conversions: {" ".join(",".join(convs) for convs in conversions)}')
    all_conversions = sorted(utils.flatten_list(conversions))

    counts_path = os.path.join(count_dir, f'{constants.COUNTS_PREFIX}_{"_".join(all_conversions)}.csv')
    logger.info(f'Reading {counts_path}')
    df_counts_uncomplemented = preprocessing.read_counts(counts_path)

    # Check that all requested read groups can be corrected.
    transcriptome_any = df_counts_uncomplemented['transcriptome'].any()
    transcriptome_all = df_counts_uncomplemented['transcriptome'].all()
    if reads != 'complete':
        if 'transcriptome' in reads and not transcriptome_any:
            raise Exception(
                'No reads are assigned to `transcriptome`, so estimation is not supported for this read group.'
            )
        if 'total' in reads and transcriptome_all:
            raise Exception(
                'All reads are assigned to `transcriptome`, so estimation is not supported for `total` read group.'
            )
        for key in set(reads).intersection(('spliced', 'unspliced')):
            if not (df_counts_uncomplemented['velocity'] == key).any():
                raise Exception(
                    f'No reads are assigned to `{key}`, so estimation is not supported for this read group.'
                )
    else:
        reads = []
        if transcriptome_any:
            reads.append('transcriptome')
        if not transcriptome_all:
            reads.append('total')
        reads += list(set(df_counts_uncomplemented['velocity'].unique()) - set(config.VELOCITY_BLACKLIST))
    logger.info(f'Estimation will be done on the following read groups: {reads}')

    gene_infos = utils.read_pickle(os.path.join(count_dir, constants.GENES_FILENAME))
    df_counts_complemented = preprocessing.complement_counts(df_counts_uncomplemented, gene_infos)

    # If cell groups are provided, change barcodes to cell groups
    if groups:
        df_counts_complemented['barcode'] = df_counts_complemented['barcode'].map(groups)
        df_counts_complemented = df_counts_complemented.dropna(subset=['barcode'])

    # Aggregate counts to construct A matrix
    aggregates_paths = {}
    for key in set(reads).union(['transcriptome'] if transcriptome_all else ['total']):
        logger.info(f'Aggregating counts for `{key}`')

        if key == 'transcriptome':
            df = df_counts_complemented[df_counts_complemented['transcriptome']]
        elif key == 'total':
            df = df_counts_complemented
        else:
            df = df_counts_complemented[df_counts_complemented['velocity'] == key]

        for convs in conversions:
            other_convs = list(set(all_conversions) - set(convs))
            aggregates_paths.setdefault(key, {})[tuple(convs)] = preprocessing.aggregate_counts(
                df[(df[other_convs] == 0).all(axis=1)] if other_convs else df,
                os.path.join(out_dir, f'A_{key}_{"_".join(convs)}.csv'),
                conversions=convs
            )

    # Estimate p_e
    p_e_path = os.path.join(out_dir, constants.P_E_FILENAME)
    if control_p_e:
        logger.info('`--p-e` provided. No background mutation rate estimation will be done.')
        df_barcodes = df_counts_complemented['barcode'].drop_duplicates().reset_index(drop=True)
        df_barcodes['p_e'] = control_p_e
        df_barcodes.to_csv(p_e_path, header=['barcode', 'p_e'], index=False)
    else:
        logger.info(f'Estimating average conversion rate in unlabeled RNA to {p_e_path}')
        if control:
            p_e_path = estimation.estimate_p_e_control(
                df_counts_complemented,
                p_e_path,
                conversions=conversions,
            )
        elif nasc:
            p_e_path = estimation.estimate_p_e_nasc(
                preprocessing.read_rates(os.path.join(count_dir, constants.RATES_FILENAME)),
                p_e_path,
                group_by=['barcode'],
            )
        else:
            p_e_path = estimation.estimate_p_e(
                df_counts_complemented,
                p_e_path,
                conversions=conversions,
                group_by=['barcode'],
            )

    if control:
        logger.info('Downstream processing skipped for controls')
        logger.info(f'Use `--p-e {p_e_path}` to run test samples')
        return

    # Estimate p_c
    p_c_paths = {}
    for convs in conversions:
        p_c_path = os.path.join(out_dir, f'{constants.P_C_PREFIX}_{"_".join(convs)}.csv')
        logger.info(f'Estimating {convs} conversion rate in labeled RNA to {p_c_path}')
        df_aggregates = preprocessing.read_aggregates(
            aggregates_paths['transcriptome' if transcriptome_all else 'total'][tuple(convs)]
        )
        p_c_paths[tuple(convs)] = estimation.estimate_p_c(
            df_aggregates,
            estimation.read_p_e(p_e_path, group_by=['barcode']),
            p_c_path,
            group_by=['barcode'],
            threshold=cell_threshold,
            n_threads=n_threads,
            nasc=nasc
        )

    # Estimate pi
    logger.info(f'Compling STAN model from {config.MODEL_PATH}')
    model = pystan.StanModel(file=config.MODEL_PATH, model_name=config.MODEL_NAME)

    pi_paths = {}
    for key in reads:
        for convs in conversions:
            pi_path = os.path.join(out_dir, f'pi_{key}_{"_".join(convs)}.csv')
            logger.info(f'Estimating fraction of labeled `{key}` RNA for conversions {convs} to {pi_path}')
            pi_paths.setdefault(key, {})[tuple(convs)] = estimation.estimate_pi(
                preprocessing.read_aggregates(aggregates_paths[key][tuple(convs)]),
                estimation.read_p_e(p_e_path, group_by=['barcode']),
                estimation.read_p_c(p_c_paths[tuple(convs)], group_by=['barcode']),
                pi_path,
                p_group_by=['barcode'],
                n_threads=n_threads,
                threshold=cell_gene_threshold,
                seed=seed,
                nasc=nasc,
                model=model,
            )
    pis = {key: {convs: estimation.read_pi(pi_paths[key][convs]) for convs in pi_paths[key]} for key in pi_paths}
    if groups:
        group_cells = {}
        for barcode, group in groups.items():
            group_cells.setdefault(group, []).append(barcode)
        pis = {(barcode, gx): value
               for key in pis for convs in pis[key] for (group, gx), value in pis[key][convs].items()
               for barcode in group_cells[group]}

    adata_path = os.path.join(out_dir, constants.ADATA_FILENAME)
    logger.info(f'Combining results into Anndata object at {adata_path}')
    df_counts = preprocessing.complement_counts(df_counts_uncomplemented, gene_infos)
    adata = utils.results_to_adata(
        df_counts[df_counts['barcode'].isin(groups.keys())] if groups else df_counts,
        conversions,
        gene_infos=gene_infos,
        pis=pis
    )
    if groups:
        adata.obs['group'] = adata.obs.index.map(groups).astype('category')
    p_es = estimation.read_p_e(p_e_path, group_by=['barcode'])
    adata.obs['p_e'] = adata.obs['group'].map(p_es) if groups else adata.obs.index.map(p_es)
    for convs in conversions:
        p_cs = estimation.read_p_c(p_c_paths[tuple(convs)], group_by=['barcode'])
        adata.obs[f'p_c_{"_".join(convs)}'] = adata.obs['group'].map(p_cs) if groups else adata.obs.index.map(p_cs)
    adata.write(adata_path)
    stats.end()
    stats.save(stats_path)
