import datetime as dt
import os
import pysam

import anndata
import pandas as pd

from . import config, constants, estimation, preprocessing, utils
from .logging import logger
from .stats import Stats


@logger.namespaced('count')
def count(
    bam_path,
    gtf_path,
    out_dir,
    strand='forward',
    umi_tag=None,
    barcode_tag=None,
    gene_tag='GX',
    barcodes=None,
    control=False,
    quality=27,
    conversions=[['TC']],
    snp_threshold=0.5,
    snp_csv=None,
    correct=None,
    cell_threshold=1000,
    cell_gene_threshold=16,
    control_p_e=None,
    p_group_by=None,
    whitelist_path=None,
    n_threads=8,
    temp_dir=None,
    re=None,
    nasc=False,
    velocity=True,
    seed=None,
):

    def redo(key):
        return re in config.RE_CHOICES[:config.RE_CHOICES.index(key) + 1]

    stats = Stats()
    stats.start()
    stats_path = os.path.join(
        out_dir, f'{constants.STATS_PREFIX}_{dt.datetime.strftime(stats.start_time, "%Y%m%d_%H%M%S_%f")}.json'
    )
    os.makedirs(out_dir, exist_ok=True)

    all_conversions = sorted(utils.flatten_list(conversions))
    correct = correct or []

    # Check memory.
    available_memory = utils.get_available_memory()
    if available_memory < config.RECOMMENDED_MEMORY:
        logger.warning(
            f'There is only {available_memory / (1024 ** 3):.2f} GB of free memory on the machine. '
            f'It is highly recommended to have at least {config.RECOMMENDED_MEMORY // (1024 ** 3)} GB '
            'free when running dynast. Continuing may cause dynast to crash with an out-of-memory error.'
        )

    # Check that BAM is sorted by coordinate. If not, run samtools sort.
    sorted_bam_path = '{}.sortedByCoord{}'.format(*os.path.splitext(bam_path))
    bam_sorted = False
    with pysam.AlignmentFile(bam_path, 'rb') as f:
        if f.header.get('HD', {}).get('SO') == 'coordinate':
            bam_sorted = True
    skip = utils.all_exists([sorted_bam_path]) or bam_sorted
    with stats.step('sort', skipped=skip) or redo('sort'):
        if not skip:
            logger.info(f'Sorting {bam_path} with samtools to {sorted_bam_path}')
            pysam.sort(bam_path, '-o', sorted_bam_path, '-@', str(n_threads))
            bam_path = sorted_bam_path

    # Check if BAM index exists and create one if it doesn't.
    bai_path = f'{bam_path}.bai'
    skip = utils.all_exists([bai_path]) and not redo('index')
    with stats.step('index', skipped=skip):
        if not skip:
            logger.info(f'Indexing {bam_path} with samtools to {bai_path}')
            pysam.index(bam_path, bai_path, '-@', str(n_threads))

    # Parse BAM and save results
    parse_dir = os.path.join(out_dir, constants.PARSE_DIR)
    conversions_path = os.path.join(parse_dir, constants.CONVERSIONS_FILENAME)
    conversions_index_path = os.path.join(parse_dir, constants.CONVERSIONS_INDEX_FILENAME)
    no_conversions_path = os.path.join(parse_dir, constants.NO_CONVERSIONS_FILENAME)
    no_index_path = os.path.join(parse_dir, constants.NO_CONVERSIONS_INDEX_FILENAME)
    genes_path = os.path.join(parse_dir, constants.GENES_FILENAME)
    transcripts_path = os.path.join(parse_dir, constants.TRANSCRIPTS_FILENAME)
    conversions_required = [
        conversions_path, conversions_index_path, no_conversions_path, no_index_path, genes_path, transcripts_path
    ]
    gene_infos = None
    transcript_infos = None
    skip = utils.all_exists(conversions_required) and not redo('parse')
    with stats.step('parse', skipped=skip), logger.namespaced_context('parse'):
        if not skip:
            re = 'parse'
            os.makedirs(parse_dir, exist_ok=True)
            logger.info('Parsing gene and transcript information from GTF')
            gene_infos, transcript_infos = preprocessing.parse_gtf(gtf_path)
            utils.write_pickle(gene_infos, genes_path)
            utils.write_pickle(transcript_infos, transcripts_path)

            logger.info(f'Parsing read conversion information from BAM to {conversions_path}')
            conversions_path, index_path, no_conversions_path, no_index_path = preprocessing.parse_all_reads(
                bam_path,
                conversions_path,
                conversions_index_path,
                no_conversions_path,
                no_index_path,
                gene_infos,
                transcript_infos,
                strand=strand,
                umi_tag=umi_tag,
                barcode_tag=barcode_tag,
                gene_tag=gene_tag,
                barcodes=barcodes,
                n_threads=n_threads,
                temp_dir=temp_dir,
                nasc=nasc,
                velocity=velocity
            )
        else:
            logger.info('Skipped')

    # Detect SNPs
    snp_dir = os.path.join(out_dir, constants.SNP_DIR)
    coverage_path = os.path.join(snp_dir, constants.COVERAGE_FILENAME)
    coverage_index_path = os.path.join(snp_dir, constants.COVERAGE_INDEX_FILENAME)
    snps_path = os.path.join(snp_dir, constants.SNPS_FILENAME)
    snps_required = [coverage_path, coverage_index_path, snps_path]
    skip = not snp_threshold or (utils.all_exists(snps_required) and not redo('snp'))
    with stats.step('snp', skipped=skip), logger.namespaced_context('snp'):
        if not skip:
            re = 'snp'
            logger.info(f'Calculating coverage and outputting to {coverage_path}')
            os.makedirs(snp_dir, exist_ok=True)
            coverage_path, coverage_index_path = preprocessing.calculate_coverage(
                bam_path, {
                    contig: set(df_part['genome_i'])
                    for contig, df_part in preprocessing.read_conversions(
                        conversions_path, usecols=['contig', 'genome_i']
                    ).drop_duplicates().groupby('contig', sort=False, observed=True)
                },
                coverage_path,
                coverage_index_path,
                umi_tag=umi_tag,
                barcode_tag=barcode_tag,
                gene_tag=gene_tag,
                barcodes=barcodes,
                n_threads=n_threads,
                temp_dir=temp_dir,
                velocity=velocity
            )

            logger.info(f'Detecting SNPs with threshold {snp_threshold} to {snps_path}')
            snps_path = preprocessing.detect_snps(
                conversions_path,
                conversions_index_path,
                coverage_path,
                coverage_index_path,
                snps_path,
                quality=quality,
                threshold=snp_threshold,
                n_threads=n_threads,
            )
        else:
            logger.info('Skipped')

    # Count conversions and calculate mutation rates
    count_dir = os.path.join(out_dir, constants.COUNT_DIR)
    counts_path = os.path.join(count_dir, f'{constants.COUNTS_PREFIX}_{"_".join(all_conversions)}.csv')
    count_required = [counts_path]
    skip = utils.all_exists(count_required) and not redo('count')
    with stats.step('count', skipped=skip), logger.namespaced_context('count'):
        if not skip:
            re = 'count'
            os.makedirs(count_dir, exist_ok=True)
            logger.info(f'Counting conversions to {conversions_path}')
            snps = utils.merge_dictionaries(
                preprocessing.read_snps(snps_path) if snp_threshold else {},
                preprocessing.read_snp_csv(snp_csv) if snp_csv else {},
                f=set.union,
                default=set,
            )
            counts_path = preprocessing.count_conversions(
                conversions_path,
                conversions_index_path,
                no_conversions_path,
                no_index_path,
                counts_path,
                barcodes=barcodes,
                snps=snps,
                quality=quality,
                conversions=all_conversions,
                n_threads=n_threads,
                temp_dir=temp_dir
            )
        else:
            logger.info('Skipped')
    if barcodes:
        count_barcodes = set(preprocessing.read_counts(counts_path, usecols=['barcode'])['barcode'])
        extra_barcodes = count_barcodes - set(barcodes)
        if extra_barcodes:
            logger.warning(
                f'{counts_path} has {len(extra_barcodes)} barcodes not in the filter list provided with `--barcodes`. '
                'Re-run dynast with `--re count` to fix this inconsistency. '
                'Otherwise, all of these extra barcodes will be included in following steps.'
            )
            barcodes = sorted(count_barcodes.union(barcodes))

    aggregates_dir = os.path.join(out_dir, constants.AGGREGATES_DIR)
    rates_path = os.path.join(aggregates_dir, constants.RATES_FILENAME)
    velocities = list(preprocessing.read_counts(counts_path, usecols=['velocity'])['velocity'].unique()
                      ) if velocity else []
    aggregates_paths = {
        key: {tuple(convs): os.path.join(aggregates_dir, f'{key}_{"_".join(convs)}.csv')
              for convs in conversions}
        for key in set(correct).intersection(velocities + ['total', 'transcriptome']
                                             ).union({'total'} if velocity else {'transcriptome'})
    } if correct else {}
    aggregates_required = utils.flatten_dict_values(aggregates_paths) + [rates_path]
    df_counts_uncomplemented = None
    df_counts_complemented = None
    skip = utils.all_exists(aggregates_required) and not redo('aggregate')
    with stats.step('aggregate', skipped=skip), logger.namespaced_context('aggregate'):
        if not skip:
            re = 'aggregate'
            logger.info(f'Computing mutation rates and outputting to {rates_path}')
            os.makedirs(aggregates_dir, exist_ok=True)
            gene_infos = gene_infos or utils.read_pickle(genes_path)
            df_counts_uncomplemented = preprocessing.read_counts(counts_path)
            df_counts_complemented = preprocessing.complement_counts(df_counts_uncomplemented, gene_infos)
            rates_path = preprocessing.calculate_mutation_rates(
                df_counts_uncomplemented if nasc else df_counts_complemented, rates_path, group_by=p_group_by
            )

            if correct:
                # Total/transcriptome reads are required to calculate p_e/p_c later on
                dfs = {
                    'total': df_counts_complemented
                } if velocity else {
                    'transcriptome': df_counts_complemented[df_counts_complemented['transcriptome']]
                }
                for key in correct:
                    if key in ('total', 'transcriptome', 'ambiguous'):
                        continue
                    df_key = df_counts_complemented[df_counts_complemented['velocity'] == key]
                    if df_key.shape[0] == 0:
                        logger.warning(f'No reads were assigned to `{key}`. No correction can be done.')
                    else:
                        dfs.update({key: df_key})

                for key, df in dfs.items():
                    for convs in conversions:
                        aggregates_path = aggregates_paths[key][tuple(convs)]
                        logger.info(
                            f'Aggregating counts for `{key}` reads for conversions {convs} to {aggregates_path}'
                        )
                        # Ignore reads that have more than one conversion of interest
                        other_convs = list(set(all_conversions) - set(convs))
                        aggregates_paths[key][tuple(convs)] = preprocessing.aggregate_counts(
                            df[(df[other_convs] == 0).all(axis=1)] if other_convs else df,
                            aggregates_path,
                            conversions=convs
                        )
        else:
            logger.info('Skipped')

    velocity_blacklist = ['unassigned', 'ambiguous']
    estimates_dir = os.path.join(out_dir, constants.ESTIMATES_DIR)
    p_e_path = os.path.join(estimates_dir, constants.P_E_FILENAME)
    p_c_paths = {tuple(convs): os.path.join(estimates_dir, f'p_c_{"_".join(convs)}.csv') for convs in conversions}
    pi_paths = {
        key: {tuple(convs): os.path.join(estimates_dir, f'{key}_{"_".join(convs)}.csv')
              for convs in conversions}
        for key in aggregates_paths.keys()
        if key not in velocity_blacklist
    }
    estimates_paths = [p_e_path] + (
        utils.flatten_dict_values(p_c_paths) + utils.flatten_dict_values(pi_paths) if not control else []
    )
    skip = (not control and not correct) or (utils.all_exists(estimates_paths) and not redo('estimate'))
    with stats.step('estimate', skipped=skip), logger.namespaced_context('estimate'):
        if not skip:
            re = 'estimate'
            os.makedirs(estimates_dir, exist_ok=True)
            if control_p_e:
                logger.info('`--p-e` provided. No background mutation rate estimation will be done.')
                if p_group_by is not None:
                    df_barcodes = preprocessing.read_counts(
                        counts_path, usecols=p_group_by
                    ).drop_duplicates().reset_index(drop=True)
                    df_barcodes['p_e'] = control_p_e
                    df_barcodes.to_csv(p_e_path, header=p_group_by + ['p_e'], index=False)
                else:
                    with open(p_e_path, 'w') as f:
                        f.write(str(control_p_e))

            else:
                logger.info(f'Estimating average mismatch rate in unlabeled RNA to {p_e_path}')
                if control:
                    gene_infos = gene_infos or utils.read_pickle(genes_path)
                    p_e_path = estimation.estimate_p_e_control(
                        df_counts_complemented if df_counts_complemented is not None else
                        preprocessing.complement_counts(preprocessing.read_counts(counts_path), gene_infos),
                        p_e_path,
                        conversions=conversions,
                    )
                elif nasc:
                    p_e_path = estimation.estimate_p_e_nasc(
                        preprocessing.read_rates(rates_path),
                        p_e_path,
                        group_by=p_group_by,
                    )
                else:
                    gene_infos = gene_infos or utils.read_pickle(genes_path)
                    p_e_path = estimation.estimate_p_e(
                        df_counts_complemented if df_counts_complemented is not None else
                        preprocessing.complement_counts(preprocessing.read_counts(counts_path), gene_infos),
                        p_e_path,
                        conversions=conversions,
                        group_by=p_group_by,
                    )

            if not control:
                for convs in conversions:
                    p_c_path = p_c_paths[tuple(convs)]
                    logger.info(
                        f'Estimating average mismatch rate in labeled RNA for conversions {convs} to {p_c_path}'
                    )
                    df_aggregates = preprocessing.read_aggregates(
                        aggregates_paths['total' if velocity else 'transcriptome'][tuple(convs)]
                    )
                    p_c_paths[tuple(convs)] = estimation.estimate_p_c(
                        df_aggregates,
                        estimation.read_p_e(p_e_path, group_by=p_group_by),
                        p_c_path,
                        group_by=p_group_by,
                        threshold=cell_threshold,
                        n_threads=n_threads,
                        nasc=nasc
                    )

                    for key, paths in aggregates_paths.items():
                        if key in velocity_blacklist or key not in correct:
                            continue
                        pi_path = pi_paths[key][tuple(convs)]
                        logger.info(f'Estimating fraction of labeled `{key}` RNA for conversions {convs} to {pi_path}')
                        pi_paths[key][tuple(convs)] = estimation.estimate_pi(
                            preprocessing.read_aggregates(paths[tuple(convs)]),
                            estimation.read_p_e(p_e_path, group_by=p_group_by),
                            estimation.read_p_c(p_c_paths[tuple(convs)], group_by=p_group_by),
                            pi_path,
                            p_group_by=p_group_by,
                            n_threads=n_threads,
                            threshold=cell_gene_threshold,
                            seed=seed,
                            nasc=nasc
                        )
        else:
            logger.info('Skipped')

    adata_path = os.path.join(out_dir, constants.ADATA_FILENAME)
    skip = control or (utils.all_exists([adata_path]) and not redo('split'))
    with stats.step('split', skipped=skip), logger.namespaced_context('split'):
        if not skip:
            re = 'split'
            logger.info(f'Combining results into an Anndata object at {adata_path}')
            gene_infos = gene_infos or utils.read_pickle(genes_path)

            logger.info('Loading results for `transcriptome` reads')
            df_counts = (
                df_counts_complemented if df_counts_complemented is not None else
                preprocessing.complement_counts(preprocessing.read_counts(counts_path), gene_infos)
            )
            df_counts_transcriptome = df_counts[df_counts['transcriptome']]
            barcodes = barcodes or sorted(df_counts['barcode'].unique())
            features = sorted(df_counts['GX'].unique())
            names = [gene_infos.get(feature, {}).get('gene_name') for feature in features]

            # Deal with transcriptome reads first because those need to be
            # in adata.X
            matrix = preprocessing.counts_to_matrix(df_counts_transcriptome, barcodes, features)
            layers = {}
            for convs in conversions:
                # Ignore reads that have other conversions
                other_convs = list(set(all_conversions) - set(convs))
                join = '_'.join(convs)

                # Counts for transcriptome reads (i.e. X_unlabeled + X_labeled = X)
                layers[f'X_n_{join}'], layers[f'X_l_{join}'] = preprocessing.split_counts(
                    df_counts_transcriptome[(df_counts_transcriptome[other_convs] == 0).all(axis=1)],
                    barcodes,
                    features,
                    conversions=convs
                )
                if 'transcriptome' in correct:
                    pis = estimation.read_pi(pi_paths['transcriptome'][tuple(convs)])
                    (
                        _,
                        layers[f'X_n_{join}_est'],
                        layers[f'X_l_{join}_est'],
                    ) = estimation.split_matrix(layers[f'X_n_{join}'] + layers[f'X_l_{join}'], pis, barcodes, features)

                if velocity:
                    logger.info('Loading results for `total` reads')

                    # Counts for all reads, regardless of whether they map to the transcriptome
                    layers[f'unlabeled_{join}'], layers[f'labeled_{join}'] = preprocessing.split_counts(
                        df_counts[(df_counts[other_convs] == 0).all(axis=1)], barcodes, features, conversions=convs
                    )
                    if 'total' in correct:
                        pis = estimation.read_pi(pi_paths['total'][tuple(convs)])
                        (
                            _,
                            layers[f'unlabeled_{join}_est'],
                            layers[f'labeled_{join}_est'],
                        ) = estimation.split_matrix(
                            layers[f'unlabeled_{join}'] + layers[f'labeled_{join}'], pis, barcodes, features
                        )

            # Velocities
            for key in velocities:
                logger.info(f'Loading results for `{key}` reads')
                df_counts_velocity = df_counts[df_counts['velocity'] == key]
                layers[key] = preprocessing.counts_to_matrix(df_counts_velocity, barcodes, features)
                if key in velocity_blacklist:
                    continue
                for convs in conversions:
                    other_convs = list(set(all_conversions) - set(convs))
                    join = '_'.join(convs)
                    layers[f'{key[0]}n_{join}'], layers[f'{key[0]}l_{join}'] = preprocessing.split_counts(
                        df_counts_velocity[(df_counts_velocity[other_convs] == 0).all(axis=1)],
                        barcodes,
                        features,
                        conversions=convs
                    )
                    if key in correct:
                        pis = estimation.read_pi(pi_paths[key][tuple(convs)])
                        (
                            _,
                            layers[f'{key[0]}n_{join}_est'],
                            layers[f'{key[0]}l_{join}_est'],
                        ) = estimation.split_matrix(
                            layers[f'{key[0]}n_{join}'] + layers[f'{key[0]}l_{join}'], pis, barcodes, features
                        )

            # Construct adata with umi counts as layers
            adata = anndata.AnnData(
                X=matrix,
                obs=pd.DataFrame(index=pd.Series(barcodes, name='barcodes')),
                var=pd.DataFrame(index=pd.Series(features, name='gene_id'), data={'gene_name': pd.Categorical(names)}),
                layers=layers
            )
            adata.write(adata_path, compression='gzip')
        else:
            logger.info('Skipped')

    if control:
        logger.info('Downstream processing skipped for controls')
        if snp_threshold:
            logger.info(f'Use `--snp-csv {snps_path}` to run test samples')
        logger.info(f'Use `--p-e {p_e_path}` for test samples')
    stats.end()
    stats.save(stats_path)
