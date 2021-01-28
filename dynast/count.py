import datetime as dt
import logging
import os
import pysam

import anndata
import pandas as pd

from . import config, constants, estimation, preprocessing, utils
from .stats import STATS

logger = logging.getLogger(__name__)


def count(
    bam_path,
    gtf_path,
    out_dir,
    strand='forward',
    umi_tag=None,
    barcode_tag=None,
    barcodes=None,
    control=False,
    filtered_only=False,
    quality=27,
    conversion='TC',
    snp_threshold=0.5,
    snp_csv=None,
    read_threshold=16,
    control_p_e=None,
    p_group_by=None,
    whitelist_path=None,
    n_threads=8,
    temp_dir=None,
    re=None,
    nasc=False,
    subset_threshold=8000,
    velocity=True,
    seed=None,
):

    def redo(key):
        return re in config.RE_CHOICES[:config.RE_CHOICES.index(key) + 1]

    STATS.start()
    stats_path = os.path.join(
        out_dir, f'{constants.STATS_PREFIX}_{dt.datetime.strftime(STATS.start_time, "%Y%m%d_%H%M%S_%f")}.json'
    )
    os.makedirs(out_dir, exist_ok=True)

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
    with STATS.step('sort', skipped=skip) or redo('sort'):
        if not skip:
            logger.info(f'Sorting {bam_path} with samtools to {sorted_bam_path}')
            pysam.sort(bam_path, '-o', sorted_bam_path, '-@', str(n_threads))
            bam_path = sorted_bam_path

    # Check if BAM index exists and create one if it doesn't.
    bai_path = f'{bam_path}.bai'
    skip = utils.all_exists([bai_path]) and not redo('index')
    with STATS.step('index', skipped=skip):
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
    with STATS.step('parse', skipped=skip):
        if not skip:
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
                barcodes=barcodes,
                n_threads=n_threads,
                temp_dir=temp_dir,
                nasc=nasc,
                velocity=velocity
            )
        else:
            logger.info('Skipped read and conversion parsing from BAM')

    # Detect SNPs
    snp_dir = os.path.join(out_dir, constants.SNP_DIR)
    coverage_path = os.path.join(snp_dir, constants.COVERAGE_FILENAME)
    coverage_index_path = os.path.join(snp_dir, constants.COVERAGE_INDEX_FILENAME)
    snps_path = os.path.join(snp_dir, constants.SNPS_FILENAME)
    snps_required = [coverage_path, coverage_index_path, snps_path]
    skip = not snp_threshold or (utils.all_exists(snps_required) and not redo('snp'))
    with STATS.step('snp', skipped=skip):
        if not skip:
            logger.info('Calculating coverage and detecting SNPs')
            os.makedirs(snp_dir, exist_ok=True)
            coverage_path, coverage_index_path = preprocessing.calculate_coverage(
                bam_path, {
                    contig: set(df_part['genome_i'])
                    for contig, df_part in preprocessing.read_conversions(
                        conversions_path, usecols=['contig', 'genome_i']
                    ).drop_duplicates().groupby('contig')
                },
                coverage_path,
                coverage_index_path,
                umi_tag=umi_tag,
                barcode_tag=barcode_tag,
                barcodes=barcodes,
                n_threads=n_threads,
                temp_dir=temp_dir,
                velocity=velocity
            )
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
            if not snp_threshold:
                logger.info('No SNP filtering will be done. Use `--snp-threshold` to detect possible SNPs.')
            else:
                logger.info('Skipped coverage calculation and SNP detection')

    # Count conversions and calculate mutation rates
    count_dir = os.path.join(out_dir, constants.COUNT_DIR)
    counts_path = os.path.join(count_dir, constants.COUNTS_FILENAME)
    count_required = [counts_path]
    skip = utils.all_exists(count_required) and not redo('count')
    with STATS.step('count', skipped=skip):
        if not skip:
            os.makedirs(count_dir, exist_ok=True)
            logger.info('Counting conversions')
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
                snps=snps,
                quality=quality,
                n_threads=n_threads,
                temp_dir=temp_dir
            )
        else:
            logger.info('Skipped conversion counting')

    aggregates_dir = os.path.join(out_dir, constants.AGGREGATES_DIR)
    rates_path = os.path.join(aggregates_dir, constants.RATES_FILENAME)
    aggregates_paths = {
        key: {
            conversion: os.path.join(aggregates_dir, key, f'{conversion}.csv')
            for conversion in preprocessing.conversion.CONVERSION_COLUMNS
        }
        for key in preprocessing.read_counts(counts_path, usecols=['velocity'])['velocity'].unique()
    } if velocity else {}
    aggregates_paths['transcriptome'] = {
        conversion: os.path.join(aggregates_dir, constants.TRANSCRIPTOME_DIR, f'{conversion}.csv')
        for conversion in preprocessing.conversion.CONVERSION_COLUMNS
    }
    aggregates_required = utils.flatten_dict_values(aggregates_paths) + [rates_path]
    df_counts_uncomplemented = None
    df_counts_complemented = None
    df_counts_transcriptome = None
    skip = utils.all_exists(aggregates_required) and not redo('aggregate')
    with STATS.step('aggregate', skipped=skip):
        if not skip:
            logger.info('Computing mutation rates')
            os.makedirs(aggregates_dir, exist_ok=True)
            df_counts_uncomplemented = preprocessing.read_counts(counts_path)
            df_counts_complemented = preprocessing.complement_counts(
                df_counts_uncomplemented, gene_infos or utils.read_pickle(genes_path)
            )
            rates_path = preprocessing.calculate_mutation_rates(
                df_counts_uncomplemented if nasc else df_counts_complemented, rates_path, group_by=p_group_by
            )

            df_counts_transcriptome = df_counts_complemented[df_counts_complemented['transcriptome']]
            dfs = {'transcriptome': df_counts_transcriptome}
            if velocity:
                dfs.update(preprocessing.split_counts_by_velocity(df_counts_complemented))
            for key, df in dfs.items():
                logger.info(f'Aggregating counts for `{key}` reads')
                velocity_aggregates_dir = os.path.join(aggregates_dir, key)
                os.makedirs(velocity_aggregates_dir, exist_ok=True)
                aggregates_paths[key] = preprocessing.aggregate_counts(df, velocity_aggregates_dir)
        else:
            logger.info('Skipped count aggregation')

    estimates_dir = os.path.join(out_dir, constants.ESTIMATES_DIR)
    p_e_path = os.path.join(estimates_dir, constants.P_E_FILENAME)
    p_c_path = os.path.join(estimates_dir, constants.P_C_FILENAME)
    estimates_paths = [p_e_path, p_c_path]
    value_columns = [conversion, conversion[0], 'count']
    skip = utils.all_exists(estimates_paths) and not redo('p')
    with STATS.step('p', skipped=skip):
        if not skip:
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
                logger.info('Estimating average mismatch rate in unlabeled RNA')
                if control:
                    p_e_path = estimation.estimate_p_e_control(
                        df_counts_complemented
                        if df_counts_complemented is not None else preprocessing.complement_counts(
                            preprocessing.read_counts(counts_path), gene_infos or utils.read_pickle(genes_path)
                        ),
                        p_e_path,
                        conversion=conversion,
                    )
                elif nasc:
                    p_e_path = estimation.estimate_p_e_nasc(
                        preprocessing.read_rates(rates_path),
                        p_e_path,
                        conversion=conversion,
                        group_by=p_group_by,
                    )
                else:
                    p_e_path = estimation.estimate_p_e(
                        df_counts_complemented
                        if df_counts_complemented is not None else preprocessing.complement_counts(
                            preprocessing.read_counts(counts_path), gene_infos or utils.read_pickle(genes_path)
                        ),
                        p_e_path,
                        conversion=conversion,
                        group_by=p_group_by,
                    )

            if not control:
                logger.info('Estimating average mismatch rate in labeled RNA')
                df_aggregates = preprocessing.merge_aggregates(
                    *[
                        preprocessing.read_aggregates(paths[conversion])
                        for key, paths in aggregates_paths.items()
                        if key != 'transcriptome'
                    ],
                    conversion=conversion
                ) if velocity else preprocessing.read_aggregates(aggregates_paths['transcriptome'][conversion])
                p_c_path = estimation.estimate_p_c(
                    df_aggregates,
                    estimation.read_p_e(p_e_path, group_by=p_group_by),
                    p_c_path,
                    group_by=p_group_by,
                    value_columns=value_columns,
                    n_threads=n_threads,
                )
            else:
                logger.info('Average mismatch rate in labeled RNA was not calculated because `--control` was provided')
        else:
            logger.info('Skipped rate estimation')

    if control:
        logger.info('Downstream processing skipped for controls')
        if snp_threshold:
            logger.info(f'Use `--snp-csv {snps_path}` to run test samples')
        logger.info(f'Use `--p-e {p_e_path}` for test samples')
        STATS.end()
        STATS.save(stats_path)
        return

    velocity_blacklist = ['unassigned', 'ambiguous']
    pi_paths = {
        key: os.path.join(estimates_dir, f'{key}.csv')
        for key in aggregates_paths.keys()
        if key not in velocity_blacklist
    }
    skip = utils.all_exists(list(pi_paths.values())) and not redo('pi')
    with STATS.step('pi', skipped=skip):
        if not skip:
            for key, paths in aggregates_paths.items():
                if key in velocity_blacklist:
                    continue
                logger.info(f'Estimating fraction of labeled RNA for `{key}` reads')
                pi_paths[key] = estimation.estimate_pi(
                    preprocessing.read_aggregates(paths[conversion]),
                    estimation.read_p_e(p_e_path, group_by=p_group_by),
                    estimation.read_p_c(p_c_path, group_by=p_group_by),
                    pi_paths[key],
                    p_group_by=p_group_by,
                    value_columns=value_columns,
                    n_threads=n_threads,
                    threshold=read_threshold,
                    subset_threshold=subset_threshold,
                    seed=seed,
                )
        else:
            logger.info('Skipped estimation of labeled RNA')

    adata_path = os.path.join(out_dir, constants.ADATA_FILENAME)
    skip = utils.all_exists([adata_path]) and not redo('split')
    with STATS.step('split', skipped=skip):
        if not skip:
            logger.info(f'Combining results into an Anndata object at {adata_path}')

            logger.info('Loading results for `transcriptome` reads')
            df_counts = (
                df_counts_uncomplemented
                if df_counts_uncomplemented is not None else preprocessing.read_counts(counts_path)
            )
            barcodes = barcodes or sorted(df_counts['barcode'].unique())
            features = sorted(df_counts['GX'].unique())

            # Deal with transcriptome reads first because those need to be
            # in adata.X
            df_counts_transcriptome = df_counts[df_counts['transcriptome']]
            matrix, matrix_unlabeled, matrix_labeled = preprocessing.split_counts(
                df_counts_transcriptome, barcodes, features, conversion=conversion
            )
            # Construct adata with umi counts as layers
            adata = anndata.AnnData(
                X=matrix,
                obs=pd.DataFrame(index=pd.Series(barcodes, name='barcodes')),
                var=pd.DataFrame(index=pd.Series(features, name='gene_id')),
                layers={
                    'X_unlabeled': matrix_unlabeled,
                    'X_labeled': matrix_labeled
                }
            )
            pis = estimation.read_pi(pi_paths['transcriptome'])
            (adata.layers['X_pi'], adata.layers['X_unlabeled_estimate'],
             adata.layers['X_labeled_estimate']) = estimation.split_matrix(adata.X, pis, barcodes, features)

            # All the other counts as layers
            for key in aggregates_paths.keys():
                if key != 'transcriptome' and key not in velocity_blacklist:
                    logger.info(f'Loading result for `{key}` reads')
                    df_counts_velocity = df_counts[df_counts['velocity'] == key]
                    # UMI counts
                    (adata.layers[key], adata.layers[f'{key}_unlabeled'],
                     adata.layers[f'{key}_labeled']) = preprocessing.split_counts(
                         df_counts_velocity, barcodes, features, conversion=conversion
                     )

                    # Estimates
                    pis = estimation.read_pi(pi_paths[key])
                    (
                        adata.layers[f'{key}_pi'], adata.layers[f'{key}_unlabeled_estimate'],
                        adata.layers[f'{key}_labeled_estimate']
                    ) = estimation.split_matrix(adata.layers[key], pis, barcodes, features)

            adata.write(adata_path, compression='gzip')
        else:
            logger.info('Skipped splitting of new and old RNA')
    STATS.end()
    STATS.save(stats_path)
