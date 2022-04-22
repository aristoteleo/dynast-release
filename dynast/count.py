import datetime as dt
import os

import ngs_tools as ngs
import pysam

from . import config, constants, preprocessing, utils
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
    conversions=frozenset([('TC',)]),
    snp_threshold=0.5,
    snp_min_coverage=1,
    snp_csv=None,
    n_threads=8,
    temp_dir=None,
    velocity=True,
    nasc=False,
    overwrite=False,
):
    stats = Stats()
    stats.start()
    stats_path = os.path.join(
        out_dir, f'{constants.STATS_PREFIX}_{dt.datetime.strftime(stats.start_time, "%Y%m%d_%H%M%S_%f")}.json'
    )
    os.makedirs(out_dir, exist_ok=True)

    all_conversions = sorted(utils.flatten_iter(conversions))

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
    if not bam_sorted:
        logger.info(f'Sorting {bam_path} with samtools to {sorted_bam_path}')
        pysam.sort(bam_path, '-o', sorted_bam_path, '-@', str(n_threads))
        bam_path = sorted_bam_path

    # Check if BAM index exists and create one if it doesn't.
    bai_path = f'{bam_path}.bai'
    if not utils.all_exists(bai_path):
        logger.info(f'Indexing {bam_path} with samtools to {bai_path}')
        pysam.index(bam_path, bai_path, '-@', str(n_threads))

    # Parse BAM and save results
    conversions_path = os.path.join(out_dir, constants.CONVERSIONS_FILENAME)
    index_path = os.path.join(out_dir, constants.CONVERSIONS_INDEX_FILENAME)
    alignments_path = os.path.join(out_dir, constants.ALIGNMENTS_FILENAME)
    genes_path = os.path.join(out_dir, constants.GENES_FILENAME)
    conversions_required = [conversions_path, index_path, alignments_path, genes_path]
    bam_splits = None
    bam_parsed = False
    if not utils.all_exists(*conversions_required) or overwrite:
        logger.info('Parsing gene and transcript information from GTF')
        gene_infos, transcript_infos = ngs.gtf.genes_and_transcripts_from_gtf(gtf_path, use_version=False)
        utils.write_pickle(gene_infos, genes_path)

        logger.info(f'Parsing read conversion information from BAM to {conversions_path}')
        conversions_path, alignments_path, index_path, bam_splits = preprocessing.parse_all_reads(
            bam_path,
            conversions_path,
            alignments_path,
            index_path,
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
            velocity=velocity,
            return_splits=True
        )
        bam_parsed = True
    else:
        logger.warning('Skipped BAM parsing because files already exist. Use `--overwrite` to re-parse the BAM.')
        gene_infos = utils.read_pickle(genes_path)

    # Check consistency of alignments
    small_alignments = preprocessing.read_alignments(alignments_path, nrows=config.BAM_PEEK_READS)
    barcode_is_na = small_alignments['barcode'] == 'NA'
    umi_is_na = small_alignments['umi'] == 'NA'
    if barcode_tag and barcode_is_na.any():
        raise Exception(
            "`--barcode-tag` was provided but existing files contain NA barcodes. "
            'Re-run `dynast count` with `--overwrite` to fix this inconsistency.'
        )
    elif not barcode_tag and (~barcode_is_na).any():
        raise Exception(
            "`--barcode-tag` was not provided but existing files contain barcodes. "
            'Re-run `dynast count` with `--overwrite` to fix this inconsistency.'
        )
    if umi_tag and umi_is_na.any():
        raise Exception(
            "`--umi-tag` was provided but existing files contain NA UMIs. "
            'Re-run `dynast count` with `--overwrite` to fix this inconsistency.'
        )
    elif not umi_tag and (~umi_is_na).any():
        raise Exception(
            "`--umi-tag` was not provided but existing files contain UMIs. "
            'Re-run `dynast count` with `--overwrite` to fix this inconsistency.'
        )

    # Save conversions
    redo_snp = False
    convs_path = os.path.join(out_dir, constants.CONVS_FILENAME)
    if utils.all_exists(convs_path):
        prev_conversions = utils.read_pickle(convs_path)
        if conversions != prev_conversions:
            logger.warning(f'Conversions changed from {prev_conversions} in previous run to {conversions}.')
            redo_snp = True
    else:
        redo_snp = True

    # Detect SNPs
    coverage_path = os.path.join(out_dir, constants.COVERAGE_FILENAME)
    snps_path = os.path.join(out_dir, constants.SNPS_FILENAME)
    snp_required = [convs_path, coverage_path, snps_path]
    if snp_threshold:
        if not utils.all_exists(*snp_required) or redo_snp or bam_parsed:
            logger.info('Selecting alignments to use for SNP detection')
            alignments = preprocessing.select_alignments(preprocessing.read_alignments(alignments_path))

            snp_conversions = set(
                all_conversions + [preprocessing.CONVERSION_COMPLEMENT[conv] for conv in all_conversions]
            )
            logger.info(f'Selecting genomic locations with {snp_conversions} conversions in forward strand.')
            df_conversions = preprocessing.read_conversions(conversions_path)

            # Subset to selected alignments.
            df_conversions = df_conversions[[
                key in alignments for key in df_conversions[['read_id', 'index']].itertuples(index=False, name=None)
            ]]

            # Subset to conversions of interest
            mask = None
            for conv in snp_conversions:
                _mask = (df_conversions['original'] == conv[0]) & (df_conversions['converted'] == conv[1])
                if mask is None:
                    mask = _mask
                else:
                    mask |= _mask
            df_conversions = df_conversions.loc[mask, ['contig', 'genome_i']]

            logger.info(f'Calculating coverage and outputting to {coverage_path}')
            coverage_path = preprocessing.calculate_coverage(
                bam_path, {
                    contig: set(df_part['genome_i'])
                    for contig, df_part in df_conversions.groupby('contig', sort=False, observed=True)
                },
                coverage_path,
                alignments=alignments,
                umi_tag=umi_tag,
                barcode_tag=barcode_tag,
                gene_tag=gene_tag,
                barcodes=barcodes,
                n_threads=n_threads,
                temp_dir=temp_dir,
                velocity=velocity,
                splits=bam_splits
            )
            coverage = preprocessing.read_coverage(coverage_path)

            logger.info(f'Detecting SNPs with threshold {snp_threshold} to {snps_path}')
            snps_path = preprocessing.detect_snps(
                conversions_path,
                index_path,
                coverage,
                snps_path,
                alignments=alignments,
                conversions=snp_conversions,
                quality=quality,
                threshold=snp_threshold,
                min_coverage=snp_min_coverage,
                n_threads=n_threads,
            )

            utils.write_pickle(conversions, convs_path)
        else:
            logger.warning(
                'Skipped SNP detection because files already exist. '
                f'Remove {convs_path} to run SNP detection again.'
            )
    else:
        utils.write_pickle(conversions, convs_path)

    # Count conversions and calculate mutation rates
    counts_path = os.path.join(out_dir, f'{constants.COUNTS_PREFIX}_{"_".join(all_conversions)}.csv')
    logger.info(f'Counting conversions to {counts_path}')
    snps = utils.merge_dictionaries(
        preprocessing.read_snps(snps_path) if snp_threshold else {},
        preprocessing.read_snp_csv(snp_csv) if snp_csv else {},
        f=set.union,
        default=set,
    )
    counts_path = preprocessing.count_conversions(
        conversions_path,
        alignments_path,
        index_path,
        counts_path,
        gene_infos,
        barcodes=barcodes,
        snps=snps,
        quality=quality,
        conversions=all_conversions,
        n_threads=n_threads,
        temp_dir=temp_dir
    )
    df_counts_uncomplemented = preprocessing.read_counts(counts_path)
    df_counts_complemented = preprocessing.complement_counts(df_counts_uncomplemented, gene_infos)

    if barcodes:
        count_barcodes = set(df_counts_uncomplemented['barcode'])
        missing_barcodes = set(barcodes) - count_barcodes
        if missing_barcodes:
            logger.warning(
                f'{len(missing_barcodes)} barcodes are missing from {counts_path}. '
                'Re-run `dynast count` with `--overwrite` to fix this inconsistency. '
                'Otherwise, all missing barcodes will be ignored. '
            )

    # Calculate mutation rates
    rates_path = os.path.join(out_dir, constants.RATES_FILENAME)
    rates_path = preprocessing.calculate_mutation_rates(
        df_counts_uncomplemented if nasc else df_counts_complemented, rates_path, group_by=['barcode']
    )

    if control:
        logger.info('Downstream processing skipped for controls')
        if snp_threshold:
            logger.info(f'Use `--snp-csv {snps_path}` to run test samples')
    else:
        adata_path = os.path.join(out_dir, constants.ADATA_FILENAME)
        logger.info(f'Combining results into Anndata object at {adata_path}')
        adata = utils.results_to_adata(df_counts_complemented, conversions, gene_infos=gene_infos)
        adata.write(adata_path)
    stats.end()
    stats.save(stats_path)
