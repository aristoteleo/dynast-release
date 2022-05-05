import datetime as dt
import os

import ngs_tools as ngs

from . import config, constants, preprocessing
from .logging import logger
from .stats import Stats


@logger.namespaced('consensus')
def consensus(
    bam_path,
    gtf_path,
    out_dir,
    strand='forward',
    umi_tag=None,
    barcode_tag=None,
    gene_tag='GX',
    quality=27,
    add_RS_RI=False,
    n_threads=8,
    temp_dir=None,
):
    stats = Stats()
    stats.start()
    stats_path = os.path.join(
        out_dir, f'{constants.STATS_PREFIX}_{dt.datetime.strftime(stats.start_time, "%Y%m%d_%H%M%S_%f")}.json'
    )
    os.makedirs(out_dir, exist_ok=True)

    # Sort and index bam.
    bam_path = preprocessing.sort_and_index_bam(
        bam_path, '{}.sortedByCoord{}'.format(*os.path.splitext(bam_path)), n_threads=n_threads
    )

    # Check BAM tags
    tags = preprocessing.get_tags_from_bam(bam_path, config.BAM_PEEK_READS, n_threads=n_threads)
    required_tags = config.BAM_REQUIRED_TAGS.copy()
    if barcode_tag:
        required_tags.append(barcode_tag)
    elif config.BAM_BARCODE_TAG in tags:
        logger.warning(
            f'BAM contains reads with {config.BAM_BARCODE_TAG} tag. Are you sure '
            f'you didn\'t mean to provide `--barcode-tag {config.BAM_BARCODE_TAG}`?'
        )
    elif config.BAM_READGROUP_TAG in tags:
        logger.warning(
            f'BAM contains reads with {config.BAM_READGROUP_TAG} tag. Are you sure '
            f'you didn\'t mean to provide `--barcode-tag {config.BAM_READGROUP_TAG}`?'
        )
    if umi_tag:
        required_tags.append(umi_tag)
    elif config.BAM_UMI_TAG in tags:
        logger.warning(
            f'BAM contains reads with {config.BAM_UMI_TAG} tag. Are you sure '
            f'you didn\'t mean to provide `--umi-tag {config.BAM_UMI_TAG}`?'
        )
    if gene_tag:
        required_tags.append(gene_tag)
    elif config.BAM_GENE_TAG in tags:
        logger.warning(
            f'BAM contains reads with {config.BAM_GENE_TAG} tag. Are you sure '
            f'you didn\'t mean to provide `--gene-tag {config.BAM_GENE_TAG}`?'
        )
    missing_tags = set(required_tags) - tags
    if missing_tags:
        raise Exception(
            f'First {config.BAM_PEEK_READS} reads in the BAM do not contain the following required tags: '
            f'{", ".join(missing_tags)}. '
        )

    # Check BAM alignments
    if preprocessing.check_bam_contains_secondary(bam_path, config.BAM_PEEK_READS, n_threads=n_threads):
        logger.warning(
            'BAM contains secondary alignments, which will be ignored. Only primary '
            'alignments are considered.'
        )
    if preprocessing.check_bam_contains_unmapped(bam_path):
        logger.warning('BAM contains unmapped reads, which will be ignored.')
    if preprocessing.check_bam_contains_duplicate(bam_path, config.BAM_PEEK_READS, n_threads=n_threads):
        logger.warning('BAM contains duplicate reads, which will be ignored.')

    logger.info('Parsing gene and transcript information from GTF')
    gene_infos, transcript_infos = ngs.gtf.genes_and_transcripts_from_gtf(gtf_path, use_version=False)

    consensus_path = os.path.join(out_dir, constants.CONSENSUS_BAM_FILENAME)
    logger.info(f'Calling consensus sequences from BAM to {consensus_path}')
    preprocessing.call_consensus(
        bam_path,
        consensus_path,
        gene_infos,
        strand=strand,
        umi_tag=umi_tag,
        barcode_tag=barcode_tag,
        gene_tag=gene_tag,
        quality=quality,
        add_RS_RI=add_RS_RI,
        temp_dir=temp_dir,
        n_threads=n_threads
    )
    stats.end()
    stats.save(stats_path)
