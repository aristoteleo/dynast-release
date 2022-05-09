# flake8: noqa
from .aggregation import (
    aggregate_counts,
    calculate_mutation_rates,
    merge_aggregates,
    read_aggregates,
    read_rates,
)
from .bam import (
    check_bam_contains_duplicate, check_bam_contains_secondary, check_bam_contains_unmapped, get_tags_from_bam,
    parse_all_reads, read_alignments, read_conversions, select_alignments, sort_and_index_bam
)
from .consensus import call_consensus
from .conversion import (
    complement_counts,
    CONVERSION_COMPLEMENT,
    count_conversions,
    deduplicate_counts,
    read_counts,
    split_counts_by_velocity,
)
from .coverage import calculate_coverage, read_coverage
from .snp import detect_snps, read_snp_csv, read_snps
