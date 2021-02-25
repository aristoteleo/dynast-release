import gzip
import os
import shutil
import tempfile
from unittest import TestCase

import pandas as pd

from dynast.technology import Technology


def tqdm_mock(iterable, *args, **kwargs):
    return iterable


def import_mock(mocked, *args):

    def _import_mock(name, *args):
        return mocked

    return _import_mock


def files_equal(file1, file2, gzipped=False):
    open_f = gzip.open if gzipped else open
    with open_f(file1, 'r') as f1, open_f(file2, 'r') as f2:
        return f1.read() == f2.read()


def dataframes_equal(file1, file2, sort_by=None):
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    if sort_by:
        df1 = df1.sort_values(sort_by).reset_index(drop=True)
        df2 = df2.sort_values(sort_by).reset_index(drop=True)
    return df1.equals(df2)


class TestMixin(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.conversions = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG']
        cls.bases = ['A', 'C', 'G', 'T']
        cls.types = ['transcriptome', 'spliced', 'unspliced', 'ambiguous']
        cls.umi_technology = Technology('umi_technology', {'--arg1': 'value1', '--arg2': 2}, None)
        cls.smartseq_technology = Technology('smartseq', {'--arg1': 'value1', '--arg2': 2}, None)

        # Paths
        cls.temp_dir = None
        cls.base_dir = os.path.dirname(os.path.abspath(__file__))
        cls.fixtures_dir = os.path.join(cls.base_dir, 'fixtures')
        cls.gtf_path = os.path.join(cls.fixtures_dir, 'not_sorted.gtf')

        ###########################
        # UMI-based (with velocity)
        ###########################
        # FASTQS
        cls.umi_fastqs = [
            os.path.join(cls.fixtures_dir, 'SRR11683995_1.fastq.gz'),
            os.path.join(cls.fixtures_dir, 'SRR11683995_2.fastq.gz')
        ]

        # Align
        cls.umi_align_dir = os.path.join(cls.fixtures_dir, 'SRR11683995_align')
        cls.umi_bam_path = os.path.join(cls.umi_align_dir, 'Aligned.sortedByCoord.out.bam')

        # Count
        cls.umi_count_dir = os.path.join(cls.fixtures_dir, 'SRR11683995_count')
        cls.umi_count_parse_dir = os.path.join(cls.umi_count_dir, '0_parse')
        cls.umi_count_count_dir = os.path.join(cls.umi_count_dir, '1_count')
        cls.umi_count_aggregate_dir = os.path.join(cls.umi_count_dir, '2_aggregate')
        cls.umi_count_estimate_dir = os.path.join(cls.umi_count_dir, '3_estimate')
        cls.umi_adata_path = os.path.join(cls.umi_count_dir, 'adata.h5ad')

        cls.umi_conversions_path = os.path.join(cls.umi_count_parse_dir, 'conversions.csv')
        cls.umi_conversions_index_path = os.path.join(cls.umi_count_parse_dir, 'conversions.idx')
        cls.umi_no_conversions_path = os.path.join(cls.umi_count_parse_dir, 'no_conversions.csv')
        cls.umi_no_conversions_index_path = os.path.join(cls.umi_count_parse_dir, 'no_conversions.idx')
        cls.umi_genes_path = os.path.join(cls.umi_count_parse_dir, 'genes.pkl.gz')
        cls.umi_transcripts_path = os.path.join(cls.umi_count_parse_dir, 'transcripts.pkl.gz')
        cls.umi_counts_path = os.path.join(cls.umi_count_count_dir, 'counts_TC.csv')
        cls.umi_rates_path = os.path.join(cls.umi_count_aggregate_dir, 'rates.csv')
        cls.umi_aggregates_path = os.path.join(cls.umi_count_aggregate_dir, 'total_TC.csv')
        cls.umi_p_e_path = os.path.join(cls.umi_count_estimate_dir, 'p_e.csv')
        cls.umi_p_c_path = os.path.join(cls.umi_count_estimate_dir, 'p_c_TC.csv')
        cls.umi_pi_path = os.path.join(cls.umi_count_estimate_dir, 'total_TC.csv')

        #########
        # Control
        #########
        # Count
        cls.control_count_dir = os.path.join(cls.fixtures_dir, 'SRR11683995_count_control')
        cls.control_count_parse_dir = os.path.join(cls.control_count_dir, '0_parse')
        cls.control_count_snp_dir = os.path.join(cls.control_count_dir, '0_snp')
        cls.control_count_count_dir = os.path.join(cls.control_count_dir, '1_count')
        cls.control_count_aggregate_dir = os.path.join(cls.control_count_dir, '2_aggregate')
        cls.control_count_estimate_dir = os.path.join(cls.control_count_dir, '3_estimate')
        cls.control_adata_path = os.path.join(cls.control_count_dir, 'adata.h5ad')

        cls.control_conversions_path = os.path.join(cls.control_count_parse_dir, 'conversions.csv')
        cls.control_conversions_index_path = os.path.join(cls.control_count_parse_dir, 'conversions.idx')
        cls.control_no_conversions_path = os.path.join(cls.control_count_parse_dir, 'no_conversions.csv')
        cls.control_no_conversions_index_path = os.path.join(cls.control_count_parse_dir, 'no_conversions.idx')
        cls.control_genes_path = os.path.join(cls.control_count_parse_dir, 'genes.pkl.gz')
        cls.control_transcripts_path = os.path.join(cls.control_count_parse_dir, 'transcripts.pkl.gz')
        cls.control_coverage_path = os.path.join(cls.control_count_snp_dir, 'coverage.csv')
        cls.control_coverage_index_path = os.path.join(cls.control_count_snp_dir, 'coverage.idx')
        cls.control_snps_path = os.path.join(cls.control_count_snp_dir, 'snps.csv')
        cls.control_counts_path = os.path.join(cls.control_count_count_dir, 'counts_TC.csv')
        cls.control_rates_path = os.path.join(cls.control_count_aggregate_dir, 'rates.csv')
        cls.control_p_e_path = os.path.join(cls.control_count_estimate_dir, 'p_e.csv')

        ################################
        # Paired (smartseq, no velocity)
        ################################
        cls.paired_fastqs = [
            os.path.join(cls.fixtures_dir, 'P1_A06_S6_R1_001.fastq.gz'),
            os.path.join(cls.fixtures_dir, 'P1_A06_S6_R2_001.fastq.gz'),
            os.path.join(cls.fixtures_dir, 'P1_A09_S9_R1_001.fastq.gz'),
            os.path.join(cls.fixtures_dir, 'P1_A09_S9_R2_001.fastq.gz'),
        ]

        # Align
        cls.paired_align_dir = os.path.join(cls.fixtures_dir, 'smartseq_align')
        cls.paired_bam_path = os.path.join(cls.paired_align_dir, 'Aligned.sortedByCoord.out.bam')

        # Count
        cls.paired_count_dir = os.path.join(cls.fixtures_dir, 'smartseq_count')
        cls.paired_count_parse_dir = os.path.join(cls.paired_count_dir, '0_parse')
        cls.paired_count_count_dir = os.path.join(cls.paired_count_dir, '1_count')
        cls.paired_count_aggregate_dir = os.path.join(cls.paired_count_dir, '2_aggregate')
        cls.paired_count_estimate_dir = os.path.join(cls.paired_count_dir, '3_estimate')
        cls.paired_adata_path = os.path.join(cls.paired_count_dir, 'adata.h5ad')

        cls.paired_conversions_path = os.path.join(cls.paired_count_parse_dir, 'conversions.csv')
        cls.paired_conversions_index_path = os.path.join(cls.paired_count_parse_dir, 'conversions.idx')
        cls.paired_no_conversions_path = os.path.join(cls.paired_count_parse_dir, 'no_conversions.csv')
        cls.paired_no_conversions_index_path = os.path.join(cls.paired_count_parse_dir, 'no_conversions.idx')
        cls.paired_genes_path = os.path.join(cls.paired_count_parse_dir, 'genes.pkl.gz')
        cls.paired_transcripts_path = os.path.join(cls.paired_count_parse_dir, 'transcripts.pkl.gz')
        cls.paired_counts_path = os.path.join(cls.paired_count_count_dir, 'counts_TC.csv')
        cls.paired_rates_path = os.path.join(cls.paired_count_aggregate_dir, 'rates.csv')
        cls.paired_aggregates_path = os.path.join(cls.paired_count_aggregate_dir, 'transcriptome_TC.csv')
        cls.paired_p_e_path = os.path.join(cls.paired_count_estimate_dir, 'p_e.csv')
        cls.paired_p_c_path = os.path.join(cls.paired_count_estimate_dir, 'p_c.csv')
        cls.paired_pi_path = os.path.join(cls.paired_count_estimate_dir, 'transcriptome_TC.csv')

        ######################################
        # Paired (smartseq, no velocity, nasc)
        ######################################
        # Align
        cls.nasc_align_dir = os.path.join(cls.fixtures_dir, 'nasc_align')
        cls.nasc_bam_path = os.path.join(cls.nasc_align_dir, 'Aligned.sortedByCoord.out.bam')

        # Count
        cls.nasc_count_dir = os.path.join(cls.fixtures_dir, 'nasc_count')
        cls.nasc_count_parse_dir = os.path.join(cls.nasc_count_dir, '0_parse')
        cls.nasc_count_count_dir = os.path.join(cls.nasc_count_dir, '1_count')
        cls.nasc_count_aggregate_dir = os.path.join(cls.nasc_count_dir, '2_aggregate')
        cls.nasc_count_estimate_dir = os.path.join(cls.nasc_count_dir, '3_estimate')
        cls.nasc_adata_path = os.path.join(cls.nasc_count_dir, 'adata.h5ad')

        cls.nasc_conversions_path = os.path.join(cls.nasc_count_parse_dir, 'conversions.csv')
        cls.nasc_conversions_index_path = os.path.join(cls.nasc_count_parse_dir, 'conversions.idx')
        cls.nasc_no_conversions_path = os.path.join(cls.nasc_count_parse_dir, 'no_conversions.csv')
        cls.nasc_no_conversions_index_path = os.path.join(cls.nasc_count_parse_dir, 'no_conversions.idx')
        cls.nasc_genes_path = os.path.join(cls.nasc_count_parse_dir, 'genes.pkl.gz')
        cls.nasc_transcripts_path = os.path.join(cls.nasc_count_parse_dir, 'transcripts.pkl.gz')
        cls.nasc_counts_path = os.path.join(cls.nasc_count_count_dir, 'counts_TC.csv')
        cls.nasc_rates_path = os.path.join(cls.nasc_count_aggregate_dir, 'rates.csv')
        cls.nasc_aggregates_path = os.path.join(cls.nasc_count_aggregate_dir, 'transcriptome_TC.csv')
        cls.nasc_p_e_path = os.path.join(cls.nasc_count_estimate_dir, 'p_e.csv')
        cls.nasc_p_c_path = os.path.join(cls.nasc_count_estimate_dir, 'p_c_TC.csv')
        cls.nasc_pi_path = os.path.join(cls.nasc_count_estimate_dir, 'transcriptome_TC.csv')

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)
