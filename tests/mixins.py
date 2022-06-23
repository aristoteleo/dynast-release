import gzip
import os
import shutil
import tempfile
from unittest import TestCase

import pandas as pd

from dynast.technology import TECHNOLOGIES_MAP


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
        cls.umi_technology = TECHNOLOGIES_MAP['dropseq']
        cls.smartseq_technology = TECHNOLOGIES_MAP['smartseq']

        # Paths
        cls.temp_dir = None
        cls.base_dir = os.path.dirname(os.path.abspath(__file__))
        cls.fixtures_dir = os.path.join(cls.base_dir, 'fixtures')
        cls.gtf_path = os.path.join(cls.fixtures_dir, 'not_sorted.gtf')
        cls.genes_path = os.path.join(cls.fixtures_dir, 'genes.pkl.gz')
        cls.transcripts_path = os.path.join(cls.fixtures_dir, 'transcripts.pkl.gz')
        cls.model_path = os.path.join(cls.fixtures_dir, 'model.pkl.gz')

        # Misc. fixtures
        cls.p_e_int_path = os.path.join(cls.fixtures_dir, 'p_e_int.csv')
        cls.p_c_int_path = os.path.join(cls.fixtures_dir, 'p_c_int.csv')
        cls.pi_int_path = os.path.join(cls.fixtures_dir, 'pi_int.csv')

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
        cls.umi_alignments_path = os.path.join(cls.umi_count_dir, 'alignments.csv')
        cls.umi_conversions_path = os.path.join(cls.umi_count_dir, 'conversions.csv')
        cls.umi_conversions_index_path = os.path.join(cls.umi_count_dir, 'conversions.idx')
        cls.umi_genes_path = cls.genes_path
        cls.umi_transcripts_path = cls.transcripts_path
        cls.umi_counts_path = os.path.join(cls.umi_count_dir, 'counts_TC.csv')
        cls.umi_rates_path = os.path.join(cls.umi_count_dir, 'rates.csv')
        cls.umi_adata_path = os.path.join(cls.umi_count_dir, 'adata.h5ad')

        # Estimate
        cls.umi_estimate_dir = os.path.join(cls.fixtures_dir, 'SRR11683995_estimate')
        cls.umi_aggregates_path = os.path.join(cls.umi_estimate_dir, 'A_total_TC.csv')
        cls.umi_p_e_path = os.path.join(cls.umi_estimate_dir, 'p_e.csv')
        cls.umi_p_c_path = os.path.join(cls.umi_estimate_dir, 'p_c_TC.csv')
        cls.umi_pi_path = os.path.join(cls.umi_estimate_dir, 'pi_total_TC.csv')

        # Estimate (alpha)
        cls.umi_alpha_estimate_dir = os.path.join(cls.fixtures_dir, 'SRR11683995_estimate_alpha')
        cls.umi_alpha_aggregates_path = os.path.join(cls.umi_alpha_estimate_dir, 'A_total_TC.csv')
        cls.umi_alpha_p_e_path = os.path.join(cls.umi_alpha_estimate_dir, 'p_e.csv')
        cls.umi_alpha_p_c_path = os.path.join(cls.umi_alpha_estimate_dir, 'p_c_TC.csv')
        cls.umi_alpha_pi_path = os.path.join(cls.umi_alpha_estimate_dir, 'pi_total_TC.csv')
        cls.umi_alpha_path = os.path.join(cls.umi_alpha_estimate_dir, 'alpha_total_TC.csv')

        #########
        # Control
        #########
        # Count
        cls.control_count_dir = os.path.join(cls.fixtures_dir, 'SRR11683995_count_control')
        cls.control_conversions_path = os.path.join(cls.control_count_dir, 'conversions.csv')
        cls.control_conversions_index_path = os.path.join(cls.control_count_dir, 'conversions.idx')
        cls.control_alignments_path = os.path.join(cls.control_count_dir, 'alignments.csv')
        cls.control_genes_path = cls.genes_path
        cls.control_transcripts_path = cls.transcripts_path
        cls.control_coverage_path = os.path.join(cls.control_count_dir, 'coverage.csv')
        cls.control_coverage_index_path = os.path.join(cls.control_count_dir, 'coverage.idx')
        cls.control_snps_path = os.path.join(cls.control_count_dir, 'snps.csv')
        cls.control_counts_path = os.path.join(cls.control_count_dir, 'counts_TC.csv')
        cls.control_rates_path = os.path.join(cls.control_count_dir, 'rates.csv')
        cls.control_p_e_path = os.path.join(cls.control_count_dir, 'p_e.csv')

        # Estimate
        cls.control_estimate_dir = os.path.join(cls.fixtures_dir, 'SRR11683995_estimate_control')
        cls.control_p_e_path = os.path.join(cls.control_estimate_dir, 'p_e.csv')

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
        cls.paired_conversions_path = os.path.join(cls.paired_count_dir, 'conversions.csv')
        cls.paired_conversions_index_path = os.path.join(cls.paired_count_dir, 'conversions.idx')
        cls.paired_alignments_path = os.path.join(cls.paired_count_dir, 'alignments.csv')
        cls.paired_genes_path = cls.genes_path
        cls.paired_transcripts_path = cls.transcripts_path
        cls.paired_counts_path = os.path.join(cls.paired_count_dir, 'counts_TC.csv')
        cls.paired_rates_path = os.path.join(cls.paired_count_dir, 'rates.csv')
        cls.paired_adata_path = os.path.join(cls.paired_count_dir, 'adata.h5ad')

        # Estimate
        cls.paired_estimate_dir = os.path.join(cls.fixtures_dir, 'smartseq_estimate')
        cls.paired_aggregates_path = os.path.join(cls.paired_estimate_dir, 'A_total_TC.csv')
        cls.paired_p_e_path = os.path.join(cls.paired_estimate_dir, 'p_e.csv')
        cls.paired_p_c_path = os.path.join(cls.paired_estimate_dir, 'p_c_TC.csv')
        cls.paired_pi_path = os.path.join(cls.paired_estimate_dir, 'pi_transcriptome_TC.csv')

        ######################################
        # Paired (smartseq, no velocity, nasc)
        ######################################
        # Align
        cls.nasc_align_dir = os.path.join(cls.fixtures_dir, 'nasc_align')
        cls.nasc_bam_path = os.path.join(cls.nasc_align_dir, 'Aligned.sortedByCoord.out.bam')

        # Count
        cls.nasc_count_dir = os.path.join(cls.fixtures_dir, 'nasc_count')
        cls.nasc_conversions_path = os.path.join(cls.nasc_count_dir, 'conversions.csv')
        cls.nasc_conversions_index_path = os.path.join(cls.nasc_count_dir, 'conversions.idx')
        cls.nasc_alignments_path = os.path.join(cls.nasc_count_dir, 'alignments.csv')
        cls.nasc_genes_path = cls.genes_path
        cls.nasc_transcripts_path = cls.transcripts_path
        cls.nasc_counts_path = os.path.join(cls.nasc_count_dir, 'counts_TC.csv')
        cls.nasc_rates_path = os.path.join(cls.nasc_count_dir, 'rates.csv')
        cls.nasc_adata_path = os.path.join(cls.nasc_count_dir, 'adata.h5ad')

        # Estimate
        cls.nasc_estimate_dir = os.path.join(cls.fixtures_dir, 'nasc_estimate')
        cls.nasc_aggregates_path = os.path.join(cls.nasc_estimate_dir, 'A_transcriptome_TC.csv')
        cls.nasc_p_e_path = os.path.join(cls.nasc_estimate_dir, 'p_e.csv')
        cls.nasc_p_c_path = os.path.join(cls.nasc_estimate_dir, 'p_c_TC.csv')
        cls.nasc_pi_path = os.path.join(cls.nasc_estimate_dir, 'pi_transcriptome_TC.csv')

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)
