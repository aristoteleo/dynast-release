import gzip
import os
import shutil
import tempfile
from unittest import TestCase

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

        # FASTQS
        cls.fastqs = [
            os.path.join(cls.fixtures_dir, 'SRR11683995_1.fastq.gz'),
            os.path.join(cls.fixtures_dir, 'SRR11683995_2.fastq.gz')
        ]

        # Align
        cls.align_dir = os.path.join(cls.fixtures_dir, 'SRR11683995_align')
        cls.bam_path = os.path.join(cls.align_dir, 'Aligned.sortedByCoord.out.bam')

        # Count
        cls.count_dir = os.path.join(cls.fixtures_dir, 'SRR11683995_count')
        cls.count_parse_dir = os.path.join(cls.count_dir, '0_parse')
        cls.count_count_dir = os.path.join(cls.count_dir, '1_count')
        cls.count_aggregate_dir = os.path.join(cls.count_dir, '2_aggregate')
        cls.count_estimate_dir = os.path.join(cls.count_dir, '3_estimate')
        cls.adata_path = os.path.join(cls.count_dir, 'adata.h5ad')

        cls.conversions_path = os.path.join(cls.count_parse_dir, 'conversions.csv')
        cls.conversions_index_path = os.path.join(cls.count_parse_dir, 'conversions.idx')
        cls.no_conversions_path = os.path.join(cls.count_parse_dir, 'no_conversions.csv')
        cls.no_conversions_index_path = os.path.join(cls.count_parse_dir, 'no_conversions.idx')
        cls.genes_path = os.path.join(cls.count_parse_dir, 'genes.pkl.gz')
        cls.transcripts_path = os.path.join(cls.count_parse_dir, 'transcripts.pkl.gz')
        cls.count_counts_path = os.path.join(cls.count_count_dir, 'counts.csv')
        cls.rates_path = os.path.join(cls.count_aggregate_dir, 'rates.csv')
        cls.aggregates_paths = {
            key: {
                conversion: os.path.join(cls.count_aggregate_dir, key, f'{conversion}.csv')
                for conversion in cls.conversions
            }
            for key in cls.types
        }
        cls.p_e_path = os.path.join(cls.count_estimate_dir, 'p_e.csv')
        cls.p_c_path = os.path.join(cls.count_estimate_dir, 'p_c.csv')
        cls.pi_paths = {key: os.path.join(cls.count_estimate_dir, f'{key}.csv') for key in cls.types}

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)
