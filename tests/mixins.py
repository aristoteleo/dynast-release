import gzip
import os
import shutil
import tempfile
from unittest import TestCase


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

        # Paths
        cls.temp_dir = None
        cls.base_dir = os.path.dirname(os.path.abspath(__file__))
        cls.fixtures_dir = os.path.join(cls.base_dir, 'fixtures')
        cls.STAR_dir = os.path.join(cls.fixtures_dir, 'STAR')
        cls.STAR_barcodes_path = os.path.join(cls.STAR_dir, 'barcodes.tsv')
        cls.STAR_features_path = os.path.join(cls.STAR_dir, 'features.tsv')
        cls.STAR_matrix_path = os.path.join(cls.STAR_dir, 'matrix.mtx')

        # BAM
        cls.bam_path = os.path.join(cls.STAR_dir, 'small.bam')
        cls.bam_conversions_paths = [
            os.path.join(cls.STAR_dir, 'conversions1.csv'),
            os.path.join(cls.STAR_dir, 'conversions2.csv'),
        ]
        cls.bam_index_paths = [
            os.path.join(cls.STAR_dir, 'conversions1.idx'),
            os.path.join(cls.STAR_dir, 'conversions2.idx'),
        ]
        cls.bam_coverage_paths = [
            os.path.join(cls.STAR_dir, 'coverage1.csv'),
            os.path.join(cls.STAR_dir, 'coverage2.csv'),
        ]
        cls.bam_conversions_path = os.path.join(cls.STAR_dir, 'conversions.csv')
        cls.bam_index_path = os.path.join(cls.STAR_dir, 'conversions.idx')
        cls.bam_coverage_path = os.path.join(cls.STAR_dir, 'coverage.csv')

        # Conversion
        cls.conversions_path = os.path.join(cls.fixtures_dir, 'conversions.csv')
        cls.index_path = os.path.join(cls.fixtures_dir, 'conversions.idx')
        cls.counts_part_path = os.path.join(cls.fixtures_dir, 'counts_part.csv')
        cls.counts_part_paths = [
            os.path.join(cls.fixtures_dir, 'counts_part1.csv'),
            os.path.join(cls.fixtures_dir, 'counts_part2.csv'),
        ]
        cls.barcodes_path = os.path.join(cls.fixtures_dir, 'barcodes.csv')
        cls.genes_path = os.path.join(cls.fixtures_dir, 'genes.csv')
        cls.counts_path = os.path.join(cls.fixtures_dir, 'counts.csv')

        # Aggregation
        cls.aggregates_dir = os.path.join(cls.fixtures_dir, 'aggregates')
        cls.rates_barcode_path = os.path.join(cls.fixtures_dir, 'rates_barcode.csv')
        cls.rates_gene_path = os.path.join(cls.fixtures_dir, 'rates_gene.csv')
        cls.rates_none_path = os.path.join(cls.fixtures_dir, 'rates_none.csv')
        cls.aggregate_paths = {
            conversion: os.path.join(cls.aggregates_dir, f'{conversion}.csv')
            for conversion in cls.conversions
        }

        # Estimation
        cls.estimates_dir = os.path.join(cls.fixtures_dir, 'estimates')
        cls.p_e_barcode_path = os.path.join(cls.estimates_dir, 'p_e_barcode.csv')
        cls.p_e_gene_path = os.path.join(cls.estimates_dir, 'p_e_gene.csv')
        cls.p_e_none_path = os.path.join(cls.estimates_dir, 'p_e_none.csv')

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)
