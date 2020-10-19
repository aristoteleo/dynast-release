import gzip
import os
import shutil
import tempfile
from unittest import TestCase


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
        cls.temp_dir = None
        cls.base_dir = os.path.dirname(os.path.abspath(__file__))
        cls.fixtures_dir = os.path.join(cls.base_dir, 'fixtures')
        cls.STAR_dir = os.path.join(cls.fixtures_dir, 'STAR')
        cls.barcodes_path = os.path.join(cls.STAR_dir, 'barcodes.tsv')
        cls.features_path = os.path.join(cls.STAR_dir, 'features.tsv')
        cls.matrix_path = os.path.join(cls.STAR_dir, 'matrix.mtx')

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

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)
