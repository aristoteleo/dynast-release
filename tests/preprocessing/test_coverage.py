import os
from unittest import mock, TestCase

import pandas as pd

import dynast.preprocessing.coverage as coverage
import dynast.utils as utils
from tests import mixins


class TestCoverage(mixins.TestMixin, TestCase):

    def test_calculate_coverage(self):
        coverage_path = os.path.join(self.temp_dir, 'coverage.csv')
        index_path = os.path.join(self.temp_dir, 'coverage.idx')
        conversions = {
            contig: set(df_part['genome_i'])
            for contig, df_part in pd.read_csv(self.umi_conversions_path, usecols=['contig', 'genome_i']
                                               ).drop_duplicates().groupby('contig')
        }
        with mock.patch('dynast.preprocessing.coverage.utils.display_progress_with_counter'):
            self.assertEqual((coverage_path, index_path),
                             coverage.calculate_coverage(
                                 self.umi_bam_path,
                                 conversions,
                                 coverage_path,
                                 index_path,
                                 umi_tag='UB',
                                 barcode_tag='CB',
                                 barcodes=None,
                                 n_threads=2,
                                 temp_dir=self.temp_dir,
                                 velocity=True,
                             ))
            self.assertTrue(mixins.files_equal(self.umi_coverage_path, coverage_path))
            self.assertEqual(utils.read_pickle(self.umi_coverage_index_path), utils.read_pickle(index_path))
