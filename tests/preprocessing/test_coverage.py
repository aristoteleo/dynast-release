import os
from unittest import mock, TestCase

import pandas as pd

import dynast.preprocessing.bam as bam
import dynast.preprocessing.coverage as coverage

from .. import mixins


class TestCoverage(mixins.TestMixin, TestCase):

    def test_read_coverage(self):
        coverage.read_coverage(self.control_coverage_path)

    def test_calculate_coverage(self):
        coverage_path = os.path.join(self.temp_dir, 'coverage.csv')
        conversions = {
            contig: set(df_part['genome_i'])
            for contig, df_part in pd.read_csv(self.control_conversions_path, usecols=['contig', 'genome_i']
                                               ).drop_duplicates().groupby('contig')
        }
        alignments = bam.select_alignments(bam.read_alignments(self.control_alignments_path))
        with mock.patch('dynast.preprocessing.coverage.utils.display_progress_with_counter'):
            self.assertEqual(
                coverage_path,
                coverage.calculate_coverage(
                    self.umi_bam_path,
                    conversions,
                    coverage_path,
                    alignments=alignments,
                    umi_tag='UB',
                    barcode_tag='CB',
                    barcodes=None,
                    n_threads=2,
                    temp_dir=self.temp_dir,
                    velocity=True,
                )
            )
            self.assertTrue(mixins.files_equal(self.control_coverage_path, coverage_path))
