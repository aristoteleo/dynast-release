import os
from unittest import mock, TestCase

import pandas as pd

import dynast.preprocessing.snp as snp

from .. import mixins


class TestSnp(mixins.TestMixin, TestCase):

    def test_read_snps(self):
        snp.read_snps(self.umi_snps_path)

    def test_read_snp_csv(self):
        snp.read_snp_csv(self.umi_snps_path)

    def test_extract_conversions(self):
        with mock.patch('dynast.preprocessing.snp.utils.display_progress_with_counter'):
            conversions = snp.extract_conversions(
                self.umi_conversions_path, self.umi_conversions_index_path, quality=35, n_threads=2
            )

        df = pd.read_csv(self.umi_conversions_path)
        conversions_truth = {}
        for (contig, genome_i), count in dict(df[df['quality'] > 35].groupby(['contig', 'genome_i']).size()).items():
            conversions_truth.setdefault(contig,
                                         {})[genome_i] = conversions_truth.setdefault(contig, {}).get(genome_i, 0) + 1
        self.assertEqual(conversions_truth, conversions)

    def test_extract_coverage(self):
        with mock.patch('dynast.preprocessing.snp.utils.display_progress_with_counter'):
            coverage = snp.extract_coverage(self.umi_coverage_path, self.umi_coverage_index_path, n_threads=2)

        df = pd.read_csv(self.umi_coverage_path)
        coverage_truth = {}
        for (contig, genome_i), count in dict(df.groupby(['contig', 'genome_i']).sum()['coverage']).items():
            coverage_truth.setdefault(contig, {})[genome_i] = count
        self.assertEqual(coverage_truth, coverage)

    def test_detect_snps(self):
        snps_path = os.path.join(self.temp_dir, 'snps.csv')
        with mock.patch('dynast.preprocessing.snp.utils.display_progress_with_counter'):
            self.assertEqual(
                snps_path,
                snp.detect_snps(
                    self.umi_conversions_path,
                    self.umi_conversions_index_path,
                    self.umi_coverage_path,
                    self.umi_coverage_index_path,
                    snps_path,
                    quality=20,
                    threshold=0.5,
                    n_threads=2
                )
            )
            self.assertTrue(mixins.files_equal(self.umi_snps_path, snps_path))
