import os
from unittest import mock, TestCase

import pandas as pd

import dynast.preprocessing.bam as bam
import dynast.preprocessing.coverage as coverage
import dynast.preprocessing.snp as snp

from .. import mixins


class TestSnp(mixins.TestMixin, TestCase):

    def test_read_snps(self):
        snp.read_snps(self.control_snps_path)

    def test_read_snp_csv(self):
        snp.read_snp_csv(self.control_snps_path)

    def test_extract_conversions(self):
        alignments = bam.select_alignments(bam.read_alignments(self.control_alignments_path))
        with mock.patch('dynast.preprocessing.snp.utils.display_progress_with_counter'):
            conversions = snp.extract_conversions(
                self.control_conversions_path,
                self.control_conversions_index_path,
                alignments=alignments,
                quality=35,
                n_threads=2,
            )

        df = pd.read_csv(self.control_conversions_path)
        conversions_truth = {}
        df_indexed = df.set_index(['read_id', 'index'])
        df_selected = df_indexed[(df_indexed['quality'] > 35) & (df_indexed.index.isin(alignments))].reset_index()
        for (contig, genome_i, original, converted), count in dict(df_selected.groupby(['contig', 'genome_i',
                                                                                        'original',
                                                                                        'converted']).size()).items():
            conversion = f'{original}{converted}'
            conversions_truth.setdefault(conversion, {}).setdefault(contig, {}).setdefault(genome_i, 0)
            conversions_truth[conversion][contig][genome_i] += 1
        self.assertEqual(conversions_truth, conversions)

    def test_detect_snps(self):
        alignments = bam.select_alignments(bam.read_alignments(self.control_alignments_path))
        cov = coverage.read_coverage(self.control_coverage_path)
        snps_path = os.path.join(self.temp_dir, 'snps.csv')
        with mock.patch('dynast.preprocessing.snp.utils.display_progress_with_counter'):
            self.assertEqual(
                snps_path,
                snp.detect_snps(
                    self.control_conversions_path,
                    self.control_conversions_index_path,
                    cov,
                    snps_path,
                    alignments=alignments,
                    conversions={'TC', 'AG'},
                    quality=27,
                    threshold=0.5,
                    n_threads=2
                )
            )
            self.assertTrue(mixins.files_equal(self.control_snps_path, snps_path))
