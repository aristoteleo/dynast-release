import os
from unittest import mock, TestCase

import pandas as pd

import dynast.preprocessing.bam as bam
import dynast.utils as utils

from .. import mixins


class TestBam(mixins.TestMixin, TestCase):

    def test_read_conversions(self):
        bam.read_conversions(self.umi_conversions_path)

    def test_read_conversions_no(self):
        bam.read_conversions(self.umi_no_conversions_path)

    def test_parse_all_reads(self):
        with mock.patch('dynast.preprocessing.bam.utils.display_progress_with_counter'):
            conversions_path = os.path.join(self.temp_dir, 'conversions.csv')
            index_path = os.path.join(self.temp_dir, 'conversions.idx')
            no_conversions_path = os.path.join(self.temp_dir, 'no_conversions.csv')
            no_index_path = os.path.join(self.temp_dir, 'no_conversions.idx')

            self.assertEqual((conversions_path, index_path, no_conversions_path, no_index_path),
                             bam.parse_all_reads(
                                 self.umi_bam_path,
                                 conversions_path,
                                 index_path,
                                 no_conversions_path,
                                 no_index_path,
                                 utils.read_pickle(self.umi_genes_path),
                                 utils.read_pickle(self.umi_transcripts_path),
                                 strand='forward',
                                 umi_tag='UB',
                                 barcode_tag='CB',
                                 gene_tag='GX',
                                 barcodes=None,
                                 n_threads=2,
                                 temp_dir=self.temp_dir,
                                 nasc=False,
                                 velocity=True
                             ))

            self.assertTrue(mixins.files_equal(self.umi_conversions_path, conversions_path))
            self.assertEqual(utils.read_pickle(self.umi_conversions_index_path), utils.read_pickle(index_path))
            self.assertTrue(mixins.files_equal(self.umi_no_conversions_path, no_conversions_path))
            self.assertEqual(utils.read_pickle(self.umi_no_conversions_index_path), utils.read_pickle(no_index_path))

    def test_parse_all_reads_no_velocity(self):
        with mock.patch('dynast.preprocessing.bam.utils.display_progress_with_counter'):
            conversions_path = os.path.join(self.temp_dir, 'conversions.csv')
            index_path = os.path.join(self.temp_dir, 'conversions.idx')
            no_conversions_path = os.path.join(self.temp_dir, 'no_conversions.csv')
            no_index_path = os.path.join(self.temp_dir, 'no_conversions.idx')

            self.assertEqual((conversions_path, index_path, no_conversions_path, no_index_path),
                             bam.parse_all_reads(
                                 self.umi_bam_path,
                                 conversions_path,
                                 index_path,
                                 no_conversions_path,
                                 no_index_path,
                                 gene_infos=utils.read_pickle(self.umi_genes_path),
                                 transcript_infos=utils.read_pickle(self.umi_transcripts_path),
                                 strand='forward',
                                 umi_tag='UB',
                                 barcode_tag='CB',
                                 gene_tag='GX',
                                 barcodes=None,
                                 n_threads=2,
                                 temp_dir=self.temp_dir,
                                 nasc=False,
                                 velocity=False
                             ))
            self.assertEqual(['unassigned'], list(pd.read_csv(conversions_path)['velocity'].unique()))
            self.assertEqual(['unassigned'], list(pd.read_csv(no_conversions_path)['velocity'].unique()))

    def test_parse_all_reads_paired(self):
        with mock.patch('dynast.preprocessing.bam.utils.display_progress_with_counter'):
            conversions_path = os.path.join(self.temp_dir, 'conversions.csv')
            index_path = os.path.join(self.temp_dir, 'conversions.idx')
            no_conversions_path = os.path.join(self.temp_dir, 'no_conversions.csv')
            no_index_path = os.path.join(self.temp_dir, 'no_conversions.idx')

            self.assertEqual((conversions_path, index_path, no_conversions_path, no_index_path),
                             bam.parse_all_reads(
                                 self.paired_bam_path,
                                 conversions_path,
                                 index_path,
                                 no_conversions_path,
                                 no_index_path,
                                 utils.read_pickle(self.paired_genes_path),
                                 utils.read_pickle(self.paired_transcripts_path),
                                 strand='unstranded',
                                 umi_tag=None,
                                 barcode_tag='RG',
                                 gene_tag='GX',
                                 barcodes=None,
                                 n_threads=2,
                                 temp_dir=self.temp_dir,
                                 nasc=False,
                                 velocity=False
                             ))

            self.assertTrue(mixins.files_equal(self.paired_conversions_path, conversions_path))
            self.assertEqual(utils.read_pickle(self.paired_conversions_index_path), utils.read_pickle(index_path))
            self.assertTrue(mixins.files_equal(self.paired_no_conversions_path, no_conversions_path))
            self.assertEqual(utils.read_pickle(self.paired_no_conversions_index_path), utils.read_pickle(no_index_path))

    def test_parse_all_reads_nasc(self):
        with mock.patch('dynast.preprocessing.bam.utils.display_progress_with_counter'):
            conversions_path = os.path.join(self.temp_dir, 'conversions.csv')
            index_path = os.path.join(self.temp_dir, 'conversions.idx')
            no_conversions_path = os.path.join(self.temp_dir, 'no_conversions.csv')
            no_index_path = os.path.join(self.temp_dir, 'no_conversions.idx')

            self.assertEqual((conversions_path, index_path, no_conversions_path, no_index_path),
                             bam.parse_all_reads(
                                 self.nasc_bam_path,
                                 conversions_path,
                                 index_path,
                                 no_conversions_path,
                                 no_index_path,
                                 utils.read_pickle(self.nasc_genes_path),
                                 utils.read_pickle(self.nasc_transcripts_path),
                                 strand='forward',
                                 umi_tag=None,
                                 barcode_tag='RG',
                                 gene_tag='GX',
                                 barcodes=None,
                                 n_threads=2,
                                 temp_dir=self.temp_dir,
                                 nasc=True,
                                 velocity=False
                             ))

            self.assertTrue(mixins.files_equal(self.nasc_conversions_path, conversions_path))
            self.assertEqual(utils.read_pickle(self.nasc_conversions_index_path), utils.read_pickle(index_path))
            self.assertTrue(mixins.files_equal(self.nasc_no_conversions_path, no_conversions_path))
            self.assertEqual(utils.read_pickle(self.nasc_no_conversions_index_path), utils.read_pickle(no_index_path))
