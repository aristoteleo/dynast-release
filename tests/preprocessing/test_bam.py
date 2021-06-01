import os
from unittest import mock, TestCase

import pandas as pd

import dynast.preprocessing.bam as bam
import dynast.utils as utils

from .. import mixins


class TestBam(mixins.TestMixin, TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestBam, cls).setUpClass()
        cls.gene_infos = utils.read_pickle(cls.genes_path)
        cls.transcript_infos = utils.read_pickle(cls.transcripts_path)

    def test_read_alignments(self):
        bam.read_alignments(self.umi_alignments_path)

    def test_read_conversions(self):
        bam.read_conversions(self.umi_conversions_path)

    def test_select_alignments(self):
        df = bam.read_alignments(self.umi_alignments_path)
        result = bam.select_alignments(df)
        df_indexed = df.set_index(['read_id', 'index'])
        df_deduplicated = df_indexed[df_indexed.index.isin(result)].reset_index()
        self.assertFalse(df_deduplicated.duplicated(['barcode', 'umi', 'GX'], keep=False).all())

    def test_select_alignments_paired(self):
        df = bam.read_alignments(self.umi_alignments_path)
        result = bam.select_alignments(df)
        df_indexed = df.set_index(['read_id', 'index'])
        df_deduplicated = df_indexed[df_indexed.index.isin(result)].reset_index()
        self.assertFalse(df_deduplicated.duplicated('read_id', keep=False).all())

    def test_parse_all_reads(self):
        with mock.patch('dynast.preprocessing.bam.utils.display_progress_with_counter'):
            conversions_path = os.path.join(self.temp_dir, 'conversions.csv')
            index_path = os.path.join(self.temp_dir, 'conversions.idx')
            alignments_path = os.path.join(self.temp_dir, 'alignments.csv')

            self.assertEqual((conversions_path, alignments_path, index_path),
                             bam.parse_all_reads(
                                 self.umi_bam_path,
                                 conversions_path,
                                 alignments_path,
                                 index_path,
                                 self.gene_infos,
                                 self.transcript_infos,
                                 strand='forward',
                                 umi_tag='UB',
                                 barcode_tag='CB',
                                 gene_tag='GX',
                                 barcodes=None,
                                 n_threads=1,
                                 temp_dir=self.temp_dir,
                                 nasc=False,
                                 velocity=True
                             ))

            self.assertTrue(mixins.files_equal(self.umi_conversions_path, conversions_path))
            self.assertEqual(utils.read_pickle(self.umi_conversions_index_path), utils.read_pickle(index_path))
            self.assertTrue(mixins.files_equal(self.umi_alignments_path, alignments_path))

    def test_parse_all_reads_no_velocity(self):
        with mock.patch('dynast.preprocessing.bam.utils.display_progress_with_counter'):
            conversions_path = os.path.join(self.temp_dir, 'conversions.csv')
            index_path = os.path.join(self.temp_dir, 'conversions.idx')
            alignments_path = os.path.join(self.temp_dir, 'alignments.csv')

            self.assertEqual((conversions_path, alignments_path, index_path),
                             bam.parse_all_reads(
                                 self.umi_bam_path,
                                 conversions_path,
                                 alignments_path,
                                 index_path,
                                 self.gene_infos,
                                 self.transcript_infos,
                                 strand='forward',
                                 umi_tag='UB',
                                 barcode_tag='CB',
                                 gene_tag='GX',
                                 barcodes=None,
                                 n_threads=1,
                                 temp_dir=self.temp_dir,
                                 nasc=False,
                                 velocity=False
                             ))
            self.assertEqual(['unassigned'], list(pd.read_csv(alignments_path)['velocity'].unique()))

    def test_parse_all_reads_paired(self):
        with mock.patch('dynast.preprocessing.bam.utils.display_progress_with_counter'):
            conversions_path = os.path.join(self.temp_dir, 'conversions.csv')
            index_path = os.path.join(self.temp_dir, 'conversions.idx')
            alignments_path = os.path.join(self.temp_dir, 'alignments.csv')

            self.assertEqual((conversions_path, alignments_path, index_path),
                             bam.parse_all_reads(
                                 self.paired_bam_path,
                                 conversions_path,
                                 alignments_path,
                                 index_path,
                                 self.gene_infos,
                                 self.transcript_infos,
                                 strand='unstranded',
                                 umi_tag=None,
                                 barcode_tag='RG',
                                 gene_tag='GX',
                                 barcodes=None,
                                 n_threads=1,
                                 temp_dir=self.temp_dir,
                                 nasc=False,
                                 velocity=False
                             ))

            self.assertTrue(mixins.files_equal(self.paired_conversions_path, conversions_path))
            self.assertEqual(utils.read_pickle(self.paired_conversions_index_path), utils.read_pickle(index_path))
            self.assertTrue(mixins.files_equal(self.paired_alignments_path, alignments_path))

    def test_parse_all_reads_nasc(self):
        with mock.patch('dynast.preprocessing.bam.utils.display_progress_with_counter'):
            conversions_path = os.path.join(self.temp_dir, 'conversions.csv')
            index_path = os.path.join(self.temp_dir, 'conversions.idx')
            alignments_path = os.path.join(self.temp_dir, 'alignments.csv')

            self.assertEqual((conversions_path, alignments_path, index_path),
                             bam.parse_all_reads(
                                 self.nasc_bam_path,
                                 conversions_path,
                                 alignments_path,
                                 index_path,
                                 self.gene_infos,
                                 self.transcript_infos,
                                 strand='forward',
                                 umi_tag=None,
                                 barcode_tag='RG',
                                 gene_tag='GX',
                                 barcodes=None,
                                 n_threads=1,
                                 temp_dir=self.temp_dir,
                                 nasc=True,
                                 velocity=False
                             ))

            self.assertTrue(mixins.files_equal(self.nasc_conversions_path, conversions_path))
            self.assertEqual(utils.read_pickle(self.nasc_conversions_index_path), utils.read_pickle(index_path))
            self.assertTrue(mixins.files_equal(self.nasc_alignments_path, alignments_path))
