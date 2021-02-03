import os
from unittest import mock, TestCase

import pandas as pd

import dynast.preprocessing.conversion as conversion
import dynast.utils as utils

from .. import mixins


class TestConversion(mixins.TestMixin, TestCase):

    def test_read_counts(self):
        conversion.read_counts(self.umi_counts_path)

    def test_complement_counts(self):
        df = pd.read_csv(self.umi_counts_path)
        gene_infos = utils.read_pickle(self.umi_genes_path)
        df_complemented = conversion.complement_counts(df, gene_infos)
        self.assertEqual(df.shape[0], df_complemented.shape[0])
        # Select one gene to investigate
        gene = list(df['GX'][df['GX'].map(lambda gx: gene_infos[gx]['strand']) == '-'])[0]
        df_gene = df[df['GX'] == gene].sort_values(['barcode', 'GX'])
        df_complemented_gene = df_complemented[df_complemented['GX'] == gene].sort_values(['barcode', 'GX'])
        self.assertEqual(df_gene.iloc[0]['A'], df_complemented_gene.iloc[0]['T'])
        self.assertEqual(df_gene.iloc[0]['G'], df_complemented_gene.iloc[0]['C'])
        self.assertEqual(df_gene.iloc[0]['TC'], df_complemented_gene.iloc[0]['AG'])
        self.assertEqual(df_gene.iloc[0]['CG'], df_complemented_gene.iloc[0]['GC'])

    def test_deduplicate_counts_transcriptome(self):
        rows = [
            ['barcode', 'umi', 'GX'] + [0] * len(conversion.COLUMNS) + [True],
            ['barcode', 'umi', 'GX'] + [0] * len(conversion.COLUMNS) + [False],
        ]
        df = pd.DataFrame(rows, columns=['barcode', 'umi', 'GX'] + conversion.COLUMNS + ['transcriptome'])
        df_deduplicated = conversion.deduplicate_counts(df)
        self.assertEqual(1, df_deduplicated.shape[0])
        self.assertTrue(df_deduplicated.iloc[0]['transcriptome'])

    def test_deduplicate_counts_alignment(self):
        rows = [
            ['barcode', 'umi', 'GX'] + [0] * len(conversion.CONVERSION_COLUMNS) + [1] * len(conversion.BASE_COLUMNS) +
            [True],
            ['barcode', 'umi', 'GX'] + [0] * len(conversion.CONVERSION_COLUMNS) + [0] * len(conversion.BASE_COLUMNS) +
            [True],
        ]
        df = pd.DataFrame(rows, columns=['barcode', 'umi', 'GX'] + conversion.COLUMNS + ['transcriptome'])
        df_deduplicated = conversion.deduplicate_counts(df)
        self.assertEqual(1, df_deduplicated.shape[0])
        self.assertEqual(1, df_deduplicated.iloc[0]['A'])

    def test_deduplicate_counts_conversion(self):
        rows = [
            ['barcode', 'umi', 'GX'] + [1] * len(conversion.CONVERSION_COLUMNS) + [0] * len(conversion.BASE_COLUMNS) +
            [True],
            ['barcode', 'umi', 'GX'] + [0] * len(conversion.CONVERSION_COLUMNS) + [0] * len(conversion.BASE_COLUMNS) +
            [True],
        ]
        df = pd.DataFrame(rows, columns=['barcode', 'umi', 'GX'] + conversion.COLUMNS + ['transcriptome'])
        df_deduplicated = conversion.deduplicate_counts(df)
        self.assertEqual(1, df_deduplicated.shape[0])
        self.assertEqual(0, df_deduplicated.iloc[0]['AC'])

    def test_split_counts_by_velocity(self):
        df = pd.read_csv(self.umi_counts_path)
        dfs = conversion.split_counts_by_velocity(df)
        self.assertEqual({'spliced', 'unspliced', 'ambiguous'}, set(dfs.keys()))
        self.assertEqual((df['velocity'] == 'spliced').sum(), dfs['spliced'].shape[0])
        self.assertEqual((df['velocity'] == 'unspliced').sum(), dfs['unspliced'].shape[0])
        self.assertEqual((df['velocity'] == 'ambiguous').sum(), dfs['ambiguous'].shape[0])

    def test_counts_to_matrix(self):
        rows = [
            ['barcode1', 'GX1'],
            ['barcode1', 'GX1'],
            ['barcode2', 'GX2'],
        ]
        df = pd.DataFrame(rows, columns=['barcode', 'GX'])
        barcodes = ['barcode1', 'barcode2']
        features = ['GX1', 'GX2']
        matrix = conversion.counts_to_matrix(df, barcodes, features)
        self.assertEqual([[2, 0], [0, 1]], matrix.toarray().tolist())

    def test_split_counts(self):
        rows = [
            ['barcode1', 'GX1', 0],
            ['barcode1', 'GX1', 1],
            ['barcode2', 'GX2', 0],
        ]
        df = pd.DataFrame(rows, columns=['barcode', 'GX', 'TC'])
        barcodes = ['barcode1', 'barcode2']
        features = ['GX1', 'GX2']
        with mock.patch('dynast.preprocessing.conversion.counts_to_matrix') as counts_to_matrix:
            self.assertEqual((counts_to_matrix.return_value, counts_to_matrix.return_value),
                             conversion.split_counts(df, barcodes, features))
            self.assertEqual(2, counts_to_matrix.call_count)

    def test_conversions(self):
        counts_path = os.path.join(self.temp_dir, 'counts.csv')
        with mock.patch('dynast.preprocessing.conversion.utils.display_progress_with_counter'):
            self.assertEqual(
                counts_path,
                conversion.count_conversions(
                    self.umi_conversions_path,
                    self.umi_conversions_index_path,
                    self.umi_no_conversions_path,
                    self.umi_no_conversions_index_path,
                    counts_path,
                    snps=None,
                    quality=20,
                    n_threads=2,
                    temp_dir=self.temp_dir
                )
            )
            self.assertTrue(mixins.files_equal(self.umi_counts_path, counts_path))

    def test_conversions_paired(self):
        counts_path = os.path.join(self.temp_dir, 'counts.csv')
        with mock.patch('dynast.preprocessing.conversion.utils.display_progress_with_counter'):
            self.assertEqual(
                counts_path,
                conversion.count_conversions(
                    self.paired_conversions_path,
                    self.paired_conversions_index_path,
                    self.paired_no_conversions_path,
                    self.paired_no_conversions_index_path,
                    counts_path,
                    snps=None,
                    quality=27,
                    n_threads=2,
                    temp_dir=self.temp_dir
                )
            )
            self.assertTrue(mixins.files_equal(self.paired_counts_path, counts_path))