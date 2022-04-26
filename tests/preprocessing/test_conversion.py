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
        df['new_column'] = 'test'
        gene_infos = utils.read_pickle(self.umi_genes_path)
        df_complemented = conversion.complement_counts(df, gene_infos)
        self.assertEqual(df.shape[0], df_complemented.shape[0])
        self.assertIn('new_column', df_complemented.columns)
        # Select one gene to investigate
        gene = list(df['GX'][df['GX'].map(lambda gx: gene_infos[gx]['strand']) == '-'])[0]
        df_gene = df[df['GX'] == gene].sort_values(['barcode', 'GX'])
        df_complemented_gene = df_complemented[df_complemented['GX'] == gene].sort_values(['barcode', 'GX'])
        self.assertEqual(df_gene.iloc[0]['A'], df_complemented_gene.iloc[0]['T'])
        self.assertEqual(df_gene.iloc[0]['G'], df_complemented_gene.iloc[0]['C'])
        self.assertEqual(df_gene.iloc[0]['TC'], df_complemented_gene.iloc[0]['AG'])
        self.assertEqual(df_gene.iloc[0]['CG'], df_complemented_gene.iloc[0]['GC'])

    def test_complement_counts_inverse(self):
        df = pd.read_csv(self.umi_counts_path)
        gene_infos = utils.read_pickle(self.umi_genes_path)
        df_inverse = conversion.complement_counts(conversion.complement_counts(df, gene_infos), gene_infos)
        pd.testing.assert_frame_equal(df, df_inverse[df.columns])

    def test_drop_multimappers(self):
        rows = [
            ['read1', 'BC', 'GX', 'unassigned', True, 0] + [0] * len(conversion.COLUMNS),
            ['read2', 'BC', 'GX', 'unassigned', False, 0] + [0] * len(conversion.COLUMNS),
        ]
        df = pd.DataFrame(
            rows, columns=['read_id', 'barcode', 'GX', 'velocity', 'transcriptome', 'score'] + conversion.COLUMNS
        )
        df['velocity'] = df['velocity'].astype('category')
        df_dropped = conversion.drop_multimappers(df)
        self.assertEqual(2, df_dropped.shape[0])
        pd.testing.assert_frame_equal(df.set_index('read_id'), df_dropped.set_index('read_id'), check_like=True)

    def test_drop_multimappers_transcriptome(self):
        rows = [
            ['read1', 'BC', 'GX', 'unassigned', True, 0] + [0] * len(conversion.COLUMNS),
            ['read1', 'BC', 'GX', 'unassigned', False, 0] + [0] * len(conversion.COLUMNS),
        ]
        df = pd.DataFrame(
            rows, columns=['read_id', 'barcode', 'GX', 'velocity', 'transcriptome', 'score'] + conversion.COLUMNS
        )
        df['velocity'] = df['velocity'].astype('category')
        df_dropped = conversion.drop_multimappers(df)
        self.assertEqual(1, df_dropped.shape[0])
        pd.testing.assert_series_equal(df.iloc[0], df_dropped.iloc[0])

    def test_drop_multimappers_multiple_genes(self):
        rows = [
            ['read1', 'BC', 'GX1', 'unassigned', False, 0] + [0] * len(conversion.COLUMNS),
            ['read1', 'BC', 'GX2', 'unassigned', False, 0] + [0] * len(conversion.COLUMNS),
        ]
        df = pd.DataFrame(
            rows, columns=['read_id', 'barcode', 'GX', 'velocity', 'transcriptome', 'score'] + conversion.COLUMNS
        )
        df['velocity'] = df['velocity'].astype('category')
        df_dropped = conversion.drop_multimappers(df)
        self.assertEqual(0, df_dropped.shape[0])

    def test_drop_multimappers_multiple_velocity(self):
        rows = [
            ['read1', 'BC', 'GX1', 'spliced', False, 0] + [0] * len(conversion.COLUMNS),
            ['read1', 'BC', 'GX1', 'unspliced', False, 0] + [0] * len(conversion.COLUMNS),
        ]
        df = pd.DataFrame(
            rows, columns=['read_id', 'barcode', 'GX', 'velocity', 'transcriptome', 'score'] + conversion.COLUMNS
        )
        df['velocity'] = df['velocity'].astype('category')
        df_dropped = conversion.drop_multimappers(df)
        self.assertEqual(1, df_dropped.shape[0])
        self.assertEqual('ambiguous', df_dropped.iloc[0]['velocity'])

    def test_deduplicate_counts_transcriptome(self):
        rows = [
            ['barcode', 'umi', 'GX', 0] + [0] * len(conversion.COLUMNS) + [True],
            ['barcode', 'umi', 'GX', 0] + [0] * len(conversion.COLUMNS) + [False],
        ]
        df = pd.DataFrame(rows, columns=['barcode', 'umi', 'GX', 'score'] + conversion.COLUMNS + ['transcriptome'])
        df_deduplicated = conversion.deduplicate_counts(df)
        self.assertEqual(1, df_deduplicated.shape[0])
        self.assertTrue(df_deduplicated.iloc[0]['transcriptome'])

    def test_deduplicate_counts_score(self):
        rows = [
            ['barcode', 'umi', 'GX', 0] + [0] * len(conversion.COLUMNS) + [False],
            ['barcode', 'umi', 'GX', 1] + [0] * len(conversion.COLUMNS) + [False],
        ]
        df = pd.DataFrame(rows, columns=['barcode', 'umi', 'GX', 'score'] + conversion.COLUMNS + ['transcriptome'])
        df_deduplicated = conversion.deduplicate_counts(df)
        self.assertEqual(1, df_deduplicated.shape[0])
        pd.testing.assert_series_equal(df.iloc[1], df_deduplicated.iloc[0], check_names=False)

    def test_deduplicate_counts_conversion(self):
        rows = [
            ['barcode', 'umi', 'GX', 0] + [1] * len(conversion.CONVERSION_COLUMNS) +
            [0] * len(conversion.BASE_COLUMNS) + [True],
            ['barcode', 'umi', 'GX', 0] + [0] * len(conversion.CONVERSION_COLUMNS) +
            [0] * len(conversion.BASE_COLUMNS) + [True],
        ]
        df = pd.DataFrame(rows, columns=['barcode', 'umi', 'GX', 'score'] + conversion.COLUMNS + ['transcriptome'])
        df_deduplicated = conversion.deduplicate_counts(df)
        self.assertEqual(1, df_deduplicated.shape[0])
        self.assertEqual(1, df_deduplicated.iloc[0]['AC'])

    def test_split_counts_by_velocity(self):
        df = pd.read_csv(self.umi_counts_path)
        dfs = conversion.split_counts_by_velocity(df)
        self.assertEqual({'spliced', 'unspliced', 'ambiguous'}, set(dfs.keys()))
        self.assertEqual((df['velocity'] == 'spliced').sum(), dfs['spliced'].shape[0])
        self.assertEqual((df['velocity'] == 'unspliced').sum(), dfs['unspliced'].shape[0])
        self.assertEqual((df['velocity'] == 'ambiguous').sum(), dfs['ambiguous'].shape[0])

    def test_conversions(self):
        counts_path = os.path.join(self.temp_dir, 'counts.csv')
        with mock.patch('dynast.preprocessing.conversion.utils.display_progress_with_counter'):
            self.assertEqual(
                counts_path,
                conversion.count_conversions(
                    self.umi_conversions_path,
                    self.umi_alignments_path,
                    self.umi_conversions_index_path,
                    counts_path,
                    utils.read_pickle(self.umi_genes_path),
                    snps=None,
                    quality=27,
                    conversions=['TC'],
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
                    self.paired_alignments_path,
                    self.paired_conversions_index_path,
                    counts_path,
                    utils.read_pickle(self.paired_genes_path),
                    snps=None,
                    quality=27,
                    n_threads=2,
                    temp_dir=self.temp_dir
                )
            )
            self.assertTrue(mixins.files_equal(self.paired_counts_path, counts_path))
