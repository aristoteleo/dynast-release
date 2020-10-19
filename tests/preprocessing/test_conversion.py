import gzip
import multiprocessing
import os
import pickle
from unittest import mock, TestCase

import dynast.preprocessing.conversion as conversion
from tests import mixins


class TestConversion(mixins.TestMixin, TestCase):

    def test_split_index(self):
        index = [(1, 5), (2, 5), (3, 5), (4, 5), (5, 5)]
        self.assertEqual([(1, 15), (4, 10)], conversion.split_index(index, n=2))

    def test_count_conversions_part(self):
        counter = multiprocessing.Value('I', 0)
        lock = multiprocessing.Lock()
        with gzip.open(self.index_path, 'rb') as f:
            idx = pickle.load(f)
        pos = idx[0][0]
        n_lines = idx[0][1] + idx[1][1]

        counts_path, barcodes, genes = conversion.count_conversions_part(
            self.conversions_path, counter, lock, pos, n_lines, quality=0, temp_dir=self.temp_dir
        )
        self.assertEqual({'CCCCCTCTGTCA': 'CAATCTACTTCT', 'AGCGCGAATTAT': 'ACTAGTACCCAG'}, barcodes)
        self.assertEqual({'ENSG00000227232': ('WASH7P', '-'), 'ENSG00000186092': ('OR4F5', '+')}, genes)
        self.assertTrue(mixins.files_equal(self.counts_part_path, counts_path))

    def test_count_coversions(self):
        barcodes_path = os.path.join(self.temp_dir, 'barcodes.csv')
        genes_path = os.path.join(self.temp_dir, 'genes.csv')
        counts_path = os.path.join(self.temp_dir, 'counts.csv')
        with gzip.open(self.index_path, 'rb') as f:
            idx = pickle.load(f)
        parts = [(idx[0][0], sum(i[1] for i in idx[:20])), (idx[20][0], sum(i[1] for i in idx[20:]))]
        with mock.patch('dynast.preprocessing.conversion.split_index') as split_index, \
            mock.patch('dynast.preprocessing.conversion.tqdm'), \
            mock.patch('dynast.preprocessing.conversion.partial') as partial, \
            mock.patch('dynast.preprocessing.conversion.multiprocessing') as mp:
            split_index.return_value = parts
            async_result = mock.MagicMock()
            async_result.ready.return_value = True
            async_result.get.return_value = zip(
                self.counts_part_paths, [{
                    'CR1': 'CB1'
                }, {
                    'CR2': 'CB2'
                }], [{
                    'GX1': ('GN1', '+')
                }, {
                    'GX2': ('GN2', '-')
                }]
            )
            mp.Pool.return_value.starmap_async.return_value = async_result

            result = conversion.count_conversions(
                self.conversions_path,
                self.index_path,
                barcodes_path,
                genes_path,
                counts_path,
                quality=0,
                temp_dir=self.temp_dir
            )
            self.assertEqual((barcodes_path, genes_path, counts_path), result)

            partial.assert_called_once_with(
                conversion.count_conversions_part,
                self.conversions_path,
                mp.Manager.return_value.Value.return_value,
                mp.Manager.return_value.Lock.return_value,
                use_corrected=False,
                quality=0,
                temp_dir=mock.ANY,
            )
            mp.Pool.return_value.starmap_async.assert_called_once_with(partial.return_value, parts)
            self.assertTrue(mixins.files_equal(self.counts_path, counts_path))
            with open(genes_path, 'r') as f:
                self.assertEqual('GX,GN,strand\nGX1,GN1,+\nGX2,GN2,-\n', f.read())
            with open(barcodes_path, 'r') as f:
                self.assertEqual('CR,CB\nCR1,CB1\nCR2,CB2\n', f.read())
