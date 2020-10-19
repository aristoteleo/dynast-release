import multiprocessing
import os
from unittest import mock, TestCase

import dynast.preprocessing.bam as bam
from tests import mixins


class TestBam(mixins.TestMixin, TestCase):

    def test_parse_read_contig(self):
        counter = multiprocessing.Value('I', 0)
        lock = multiprocessing.Lock()
        contig = '1'
        conversions_path, index_path, coverage_path = bam.parse_read_contig(
            self.bam_path, counter, lock, contig, temp_dir=self.temp_dir
        )

        self.assertTrue(mixins.files_equal(self.bam_conversions_paths[0], conversions_path))
        self.assertTrue(mixins.files_equal(self.bam_index_paths[0], index_path, gzipped=True))
        self.assertTrue(mixins.files_equal(self.bam_coverage_paths[0], coverage_path))

    def test_parse_all_reads(self):
        conversions_path = os.path.join(self.temp_dir, 'conversions.csv')
        index_path = os.path.join(self.temp_dir, 'conversions.idx')
        coverage_path = os.path.join(self.temp_dir, 'coverage_path')
        with mock.patch('dynast.preprocessing.bam.parse_read_contig'), \
            mock.patch('dynast.preprocessing.bam.tqdm'), \
            mock.patch('dynast.preprocessing.bam.partial') as partial, \
            mock.patch('dynast.preprocessing.bam.multiprocessing') as mp:
            async_result = mock.MagicMock()
            async_result.ready.return_value = True
            async_result.get.return_value = zip(
                self.bam_conversions_paths, self.bam_index_paths, self.bam_coverage_paths
            )
            mp.Pool.return_value.map_async.return_value = async_result

            result = bam.parse_all_reads(self.bam_path, conversions_path, index_path, coverage_path)
            self.assertEqual((conversions_path, index_path, coverage_path), result)

            partial.assert_called_once_with(
                bam.parse_read_contig,
                self.bam_path,
                mp.Manager.return_value.Value.return_value,
                mp.Manager.return_value.Lock.return_value,
                temp_dir=mock.ANY
            )
            mp.Pool.return_value.map_async.assert_called_once_with(partial.return_value, ['1', '2'])

            self.assertTrue(mixins.files_equal(self.bam_conversions_path, conversions_path))
            self.assertTrue(mixins.files_equal(self.bam_index_path, index_path, gzipped=True))
            self.assertTrue(mixins.files_equal(self.bam_coverage_path, coverage_path))
