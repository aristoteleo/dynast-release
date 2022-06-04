from unittest import TestCase, mock

import numpy as np
import pandas as pd

import dynast.utils as utils
from tests import mixins


class TestUtils(mixins.TestMixin, TestCase):

    def test_get_STAR_binary_path(self):
        self.assertTrue(utils.get_STAR_binary_path().endswith('STAR'))

    def test_get_STAR_version(self):
        self.assertTrue('.' in utils.get_STAR_version())

    def test_combine_arguments(self):
        args = {'--arg1': 'value1', '--arg2': ['value2', 'value3'], '--arg3': ['value4'], '--arg4': 'value5'}
        additional = {'--arg1': 'value6', '--arg2': ['value7'], '--arg3': 'value8', '--arg5': 'value9'}
        self.assertEqual({
            '--arg1': 'value6',
            '--arg2': ['value2', 'value3', 'value7'],
            '--arg3': 'value8',
            '--arg4': 'value5',
            '--arg5': 'value9'
        }, utils.combine_arguments(args, additional))

    def test_arguments_to_list(self):
        args = {'--arg1': 'value1', '--arg2': ['value2', 'value3']}
        self.assertEqual(['--arg1', 'value1', '--arg2', 'value2', 'value3'], utils.arguments_to_list(args))

    def test_get_file_descriptor_limit_windows(self):
        import_mock = mock.MagicMock()
        with mock.patch('dynast.utils.config.PLATFORM', 'windows'), \
            mock.patch('builtins.__import__', side_effect=mixins.import_mock(import_mock)):
            self.assertEqual(import_mock._getmaxstdio.return_value, utils.get_file_descriptor_limit())

    def test_get_file_descriptor_limit_unix(self):
        import_mock = mock.MagicMock()
        with mock.patch('dynast.utils.config.PLATFORM', 'linux'), \
            mock.patch('builtins.__import__', side_effect=mixins.import_mock(import_mock)):
            self.assertEqual(import_mock.getrlimit.return_value[0], utils.get_file_descriptor_limit())

    def test_get_max_file_descriptor_limit_windows(self):
        with mock.patch('dynast.utils.config.PLATFORM', 'windows'):
            self.assertEqual(8192, utils.get_max_file_descriptor_limit())

    def test_get_max_file_descriptor_limit_unix(self):
        import_mock = mock.MagicMock()
        with mock.patch('dynast.utils.config.PLATFORM', 'linux'), \
            mock.patch('builtins.__import__', side_effect=mixins.import_mock(import_mock)):
            self.assertEqual(import_mock.getrlimit.return_value[1], utils.get_max_file_descriptor_limit())

    def test_increase_file_descriptor_limit_windows(self):
        import_mock = mock.MagicMock()
        with mock.patch('dynast.utils.config.PLATFORM', 'windows'), \
            mock.patch('builtins.__import__', side_effect=mixins.import_mock(import_mock)):

            import_mock._getmaxstdio.return_value = 'old'
            with utils.increase_file_descriptor_limit('new'):
                import_mock._setmaxstdio.assert_called_once_with('new')
            import_mock._setmaxstdio.assert_has_calls([mock.call('old')])

    def test_increase_file_descriptor_limit_unix(self):
        import_mock = mock.MagicMock()
        with mock.patch('dynast.utils.config.PLATFORM', 'linux'), \
            mock.patch('builtins.__import__', side_effect=mixins.import_mock(import_mock)):

            import_mock.getrlimit.return_value = ('old', 'max')
            with utils.increase_file_descriptor_limit('new'):
                import_mock.setrlimit.assert_called_once_with(mock.ANY, ('new', 'max'))
            import_mock.setrlimit.assert_has_calls([mock.call(mock.ANY, ('old', 'max'))])

    def test_get_available_memory(self):
        with mock.patch('dynast.utils.psutil.virtual_memory') as vm:
            self.assertEqual(vm.return_value.available, utils.get_available_memory())

    def test_make_pool_with_counter(self):
        with mock.patch('dynast.utils.multiprocessing.Manager') as Manager, \
            mock.patch('dynast.utils.multiprocessing.Pool') as Pool:
            self.assertEqual(
                (Pool.return_value, Manager.return_value.Value.return_value, Manager.return_value.Lock.return_value),
                utils.make_pool_with_counter(8)
            )
            Manager.return_value.Value.assert_called_once_with('I', 0)
            Pool.assert_called_once_with(8)

    def test_display_progress_with_counter(self):
        counter = mock.MagicMock()
        total = mock.MagicMock()
        async_result = mock.MagicMock()
        async_result.ready.return_value = True
        with mock.patch('dynast.utils.ngs.progress.progress'):
            utils.display_progress_with_counter(counter, total, async_result)

    def test_as_completed_with_progress(self):
        futures = [mock.MagicMock(), mock.MagicMock()]
        with mock.patch('dynast.utils.ngs.progress.progress'), \
            mock.patch('dynast.utils.as_completed') as as_completed:
            as_completed.return_value = ['future1', 'future2']
            self.assertEqual(['future1', 'future2'], list(utils.as_completed(futures)))
            as_completed.assert_called_once_with(futures)

    def test_split_index(self):
        index = [(0, 5), (1, 5), (2, 5), (3, 5), (4, 5)]
        self.assertEqual([[(0, 5), (1, 5), (2, 5)], [(3, 5), (4, 5)]], utils.split_index(index, n=2))

    def test_downsample_counts_with_count(self):
        df = pd.DataFrame([
            ['a', 'b', 'c'],
            ['d', 'e', 'f'],
            ['g', 'h', 'i'],
        ], columns=['col1', 'col2', 'col3'])
        df_downsampled = utils.downsample_counts(df, count=1, seed=0)
        pd.testing.assert_frame_equal(df.iloc[[2]], df_downsampled)

    def test_downsample_counts_with_proportion(self):
        df = pd.DataFrame([
            ['a', 'b', 'c'],
            ['d', 'e', 'f'],
        ], columns=['col1', 'col2', 'col3'])
        df_downsampled = utils.downsample_counts(df, proportion=0.5, seed=0)
        pd.testing.assert_frame_equal(df.iloc[[1]], df_downsampled)

    def test_downsample_counts_group_by(self):
        df = pd.DataFrame([
            ['a', 'b', 'c'],
            ['a', 'e', 'f'],
            ['b', 'e', 'f'],
        ], columns=['cell', 'col2', 'col3'])
        df_downsampled = utils.downsample_counts(df, count=1, seed=0, group_by=['cell'])
        pd.testing.assert_frame_equal(df.iloc[[1, 2]], df_downsampled)

    def test_counts_to_matrix(self):
        rows = [
            ['barcode1', 'GX1'],
            ['barcode1', 'GX1'],
            ['barcode2', 'GX2'],
        ]
        df = pd.DataFrame(rows, columns=['barcode', 'GX'])
        barcodes = ['barcode1', 'barcode2']
        features = ['GX1', 'GX2']
        matrix = utils.counts_to_matrix(df, barcodes, features)
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
        with mock.patch('dynast.utils.counts_to_matrix') as counts_to_matrix:
            self.assertEqual((counts_to_matrix.return_value, counts_to_matrix.return_value),
                             utils.split_counts(df, barcodes, features))
            self.assertEqual(2, counts_to_matrix.call_count)

    def test_split_matrix(self):
        matrix = np.array([[1, 2], [4, 5]])
        barcodes = ['bc1', 'bc2']
        features = ['gx1', 'gx2']
        pis = {('bc1', 'gx1'): 0.5, ('bc2', 'gx1'): 0.25}
        pi_mask, unlabeled_matrix, labeled_matrix = utils.split_matrix(matrix, pis, barcodes, features)
        self.assertTrue(np.array_equal([[True, False], [True, False]], pi_mask.A))
        self.assertTrue(np.array_equal([[0.5, 0], [1, 0]], labeled_matrix.A))
        self.assertTrue(np.array_equal([[0.5, 0], [3, 0]], unlabeled_matrix.A))
