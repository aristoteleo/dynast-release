import gzip
import os
import pickle
import subprocess as sp
import tempfile
import uuid
from unittest import mock, TestCase

import numpy as np
import pandas as pd

import dynast.utils as utils
from tests import mixins


class TestUtils(mixins.TestMixin, TestCase):

    def test_run_executable(self):
        p = utils.run_executable(['echo', 'TEST'], stream=False)
        self.assertEqual(p.stdout.read(), 'TEST\n')

    def test_run_exectuable_raises_exception(self):
        with self.assertRaises(sp.SubprocessError):
            utils.run_executable(['bash', 'nonexistent option'])

    def test_run_exectuable_with_returncode(self):
        utils.run_executable(['bash', 'nonexistent option'], returncode=127)

    def test_run_executable_no_wait(self):
        with mock.patch('dynast.utils.sp') as sp_mock:
            sp_mock.Popen().returncode = 0
            utils.run_executable(['echo', 'TEST'], wait=False)
            sp_mock.Popen().poll.assert_not_called()

    def test_run_executable_with_stream(self):
        with mock.patch('dynast.utils.logger.debug') as debug_mock:
            utils.run_executable(['echo', 'TEST'], stream=True)
            debug_mock.assert_has_calls([mock.call('TEST')])

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

    def test_open_as_text_textfile(self):
        path = os.path.join(tempfile.gettempdir(), f'{uuid.uuid4()}.txt')
        with utils.open_as_text(path, 'w') as f:
            f.write('TESTING')
        self.assertTrue(os.path.exists(path))
        with utils.open_as_text(path, 'r') as f:
            self.assertEqual(f.read(), 'TESTING')

    def test_open_as_text_gzip(self):
        path = os.path.join(tempfile.gettempdir(), f'{uuid.uuid4()}.gz')
        with utils.open_as_text(path, 'w') as f:
            f.write('TESTING')
        self.assertTrue(os.path.exists(path))
        with utils.open_as_text(path, 'r') as f:
            self.assertEqual(f.read(), 'TESTING')

    def test_decompress_gzip(self):
        filename = str(uuid.uuid4())
        gzip_path = os.path.join(tempfile.gettempdir(), f'{filename}.gz')
        out_path = os.path.join(tempfile.gettempdir(), filename)
        with gzip.open(gzip_path, 'wt') as f:
            f.write('TESTING\nTEST')
        self.assertEqual(out_path, utils.decompress_gzip(gzip_path, out_path))
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, 'r') as f:
            self.assertEqual('TESTING\nTEST', f.read())

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

    def test_flatten_dict_values(self):
        d = {'key1': 'value1', 'key2': {'key3': {'key4': 'value2', 'key5': 'value3'}, 'key6': 'value4'}}
        self.assertEqual(['value1', 'value2', 'value3', 'value4'], utils.flatten_dict_values(d))

    def test_mkstemp(self):
        path = utils.mkstemp()
        self.assertTrue(os.path.exists(path))

    def test_get_available_memory(self):
        with mock.patch('dynast.utils.psutil.virtual_memory') as vm:
            self.assertEqual(vm.return_value.available, utils.get_available_memory())

    def test_all_exists(self):
        file1 = os.path.join(self.temp_dir, 'file1')
        file2 = os.path.join(self.temp_dir, 'file2')
        file3 = os.path.join(self.temp_dir, 'file3')
        with open(file1, 'w') as f:
            f.write('')
        with open(file2, 'w') as f:
            f.write('')
        self.assertTrue(utils.all_exists([file1, file2]))
        self.assertFalse(utils.all_exists([file1, file2, file3]))

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
        with mock.patch('dynast.utils.tqdm'):
            utils.display_progress_with_counter(counter, total, async_result)

    def test_as_completed_with_progress(self):
        futures = [mock.MagicMock(), mock.MagicMock()]
        with mock.patch('dynast.utils.tqdm'), \
            mock.patch('dynast.utils.as_completed') as as_completed:
            as_completed.return_value = ['future1', 'future2']
            self.assertEqual(['future1', 'future2'], list(utils.as_completed(futures)))
            as_completed.assert_called_once_with(futures)

    def test_flatten_dictionary(self):
        d = {'a': 'b', 'c': 'd', 'e': {'f': 'g', 'h': {'i': 'j'}}}
        self.assertEqual([(('a',), 'b'), (('c',), 'd'), (('e', 'f'), 'g'), (('e', 'h', 'i'), 'j')],
                         list(utils.flatten_dictionary(d)))

    def test_flatten_list(self):
        lst = [1, [2, 3], [4, [5, 6]]]
        self.assertEqual([1, 2, 3, 4, 5, 6], list(utils.flatten_list(lst)))

    def test_merge_dictionaries(self):
        d1 = {'a': 'b', 'c': {'d': 'e'}, 'f': 'g'}
        d2 = {'a': 'h', 'c': {'i': 'j'}}
        self.assertEqual({
            'a': 'bh',
            'c': {
                'd': 'eX',
                'i': 'Xj'
            },
            'f': 'gX'
        }, utils.merge_dictionaries(d1, d2, default='X'))

    def test_write_pickle(self):
        path = os.path.join(self.temp_dir, 'pkl')
        utils.write_pickle('test', path)
        with gzip.open(path, 'rb') as f:
            self.assertEqual('test', pickle.load(f))

    def test_read_pickle(self):
        path = os.path.join(self.temp_dir, 'pkl')
        with gzip.open(path, 'wb') as f:
            pickle.dump('test', f)
        self.assertEqual('test', utils.read_pickle(path))

    def test_split_index(self):
        index = [(0, 5), (1, 5), (2, 5), (3, 5), (4, 5)]
        self.assertEqual([(0, 15), (3, 10)], utils.split_index(index, n=2))

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
