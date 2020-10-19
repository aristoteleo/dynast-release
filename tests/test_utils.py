import gzip
import os
import subprocess as sp
import tempfile
import uuid
from unittest import mock, TestCase

import numpy as np

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

    def test_read_STAR_count_matrix(self):
        adata = utils.read_STAR_count_matrix(self.STAR_barcodes_path, self.STAR_features_path, self.STAR_matrix_path)
        self.assertEqual(2, np.count_nonzero(adata.X))
        self.assertEqual(3, adata.X[0, 0])
        self.assertEqual(5, adata.X[3, 2])
        self.assertEqual(['AA', 'AC', 'AG', 'AT', 'CA'], list(adata.obs.barcode))
        self.assertEqual(['G1', 'G2', 'G3', 'G4', 'G5'], list(adata.var.gene_id))
        self.assertEqual(['N1', 'N2', 'N3', 'N4', 'N5'], list(adata.var.gene_name))
