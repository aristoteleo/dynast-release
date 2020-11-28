import gzip
import logging
import multiprocessing
import os
import pickle
import shutil
import subprocess as sp
import tempfile
import time
from contextlib import contextmanager
from operator import add

import anndata
import numpy as np
import pandas as pd
import psutil
import scipy.io
from tqdm import tqdm

from . import config

logger = logging.getLogger(__name__)


class TqdmLoggingHandler(logging.Handler):
    """Custom logging handler so that logging does not affect progress bars.
    """

    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception:
            self.handleError(record)


class suppress_stdout_stderr:
    """A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).
    https://github.com/facebook/prophet/issues/223
    """

    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


def run_executable(
    command,
    stdin=None,
    stdout=sp.PIPE,
    stderr=sp.PIPE,
    wait=True,
    stream=True,
    quiet=False,
    returncode=0,
    alias=True,
):
    """Execute a single shell command.

    :param command: a list representing a single shell command
    :type command: list
    :param stdin: object to pass into the `stdin` argument for `subprocess.Popen`,
                  defaults to `None`
    :type stdin: stream, optional
    :param stdout: object to pass into the `stdout` argument for `subprocess.Popen`,
                  defaults to `subprocess.PIPE`
    :type stdout: stream, optional
    :param stderr: object to pass into the `stderr` argument for `subprocess.Popen`,
                  defaults to `subprocess.PIPE`
    :type stderr: stream, optional
    :param wait: whether to wait until the command has finished, defaults to `True`
    :type wait: bool, optional
    :param stream: whether to stream the output to the command line, defaults to `True`
    :type stream: bool, optional
    :param quiet: whether to not display anything to the command line and not check the return code,
                  defaults to `False`
    :type quiet: bool, optional
    :param returncode: the return code expected if the command runs as intended,
                       defaults to `0`
    :type returncode: int, optional
    :param alias: whether to use the basename of the first element of `command`,
                  defaults to `True`
    :type alias: bool, optional

    :raise: `subprocess.CalledProcessError` if not `quiet` and the process
            exited with an exit code != `exitcode`

    :return: the spawned process
    :rtype: subprocess.Process
    """
    command = [str(c) for c in command]
    if not quiet:
        c = command.copy()
        if alias:
            c[0] = os.path.basename(c[0])
        logger.debug(' '.join(c))
    p = sp.Popen(
        command,
        stdin=stdin,
        stdout=stdout,
        stderr=stderr,
        universal_newlines=wait,
        bufsize=1 if wait else -1,
    )

    # Wait if desired.
    if wait:
        out = []
        while p.poll() is None:
            if stream and not quiet:
                for _line in p.stdout:
                    line = _line.strip()
                    out.append(line)
                    logger.debug(line)
                for _line in p.stderr:
                    line = _line.strip()
                    out.append(line)
                    logger.debug(line)
            else:
                time.sleep(1)

        if not quiet and p.returncode != returncode:
            logger.error('\n'.join(out))
            raise sp.CalledProcessError(p.returncode, ' '.join(command))

    return p


def open_as_text(path, mode):
    """Open a textfile or gzip file in text mode.

    :param path: path to textfile or gzip
    :type path: str
    :param mode: mode to open the file, either `w` for write or `r` for read
    :type mode: str

    :return: file object
    :rtype: file object
    """
    return gzip.open(path, mode + 't') if path.endswith('.gz') else open(path, mode)


def decompress_gzip(gzip_path, out_path):
    """Decompress a gzip file to provided file path.

    :param gzip_path: path to gzip file
    :type gzip_path: str
    :param out_path: path to decompressed file
    :type out_path: str

    :return: path to decompressed file
    :rtype: str
    """
    with gzip.open(gzip_path, 'rb') as f, open(out_path, 'wb') as out:
        shutil.copyfileobj(f, out)
    return out_path


def get_file_descriptor_limit():
    """Get the current value for the maximum number of open file descriptors
    in a platform-dependent way.

    :return: the current value of the maximum number of open file descriptors.
    :rtype: int
    """
    if config.PLATFORM == 'windows':
        import win32file
        return win32file._getmaxstdio()
    else:
        import resource
        return resource.getrlimit(resource.RLIMIT_NOFILE)[0]


def get_max_file_descriptor_limit():
    """Get the maximum allowed value for the maximum number of open file
    descriptors. Note that for Windows, there is not an easy way to get this,
    as it requires reading from the registry. So, we just return the maximum for
    a vanilla Windows installation, which is 8192.
    https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/setmaxstdio?view=vs-2019

    :return: maximum allowed value for the maximum number of open file descriptors
    :rtype: int
    """
    if config.PLATFORM == 'windows':
        return 8192
    else:
        import resource
        return resource.getrlimit(resource.RLIMIT_NOFILE)[1]


@contextmanager
def increase_file_descriptor_limit(limit):
    """Context manager that can be used to temporarily increase the maximum
    number of open file descriptors for the current process. The original
    value is restored when execution exits this function.

    This is required when running STAR with many threads.

    :param limit: maximum number of open file descriptors will be increased to
                  this value for the duration of the context
    :type limit: int
    """
    old = None
    if config.PLATFORM == 'windows':
        import win32file
        try:
            old = win32file._getmaxstdio()
            win32file._setmaxstdio(limit)
            yield
        finally:
            if old is not None:
                win32file._setmaxstdio(old)
    else:
        import resource
        try:
            old = resource.getrlimit(resource.RLIMIT_NOFILE)
            resource.setrlimit(resource.RLIMIT_NOFILE, (limit, old[1]))
            yield
        finally:
            if old is not None:
                resource.setrlimit(resource.RLIMIT_NOFILE, old)


def flatten_dict_values(d):
    """Extract all values from a nested dictionary.

    :param d: nested dictionary from which to extract values from
    :type d: dict

    :return: all values from the dictionary as a list
    :rtype: list
    """
    if isinstance(d, dict):
        flattened = []
        for k, v in d.items():
            if isinstance(v, dict):
                flattened.extend(flatten_dict_values(v))
            else:
                flattened.append(v)
        return flattened
    else:
        return [d]


def mkstemp(dir=None, delete=False):
    """Wrapper for `tempfile.mkstemp` that automatically closes the OS-level
    file descriptor. This function behaves like `tempfile.mkdtemp` but for
    files.

    :param dir: directory to create the temporary file. This value is passed as
                the `dir` kwarg of `tempfile.mkstemp`, defaults to `None`
    :type dir: str, optional
    :param delete: whether to delete the temporary file before returning,
                   defaults to `False`
    :type delete: bool, optional

    :return: path to the temporary file
    :rtype: str
    """
    fd, path = tempfile.mkstemp(dir=dir)
    os.close(fd)
    if delete:
        os.remove(path)
    return path


def get_available_memory():
    """Get total amount of available memory (total memory - used memory) in bytes.

    :return: available memory in bytes
    :rtype: int
    """
    return psutil.virtual_memory().available


def all_exists(paths):
    """Check if all provided paths exist.

    :param paths: list of paths
    :type paths: list

    :return: `True` iff all paths exist
    :rtype: bool
    """
    return all(os.path.exists(path) for path in paths)


def read_STAR_count_matrix(barcodes_path, features_path, matrix_path):
    """Given the barcodes, features, and matrix paths of a STAR count matrix,
    read it as an AnnData object.

    :param barcodes_path: path to barcodes.tsv
    :type barcodes_path: str
    :param features_path: path to features.tsv
    :type features_path: str
    :param matrix_path: path to matrix.mtx
    :type matrix_path: str

    :return: cells x genes AnnData matrix
    :rtype: AnnData
    """
    df_barcodes = pd.read_csv(barcodes_path, names=['barcode'])
    df_features = pd.read_csv(features_path, names=['gene_id', 'gene_name'], sep='\t', usecols=[0, 1])
    matrix = scipy.io.mmread(matrix_path).astype(np.uint32).T.toarray()

    adata = anndata.AnnData(matrix, obs=df_barcodes, var=df_features)
    adata.obs.index = adata.obs.index.astype(str)
    adata.var.index = adata.var.index.astype(str)
    return adata


def make_pool_with_counter(n_threads):
    manager = multiprocessing.Manager()
    counter = manager.Value('I', 0)
    lock = manager.Lock()
    pool = multiprocessing.Pool(n_threads)
    return pool, counter, lock


def display_progress_with_counter(async_result, counter, total):
    with tqdm(total=total, ascii=True, unit_scale=True) as pbar:
        previous_progress = 0
        while not async_result.ready():
            time.sleep(0.05)
            progress = counter.value
            pbar.update(progress - previous_progress)
            previous_progress = progress


def merge_dictionaries(d1, d2, f=add, default=0):
    merged = d1.copy()
    for key, value2 in d2.items():
        if key in d1 and isinstance(d1[key], dict) != isinstance(value2, dict):
            raise Exception(f'Inconsistent key {key}')

        value1 = d1.get(key, {} if isinstance(value2, dict) else (default() if callable(default) else default))
        if isinstance(value1, dict) and isinstance(value2, dict):
            merged[key] = merge_dictionaries(value1, value2, f=f, default=default)
        else:
            merged[key] = f(value1, value2)
    return merged


def flatten_dictionary(d, keys=None):
    keys = keys or []
    for k, v in d.items():
        new_keys = keys + [k]
        if isinstance(v, dict):
            yield from flatten_dictionary(v, new_keys)
        else:
            yield new_keys, v


def write_pickle(obj, path, *args, **kwargs):
    with gzip.open(path, 'wb') as f:
        pickle.dump(obj, f, *args, **kwargs)
    return path


def read_pickle(path):
    with gzip.open(path, 'rb') as f:
        return pickle.load(f)
