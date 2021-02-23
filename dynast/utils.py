import gzip
import multiprocessing
import os
import pickle
import shutil
import subprocess as sp
import tempfile
import time
from concurrent.futures import as_completed
from contextlib import contextmanager
from operator import add

import psutil
from tqdm import tqdm

from . import config
from .logging import logger


class UnsupportedOSException(Exception):
    pass


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


def get_STAR_binary_path():
    """Get the path to the platform-dependent STAR binary included with
    the installation.

    :return: path to the binary
    :rtype: str
    """
    bin_filename = 'STAR.exe' if config.PLATFORM == 'windows' else 'STAR'
    path = os.path.join(config.BINS_DIR, config.PLATFORM, 'STAR', bin_filename)
    if not os.path.exists(path):
        raise UnsupportedOSException(f'This operating system ({config.PLATFORM}) is not supported.')
    return path


def get_STAR_version():
    """Get the provided STAR version.

    :return: version string
    :rtype: str
    """
    p = run_executable([get_STAR_binary_path(), '--version'], quiet=True, returncode=1)
    version = p.stdout.read().strip()
    return version


def combine_arguments(args, additional):
    """Combine two dictionaries representing command-line arguments.

    Any duplicate keys will be merged according to the following procedure:
    1. If the value in both dictionaries are lists, the two lists are combined.
    2. Otherwise, the value in the first dictionary is OVERWRITTEN.

    :param args: original command-line arguments
    :type args: dictionary
    :param additional: additional command-line arguments
    :type additional: dictionary

    :return: combined command-line arguments
    :rtype: dictionary
    """
    new_args = args.copy()

    for key, value in additional.items():
        if key in new_args:
            if isinstance(value, list) and isinstance(new_args[key], list):
                new_args[key] += value
            else:
                new_args[key] = value
        else:
            new_args[key] = value

    return new_args


def arguments_to_list(args):
    """Convert a dictionary of command-line arguments to a list.

    :param args: command-line arguments
    :type args: dictionary

    :return: list of command-line arguments
    :rtype: list
    """
    arguments = []
    for key, value in args.items():
        arguments.append(key)
        if isinstance(value, list):
            arguments.extend(value)
        else:
            arguments.append(value)
    return arguments


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


def make_pool_with_counter(n_threads):
    """Create a new Process pool with a shared progress counter.

    :param n_threads: number of processes
    :type n_threads: int

    :return: (Process pool, progress counter, lock)
    :rtype: (multiprocessing.Pool, multiprocessing.Value, multiprocessing.Lock)
    """
    manager = multiprocessing.Manager()
    counter = manager.Value('I', 0)
    lock = manager.Lock()
    pool = multiprocessing.Pool(n_threads)
    return pool, counter, lock


def display_progress_with_counter(counter, total, *async_results, desc=None):
    """Display TQDM progress bar for displaying multiprocessing progress.

    :param counter: progress counter
    :type counter: multiprocessing.Value
    :param total: maximum number of units of processing
    :type total: int
    :param *async_results: multiprocessing results to monitor. These are used to
                           determine when all processes are done.
    :type *async_results: multiprocessing.pool.AsyncResult
    :param desc: progress bar description, defaults to `None`
    :type desc: str, optional
    """
    with tqdm(total=total, ascii=True, unit_scale=True, smoothing=0.1, desc=desc) as pbar:
        previous_progress = 0
        while any(not async_result.ready() for async_result in async_results):
            time.sleep(0.1)
            progress = counter.value
            pbar.update(progress - previous_progress)
            previous_progress = progress


def as_completed_with_progress(futures):
    """Wrapper around `concurrent.futures.as_completed` that displays a progress bar.

    :param futures: iterator of `concurrent.futures.Future` objects
    :type futures: iterable
    """
    with tqdm(total=len(futures), ascii=True, smoothing=0.1) as pbar:
        for future in as_completed(futures):
            yield future
            pbar.update(1)


def flatten_dictionary(d, keys=None):
    """Generator that flattens the given dictionary into 2-element tuples
    containing keys and values. For nested dictionaries, the keys are
    appended into a tuple.

    :param d: dictionary to flatten
    :type d: dictionary
    :param keys: previous keys, defaults to `None`. Used exclusively for recursion.
    :type keys: tuple, optional

    :return: flattened dictionary as (keys, value)
    :rtype: (tuple, object)
    """
    keys = keys or tuple()
    for k, v in d.items():
        new_keys = keys + (k,)
        if isinstance(v, dict):
            yield from flatten_dictionary(v, new_keys)
        else:
            yield new_keys, v


def flatten_list(lst):
    """Generator that flattens the given list.

    :param lst: list to flatten
    :type lst: list

    :return: flattened list
    :rtype: list
    """
    if not isinstance(lst, list):
        yield lst

    for element in lst:
        if not isinstance(element, list):
            yield element
        else:
            yield from flatten_list(element)


def merge_dictionaries(d1, d2, f=add, default=0):
    """Merge two dictionaries, applying an arbitrary function `f` to duplicate keys.
    Dictionaries may be nested.

    :param d1: first dictionary
    :type d1: dictionary
    :param d2: second dictionary
    :type d2: dictionary
    :param f: merge function. This function should take two arguments and return one,
              defaults to `+`
    :type f: callable, optional
    :param default: default value or callable to use for keys not present in either
                    dictionary, defaults to `0`
    :type default: object, optional

    :return: merged dictionary
    :rtype: dictionary
    """

    def setdefault_nested(d, t, value):
        inner = d
        for k in t[:-1]:
            inner = inner.setdefault(k, {})
        return inner.setdefault(t[-1], value)

    def get_nested(d, t, default=None):
        inner = d
        for k in t[:-1]:
            inner = inner.get(k, {})
        return inner.get(t[-1], default() if callable(default) else default)

    # Extract all keys
    d1_keys = [key for key, value in flatten_dictionary(d1)]
    d2_keys = [key for key, value in flatten_dictionary(d2)]
    keys = list(set(d1_keys + d2_keys))

    merged = {}
    for key in sorted(keys):
        setdefault_nested(merged, key, f(get_nested(d1, key, default), get_nested(d2, key, default)))

    return merged


def write_pickle(obj, path, *args, **kwargs):
    """Pickle a Python object and compress with Gzip.

    Any additional arguments and keyword arguments are passed to `pickle.dump`.

    :param obj: object to pickle
    :type obj: object
    :param path: path to save pickle
    :type path: str

    :return: saved pickle path
    :rtype: str
    """
    with gzip.open(path, 'wb') as f:
        pickle.dump(obj, f, *args, **kwargs)
    return path


def read_pickle(path):
    """Load a Python pickle that was compressed with Gzip.

    :param path: path to pickle
    :type path: str

    :return: unpickled object
    :rtype: object
    """
    with gzip.open(path, 'rb') as f:
        return pickle.load(f)


def split_index(index, n=8):
    """Split a conversions index, which is a list of tuples (file position,
    number of lines), one for each read, into `n` approximately equal parts.
    This function is used to split the conversions CSV for multiprocessing.

    :param index: index
    :type index: list
    :param n: number of splits, defaults to `8`
    :type n: int, optional

    :return: list of parts, where each part is a (file position, number of lines) tuple
    :rtype: list
    """
    n_lines = sum(idx[1] for idx in index)
    target = (n_lines // n) + 1  # add one to prevent underflow

    # Split the index to "approximately" equal parts
    parts = []
    start_pos = None
    current_size = 0
    for pos, size in index:
        if start_pos is None:
            start_pos = pos
        current_size += size

        if current_size >= target:
            parts.append((start_pos, current_size))
            start_pos = None
            current_size = 0
    if current_size > 0:
        parts.append((start_pos, current_size))

    return parts
