import gzip
import logging
import os
import shutil
import subprocess as sp
import tempfile
import time
from contextlib import contextmanager

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
    """Get the maximum number of open file descriptors in a platform-dependent
    way.
    """
    if config.PLATFORM == 'windows':
        import win32file
        return win32file._getmaxstdio()
    else:
        import resource
        return resource.getrlimit(resource.RLIMIT_NOFILE)[0]


def get_max_file_descriptor_limit():
    """
    """
    if config.PLATFORM == 'windows':
        return 2048
    else:
        import resource
        return resource.getrlimit(resource.RLIMIT_NOFILE)[1]


@contextmanager
def increase_file_descriptor_limit(limit):
    """Context manager that can be used to temporarily increase the maximum
    number of open file descriptors for the current process.

    This is required when running STAR with many threads.
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


def mkstemp(dir=None):
    fd, path = tempfile.mkstemp(dir=dir)
    os.close(fd)
    return path
