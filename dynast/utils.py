import concurrent
import functools
import multiprocessing
import os
import struct
import sys
import time
from concurrent.futures import as_completed
from contextlib import contextmanager
from typing import Any, Dict, FrozenSet, Iterable, List, Optional, Tuple, Union

import anndata
import ngs_tools as ngs
import numpy as np
import pandas as pd
import psutil
from scipy import sparse

from . import config
from .logging import logger

# As of 0.0.1, these are provided by ngs_tools, but keep these here for now because
# they are imported from this file in many places.
run_executable = ngs.utils.run_executable
open_as_text = ngs.utils.open_as_text
decompress_gzip = ngs.utils.decompress_gzip
flatten_dict_values = ngs.utils.flatten_dict_values
mkstemp = ngs.utils.mkstemp
all_exists = ngs.utils.all_exists
flatten_dictionary = ngs.utils.flatten_dictionary
flatten_iter = ngs.utils.flatten_iter
merge_dictionaries = ngs.utils.merge_dictionaries
write_pickle = ngs.utils.write_pickle
read_pickle = ngs.utils.read_pickle


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


def get_STAR_binary_path() -> str:
    """Get the path to the platform-dependent STAR binary included with
    the installation.

    Returns:
        Path to the binary
    """
    bin_filename = 'STAR.exe' if config.PLATFORM == 'windows' else 'STAR'
    path = os.path.join(config.BINS_DIR, config.PLATFORM, 'STAR', bin_filename)
    if not os.path.exists(path):
        raise UnsupportedOSException(f'This operating system ({config.PLATFORM}) is not supported.')
    return path


def get_STAR_version() -> str:
    """Get the provided STAR version.

    Returns:
        Version string
    """
    p, stdout, stderr = run_executable([get_STAR_binary_path(), '--version'], quiet=True, returncode=1)
    version = stdout.strip()
    return version


def combine_arguments(args: Dict[str, Any], additional: Dict[str, Any]) -> Dict[str, Any]:
    """Combine two dictionaries representing command-line arguments.

    Any duplicate keys will be merged according to the following procedure:
    1. If the value in both dictionaries are lists, the two lists are combined.
    2. Otherwise, the value in the first dictionary is OVERWRITTEN.

    Args:
        args: Original command-line arguments
        additional: Additional command-line arguments

    Returns:
        Combined command-line arguments
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


def arguments_to_list(args: Dict[str, Any]) -> List[Any]:
    """Convert a dictionary of command-line arguments to a list.

    Args:
        args: Command-line arguments

    Returns:
        List of command-line arguments
    """
    arguments = []
    for key, value in args.items():
        arguments.append(key)
        if isinstance(value, list):
            arguments.extend(value)
        else:
            arguments.append(value)
    return arguments


def get_file_descriptor_limit() -> int:
    """Get the current value for the maximum number of open file descriptors
    in a platform-dependent way.

    Returns:
        The current value of the maximum number of open file descriptors.
    """
    if config.PLATFORM == 'windows':
        import win32file
        return win32file._getmaxstdio()
    else:
        import resource
        return resource.getrlimit(resource.RLIMIT_NOFILE)[0]


def get_max_file_descriptor_limit() -> int:
    """Get the maximum allowed value for the maximum number of open file
    descriptors.

    Note that for Windows, there is not an easy way to get this,
    as it requires reading from the registry. So, we just return the maximum for
    a vanilla Windows installation, which is 8192.
    https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/setmaxstdio?view=vs-2019

    Similarly, on MacOS, we return a hardcoded 10240.

    Returns:
        Maximum allowed value for the maximum number of open file descriptors
    """
    if config.PLATFORM == 'windows':
        return 8192
    elif config.PLATFORM == 'darwin':
        return 10240
    else:
        import resource
        return resource.getrlimit(resource.RLIMIT_NOFILE)[1]


@contextmanager
def increase_file_descriptor_limit(limit: int):
    """Context manager that can be used to temporarily increase the maximum
    number of open file descriptors for the current process. The original
    value is restored when execution exits this function.

    This is required when running STAR with many threads.

    Args:
        limit: Maximum number of open file descriptors will be increased to
            this value for the duration of the context
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


def get_available_memory() -> int:
    """Get total amount of available memory (total memory - used memory) in bytes.

    Returns:
        Available memory in bytes
    """
    return psutil.virtual_memory().available


def make_pool_with_counter(n_threads: int) -> Tuple[multiprocessing.Pool, multiprocessing.Value, multiprocessing.Lock]:
    """Create a new Process pool with a shared progress counter.

    Args:
        n_threads: Number of processes

    Returns:
        Tuple of (Process pool, progress counter, lock)
    """
    manager = multiprocessing.Manager()
    counter = manager.Value('I', 0)
    lock = manager.Lock()
    pool = multiprocessing.Pool(n_threads)
    return pool, counter, lock


def display_progress_with_counter(
    counter: multiprocessing.Value, total: int, *async_results, desc: Optional[str] = None
):
    """Display progress bar for displaying multiprocessing progress.

    Args:
        counter: Progress counter
        total: Maximum number of units of processing
        *async_results: Multiprocessing results to monitor. These are used to
            determine when all processes are done.
        desc: Progress bar description
    """
    with ngs.progress.progress(total=total, unit_scale=True, desc=desc) as pbar:
        previous_progress = 0
        while any(not async_result.ready() for async_result in async_results):
            time.sleep(0.01)
            progress = counter.value
            pbar.update(progress - previous_progress)
            pbar.refresh()
            previous_progress = progress


def as_completed_with_progress(futures: Iterable[concurrent.futures.Future]):
    """Wrapper around `concurrent.futures.as_completed` that displays a progress bar.

    Args:
        Iterator of `concurrent.futures.Future` objects
    """
    with ngs.progress.progress(total=len(futures)) as pbar:
        for future in as_completed(futures):
            yield future
            pbar.update(1)


def split_index(index: List[Tuple[int, int, int]], n: int = 8) -> List[List[Tuple[int, int, int]]]:
    """Split a conversions index, which is a list of tuples (file position,
    number of lines, alignment position), one for each read, into `n`
    approximately equal parts. This function is used to split the conversions
    CSV for multiprocessing.

    Args:
        index: index
        n: Number of splits, defaults to `8`

    Returns:
        List of parts
    """
    n_lines = sum(idx[1] for idx in index)
    target = (n_lines // n) + 1  # add one to prevent underflow

    # Split the index to "approximately" equal parts
    parts = []
    current_part = []
    current_size = 0
    for tup in index:
        current_part.append(tup)
        current_size += tup[1]

        if current_size >= target:
            parts.append(current_part)
            current_size = 0
            current_part = []
    if current_part:
        parts.append(current_part)

    return parts


def downsample_counts(
    df_counts: pd.DataFrame,
    proportion: Optional[float] = None,
    count: Optional[int] = None,
    seed: Optional[int] = None,
    group_by: Optional[List[str]] = None
) -> pd.DataFrame:
    """Downsample the given counts dataframe according to the ``proportion`` or
    ``count`` arguments. One of these two must be provided, but not both. The dataframe
    is assumed to be UMI-deduplicated.

    Args:
        df_counts: Counts dataframe
        proportion: Proportion of reads (UMIs) to keep
        count: Absolute number of reads (UMIs) to keep
        seed: Random seed
        group_by: Columns in the counts dataframe to use to group entries.
            When this is provided, UMIs are no longer sampled at random, but instead
            grouped by this argument, and only groups that have more than ``count`` UMIs
            are downsampled.

    Returns:
        Downsampled counts dataframe
    """
    rng = np.random.default_rng(seed)
    if not group_by:
        if bool(proportion) == bool(count):
            raise Exception('Only one of `proportion` or `count` must be provided.')

        n_keep = int(df_counts.shape[0] * proportion) if proportion is not None else count
        return df_counts.iloc[rng.choice(df_counts.shape[0], n_keep, shuffle=False, replace=False)]
    else:
        if not count:
            raise Exception('`count` must be provided when using `group_by`')

        dfs = []
        for key, df_group in df_counts.groupby(group_by, sort=False, observed=True):
            if df_group.shape[0] > count:
                df_group = df_group.iloc[rng.choice(df_group.shape[0], count, shuffle=False, replace=False)]
            dfs.append(df_group)
        return pd.concat(dfs)


def counts_to_matrix(
    df_counts: pd.DataFrame,
    barcodes: List[str],
    features: List[str],
    barcode_column: str = 'barcode',
    feature_column: str = 'GX'
) -> sparse.csr_matrix:
    """Convert a counts dataframe to a sparse counts matrix.

    Counts are assumed to be appropriately deduplicated.

    Args:
        df_counts: Counts dataframe
        barcodes: List of barcodes that will map to the rows
        features: List of features (i.e. genes) that will map to the columns
        barcode_column: Column in counts dataframe to use as barcodes, defaults to `barcode`
        feature_column: Column in counts dataframe to use as features, defaults to `GX`

    Returns:
        Sparse counts matrix
    """
    # Transform to index for fast lookup
    barcode_indices = {barcode: i for i, barcode in enumerate(barcodes)}
    feature_indices = {feature: i for i, feature in enumerate(features)}

    matrix = sparse.lil_matrix((len(barcodes), len(features)), dtype=np.float32)
    for (barcode, feature), count in df_counts.groupby([barcode_column, feature_column], sort=False,
                                                       observed=True).size().items():
        matrix[barcode_indices[barcode], feature_indices[feature]] = count

    return matrix.tocsr()


def split_counts(
    df_counts: pd.DataFrame,
    barcodes: List[str],
    features: List[str],
    barcode_column: str = 'barcode',
    feature_column: str = 'GX',
    conversions: FrozenSet[str] = frozenset({'TC'})
) -> Tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Split counts dataframe into two count matrices by a column.

    Args:
        df_counts: Counts dataframe
        barcodes: List of barcodes that will map to the rows
        features: List of features (i.e. genes) that will map to the columns
        barcode_column: Column in counts dataframe to use as barcodes
        feature_column: Column in counts dataframe to use as features
        conversions: Conversion(s) in question

    Returns:
        count matrix of `conversion==0`, count matrix of `conversion>0`
    """
    matrix_unlabeled = counts_to_matrix(
        df_counts[(df_counts[list(conversions)] == 0).all(axis=1)],
        barcodes,
        features,
        barcode_column=barcode_column,
        feature_column=feature_column
    )
    matrix_labeled = counts_to_matrix(
        df_counts[(df_counts[list(conversions)] > 0).any(axis=1)],
        barcodes,
        features,
        barcode_column=barcode_column,
        feature_column=feature_column
    )
    return matrix_unlabeled, matrix_labeled


def split_matrix(
    matrix: Union[np.ndarray, sparse.spmatrix], pis: Dict[Tuple[str, str], float], barcodes: List[str],
    features: List[str]
) -> Tuple[sparse.csr_matrix, sparse.csr_matrix, sparse.csr_matrix]:
    """Split the given matrix based on provided fraction of new RNA.

    Args:
        matrix: Matrix to split
        pis: Dictionary containing pi estimates
        barcodes: All barcodes
        features: All features (i.e. genes)

    Returns:
        matrix of pi masks, matrix of unlabeled RNA, matrix of labeled RNA
    """
    unlabeled_matrix = sparse.lil_matrix((len(barcodes), len(features)), dtype=np.float32)
    labeled_matrix = sparse.lil_matrix((len(barcodes), len(features)), dtype=np.float32)
    pi_mask = sparse.lil_matrix((len(barcodes), len(features)), dtype=bool)
    barcode_indices = {barcode: i for i, barcode in enumerate(barcodes)}
    feature_indices = {feature: i for i, feature in enumerate(features)}

    for (barcode, gx), pi in pis.items():
        try:
            pi = float(pi)
        except ValueError:
            continue
        row, col = barcode_indices[barcode], feature_indices[gx]
        val = matrix[row, col]
        unlabeled_matrix[row, col] = val * (1 - pi)
        labeled_matrix[row, col] = val * pi
        pi_mask[row, col] = True

    return pi_mask.tocsr(), unlabeled_matrix.tocsr(), labeled_matrix.tocsr()


def results_to_adata(
    df_counts: pd.DataFrame,
    conversions: FrozenSet[FrozenSet[str]] = frozenset({frozenset({'TC'})}),
    gene_infos: Optional[dict] = None,
    pis: Optional[Dict[Tuple[str, str], float]] = None
) -> anndata.AnnData:
    """Compile all results to a single anndata.

    Args:
        df_counts: Counts dataframe, with complemented reverse strand bases
        conversions: Conversion(s) in question
        gene_infos: Dictionary containing gene information. If this is not provided,
            the function assumes gene names are already in the Counts dataframe.
        pis: Dictionary of estimated pis

    Returns:
        Anndata containing all results
    """
    pis = pis or {}
    all_conversions = sorted(flatten_iter(conversions))
    transcriptome_exists = df_counts['transcriptome'].any()
    transcriptome_only = df_counts['transcriptome'].all()
    velocities = df_counts['velocity'].unique()
    barcodes = sorted(df_counts['barcode'].unique())
    features = sorted(df_counts['GX'].unique())

    df_counts_transcriptome = df_counts[df_counts['transcriptome']]
    matrix = counts_to_matrix(df_counts_transcriptome, barcodes, features)
    layers = {}

    # Transcriptome reads
    if transcriptome_exists:
        for convs in conversions:
            convs = sorted(convs)
            # Ignore reads that have other conversions
            other_convs = list(set(all_conversions) - set(convs))
            join = '_'.join(convs)
            # Counts for transcriptome reads (i.e. X_unlabeled + X_labeled = X)
            layers[f'X_n_{join}'], layers[f'X_l_{join}'] = split_counts(
                df_counts_transcriptome[(df_counts_transcriptome[other_convs] == 0).all(axis=1)],
                barcodes,
                features,
                conversions=convs
            )
            pi = pis.get('transcriptome', {}).get(tuple(convs))
            if pi is not None:
                (
                    _,
                    layers[f'X_n_{join}_est'],
                    layers[f'X_l_{join}_est'],
                ) = split_matrix(layers[f'X_n_{join}'] + layers[f'X_l_{join}'], pi, barcodes, features)
    else:
        logger.warning('No reads were assigned to `transcriptome`')

    # Total reads
    if not transcriptome_only:
        layers['total'] = counts_to_matrix(df_counts, barcodes, features)
        for convs in conversions:
            convs = sorted(convs)
            other_convs = list(set(all_conversions) - set(convs))
            join = '_'.join(convs)
            layers[f'unlabeled_{join}'], layers[f'labeled_{join}'] = split_counts(
                df_counts[(df_counts[other_convs] == 0).all(axis=1)], barcodes, features, conversions=convs
            )
            pi = pis.get('total', {}).get(tuple(convs))
            if pi is not None:
                (
                    _,
                    layers[f'unlabeled_{join}_est'],
                    layers[f'labeled_{join}_est'],
                ) = split_matrix(layers[f'unlabeled_{join}'] + layers[f'labeled_{join}'], pi, barcodes, features)

    # Velocity reads
    for key in velocities:
        if key == 'unassigned':
            continue
        df_counts_velocity = df_counts[df_counts['velocity'] == key]
        layers[key] = counts_to_matrix(df_counts_velocity, barcodes, features)
        if key in config.VELOCITY_BLACKLIST:
            continue

        for conv in conversions:
            other_convs = list(set(all_conversions) - set(convs))
            join = '_'.join(convs)
            layers[f'{key[0]}n_{join}'], layers[f'{key[0]}l_{join}'] = split_counts(
                df_counts_velocity[(df_counts_velocity[other_convs] == 0).all(axis=1)],
                barcodes,
                features,
                conversions=convs
            )
            pi = pis.get(key, {}).get(tuple(convs))
            if pi is not None:
                (
                    _,
                    layers[f'{key[0]}n_{join}_est'],
                    layers[f'{key[0]}l_{join}_est'],
                ) = split_matrix(layers[f'{key[0]}n_{join}'] + layers[f'{key[0]}l_{join}'], pi, barcodes, features)

    # Construct anndata
    add_names = gene_infos is not None
    var = pd.DataFrame(index=pd.Series(features, name='gene_id' if add_names else 'gene_name'))
    if add_names:
        var['gene_name'] = pd.Categorical([gene_infos.get(feature, {}).get('gene_name') for feature in features])
    return anndata.AnnData(
        X=matrix, obs=pd.DataFrame(index=pd.Series(barcodes, name='barcode')), var=var, layers=layers
    )


def patch_mp_connection_bpo_17560():
    """Apply PR-10305 / bpo-17560 connection send/receive max size update

    See the original issue at https://bugs.python.org/issue17560 and
    https://github.com/python/cpython/pull/10305 for the pull request.

    This only supports Python versions 3.3 - 3.7, this function
    does nothing for Python versions outside of that range.

    Taken from https://stackoverflow.com/a/47776649

    """
    patchname = "Multiprocessing connection patch for bpo-17560"
    if not (3, 3) < sys.version_info < (3, 8):
        return

    logger.debug(f'Applying {patchname}')

    from multiprocessing.connection import Connection

    orig_send_bytes = Connection._send_bytes
    orig_recv_bytes = Connection._recv_bytes
    if (orig_send_bytes.__code__.co_filename == __file__ and orig_recv_bytes.__code__.co_filename == __file__):
        logger.info(f'{patchname} already applied, skipping')
        return

    @functools.wraps(orig_send_bytes)
    def send_bytes(self, buf):
        n = len(buf)
        if n > 0x7fffffff:
            pre_header = struct.pack("!i", -1)
            header = struct.pack("!Q", n)
            self._send(pre_header)
            self._send(header)
            self._send(buf)
        else:
            orig_send_bytes(self, buf)

    @functools.wraps(orig_recv_bytes)
    def recv_bytes(self, maxsize=None):
        buf = self._recv(4)
        size, = struct.unpack("!i", buf.getvalue())
        if size == -1:
            buf = self._recv(8)
            size, = struct.unpack("!Q", buf.getvalue())
        if maxsize is not None and size > maxsize:
            return None
        return self._recv(size)

    Connection._send_bytes = send_bytes
    Connection._recv_bytes = recv_bytes


def collapse_anndata(adata: anndata.AnnData, by: Optional[str] = None) -> anndata.AnnData:
    """Collapse the given Anndata by summing duplicate rows. The `by` argument
    specifies which column to use. If not provided, the index is used.

    Note:
        This function also collapses any existing layers. Additionally, the
        returned AnnData will have the values used to collapse as the index.

    Args:
        adata: The Anndata to collapse
        by: The column to collapse by. If not provided, the index is used. When
            this column contains missing values (i.e. nan or None), these
            columns are removed.

    Returns:
        A new collapsed Anndata object. All matrices are sparse, regardless of
        whether or not they were in the input Anndata.
    """
    var = adata.var
    if by is not None:
        var = var.set_index(by)
    na_mask = var.index.isna()
    adata = adata[:, ~na_mask].copy()
    adata.var = var[~na_mask]

    if not any(adata.var.index.duplicated()):
        return adata

    var_indices = {}
    for i, index in enumerate(adata.var.index):
        var_indices.setdefault(index, []).append(i)

    # Convert all original matrices to csc for fast column operations
    X = sparse.csc_matrix(adata.X)
    layers = {layer: sparse.csc_matrix(adata.layers[layer]) for layer in adata.layers}
    new_index = []
    # lil_matrix is efficient for row-by-row construction
    new_X = sparse.lil_matrix((len(var_indices), adata.shape[0]))
    new_layers = {layer: new_X.copy() for layer in adata.layers}
    for i, (index, indices) in enumerate(var_indices.items()):
        new_index.append(index)
        new_X[i] = X[:, indices].sum(axis=1).flatten()
        for layer in layers.keys():
            new_layers[layer][i] = layers[layer][:, indices].sum(axis=1).flatten()

    return anndata.AnnData(
        X=new_X.T.tocsr(),
        layers={layer: new_layers[layer].T.tocsr()
                for layer in new_layers},
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=pd.Series(new_index, name=adata.var.index.name)),
    )
