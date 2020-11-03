import gzip
import pickle


def write_index(index, index_path):
    with gzip.open(index_path, 'wb') as f:
        pickle.dump(index, f, protocol=4)
    return index_path


def read_index(index_path):
    with gzip.open(index_path, 'rb') as f:
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
