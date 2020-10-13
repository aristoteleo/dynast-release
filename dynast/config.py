import os
import platform

PACKAGE_PATH = os.path.dirname(__file__)
PLATFORM = platform.system().lower()
BINS_DIR = 'bins'

RE_CHOICES = ['align', 'parse', 'count', 'aggregate', 'estimate']
RECOMMENDED_MEMORY = 16 * (1024**3)  # 16 GB
STAR_SOLO_OPTIONS = [
    '--outSAMattributes',
    'NH',
    'HI',
    'AS',
    'NM',
    'nM',
    'MD',
    'CR',
    'CY',
    'UR',
    'UY',
    'GX',
    'GN',
    'CB',
    'UB',
    'sM',
    'sS',
    'sQ',
    '--outSAMtype',
    'BAM',
    'SortedByCoordinate',
    '--outSAMmultNmax',
    1,
    '--soloFeatures',
    'Gene',
    'Velocyto',
]


class UnsupportedOSException(Exception):
    pass


def get_STAR_binary_path():
    """Get the path to the platform-dependent STAR binary included with
    the installation.
    :return: path to the binary
    :rtype: str
    """
    bin_filename = 'STAR.exe' if PLATFORM == 'windows' else 'STAR'
    path = os.path.join(PACKAGE_PATH, BINS_DIR, PLATFORM, 'STAR', bin_filename)
    if not os.path.exists(path):
        raise UnsupportedOSException(f'This operating system ({PLATFORM}) is not supported.')
    return path
