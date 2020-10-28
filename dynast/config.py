import os
import platform

PACKAGE_PATH = os.path.dirname(__file__)
PLATFORM = platform.system().lower()
BINS_DIR = os.path.join(PACKAGE_PATH, 'bins')
WHITELIST_DIR = os.path.join(PACKAGE_PATH, 'whitelists')
MODELS_DIR = os.path.join(PACKAGE_PATH, 'models')
MODEL_PATH = os.path.join(MODELS_DIR, 'pi.stan')
MODEL_NAME = 'pi'

RE_CHOICES = ['align', 'parse', 'snp', 'count', 'aggregate', 'p', 'pi', 'split']
GROUP_CHOICES = ['barcode', 'GX']
RECOMMENDED_MEMORY = 16 * (1024**3)  # 16 GB
STAR_SOLO_OPTIONS = [
    '--outSAMtype',
    'BAM',
    'SortedByCoordinate',
    '--outSAMmultNmax',
    1,
    '--outFilterScoreMinOverLread',
    0.3,
    '--outFilterMatchNminOverLread',
    0.3,
    '--outSAMattributes',
    'NH',
    'HI',
    'AS',
    'NM',
    'nM',
    'MD',
    'GX',
    'GN',
]
STAR_SOLO_BAM_TAGS = ['CR', 'CY', 'UR', 'UY', 'CB', 'UB', 'sS', 'sQ', 'sM']


class UnsupportedOSException(Exception):
    pass


def get_STAR_binary_path():
    """Get the path to the platform-dependent STAR binary included with
    the installation.

    :return: path to the binary
    :rtype: str
    """
    bin_filename = 'STAR.exe' if PLATFORM == 'windows' else 'STAR'
    path = os.path.join(BINS_DIR, PLATFORM, 'STAR', bin_filename)
    if not os.path.exists(path):
        raise UnsupportedOSException(f'This operating system ({PLATFORM}) is not supported.')
    return path
