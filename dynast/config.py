import os
import platform

from .utils import run_executable

PACKAGE_PATH = os.path.dirname(__file__)
PLATFORM = platform.system().lower()
BINS_DIR = os.path.join(PACKAGE_PATH, 'bins')
WHITELIST_DIR = os.path.join(PACKAGE_PATH, 'whitelists')
MODELS_DIR = os.path.join(PACKAGE_PATH, 'models')
MODEL_PATH = os.path.join(MODELS_DIR, 'pi.stan')
MODEL_NAME = 'pi'

RE_CHOICES = ['align', 'index', 'parse', 'snp', 'count', 'aggregate', 'p', 'pi', 'split']
GROUP_CHOICES = ['barcode', 'GX']
RECOMMENDED_MEMORY = 16 * (1024**3)  # 16 GB

# Common arguments for all STAR runs
STAR_ARGUMENTS = {
    '--outSAMmultNmax': 1,
    '--outSAMtype': ['BAM', 'SortedByCoordinate'],
    '--outSAMattributes': ['NH', 'HI', 'AS', 'NM', 'nM', 'MD', 'GX', 'GN'],
    '--bamRemoveDuplicatesType': 'UniqueIdentical',
    '--outFilterScoreMinOverLread': 0.3,
    '--outFilterMatchNminOverLread': 0.3,
}

# Additional arguments for STARsolo runs
STAR_SOLO_ARGUMENTS = {
    '--outSAMattributes': ['CR', 'CY', 'UR', 'UY', 'CB', 'UB', 'sS', 'sQ', 'sM'],
    '--soloFeatures': ['Gene', 'Velocyto'],
    '--soloCBwhitelist': 'None',
}

# Additional arguments for NASC-seq
NASC_ARGUMENTS = {
    '--alignSJoverhangMin': 1000,
    '--alignSJDBoverhangMin': 1,
    '--outFilterMismatchNoverReadLmax': 1,
    '--outFilterMismatchNmax': 10,
    '--outFilterMismatchNoverLmax': 0.1,
    '--outFilterScoreMinOverLread': 0.66,
    '--outFilterMatchNminOverLread': 0.66,
    '--scoreDelOpen': -10000,
    '--scoreInsOpen': -10000,
}


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
    """
    arguments = []
    for key, value in args.items():
        arguments.append(key)
        if isinstance(value, list):
            arguments.extend(value)
        else:
            arguments.append(value)
    return arguments
