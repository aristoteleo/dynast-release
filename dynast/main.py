import argparse
import logging
import os
import shutil
import sys
import warnings

from . import __version__
from .config import get_STAR_binary_path, RE_CHOICES
from .count import count
from .preprocessing import CONVERSION_COLUMNS
from .ref import ref
from .technology import TECHNOLOGIES_MAP

logger = logging.getLogger(__name__)


def setup_ref_args(parser, parent):
    """Helper function to set up a subparser for the `ref` command.

    :param parser: argparse parser to add the `ref` command to
    :type args: argparse.ArgumentParser
    :param parent: argparse parser parent of the newly added subcommand.
                   used to inherit shared commands/flags
    :type args: argparse.ArgumentParser

    :return: the newly added parser
    :rtype: argparse.ArgumentParser
    """
    parser_ref = parser.add_parser(
        'ref',
        description='Build a STAR index from a reference',
        help='Build a STAR index from a reference',
        parents=[parent],
    )
    parser_ref._actions[0].help = parser_ref._actions[0].help.capitalize()

    required_ref = parser_ref.add_argument_group('required arguments')
    required_ref.add_argument(
        '-i',
        metavar='INDEX',
        help='Path to the directory where the STAR index will be generated',
        type=str,
        required=True
    )
    parser_ref.add_argument(
        '-m', metavar='MEMORY', help='Maximum memory used, in GB (default: 16)', type=int, default=16
    )
    parser_ref.add_argument(
        'fasta',
        help='Genomic FASTA file',
        type=str,
    )
    parser_ref.add_argument(
        'gtf',
        help='Reference GTF file',
        type=str,
    )

    return parser_ref


def setup_count_args(parser, parent):
    """Helper function to set up a subparser for the `count` command.

    :param parser: argparse parser to add the `count` command to
    :type args: argparse.ArgumentParser
    :param parent: argparse parser parent of the newly added subcommand.
                   used to inherit shared commands/flags
    :type args: argparse.ArgumentParser

    :return: the newly added parser
    :rtype: argparse.ArgumentParser
    """
    parser_count = parser.add_parser(
        'count',
        description='Quantify new RNA',
        help='Quantify new RNA',
        parents=[parent],
    )
    parser_count._actions[0].help = parser_count._actions[0].help.capitalize()

    required_count = parser_count.add_argument_group('required arguments')
    required_count.add_argument(
        '-i', metavar='INDEX', help='Path to the directory where the STAR index is located', type=str, required=True
    )
    parser_count.add_argument(
        '-o',
        metavar='OUT',
        help='Path to output directory (default: current directory)',
        type=str,
        default='.',
    )
    parser_count.add_argument(
        '-x',
        metavar='TECHNOLOGY',
        help=f'Single-cell technology used. Available choices are: {",".join(TECHNOLOGIES_MAP.keys())}',
        type=str,
        required=True,
        choices=TECHNOLOGIES_MAP.keys()
    )
    parser_count.add_argument(
        '-w',
        metavar='WHITELIST',
        help=('Path to file of whitelisted barcodes to correct to. '
              'If not provided, all barcodes are used.'),
        type=str,
        default=None
    )
    parser_count.add_argument(
        '--quality',
        metavar='QUALITY',
        help=(
            'Base quality threshold. Only bases with PHRED quality greater than '
            'this value will be considered when counting conversions. (default: 27)'
        ),
        type=int,
        default=27
    )
    parser_count.add_argument(
        '--re',
        metavar='RE',
        help=(
            'Re-do a step in the pipeline. Available choices are: `align`, '
            '`parse`, `count`, `aggregate`, `estimate`.'
        ),
        type=str,
        choices=RE_CHOICES,
        default=None
    )
    parser_count.add_argument(
        '--use-corrected-barcodes',
        help=(
            'Use STAR corrected barcodes instead of raw barcodes. This corresponds to '
            'the `CB` tags in the STAR alignments BAM, as opposed to the `CR` tags. '
            'All `barcode` columns in CSV files will correspond to '
            'corrected barcodes instead of raw barcodes. This option is only '
            'available if a whitelist is provided with the `-w` argument and is '
            'highly recommended.'
        ),
        action='store_true'
    )
    parser_count.add_argument(
        '--p-group-by',
        metavar='GROUPBY',
        help=(
            'Comma-delimited column names to group by when calculating p_e and p_c estimates. '
            'Available choices are: `barcode`, `GX`. `barcode` corresponds to '
            'either raw or corrected barcodes, depending on whether '
            '--use-corrected-barcodes was specified. `GX` corresponds to genes. '
            '(default: `barcode`)'
        ),
        type=str,
        default='barcode',
    )
    parser_count.add_argument(
        '--pi-group-by',
        metavar='GROUPBY',
        help=('Same as `--p-group-by`, but for pi estimation. (default: `barcode,GX`)'
              '(default: `barcode,GX`)'),
        type=str,
        default='barcode,GX',
    )
    parser_count.add_argument(
        '--conversion',
        metavar='CONVERSION',
        help=argparse.SUPPRESS,
        type=str,
        choices=CONVERSION_COLUMNS,
        default='TC',
    )
    parser_count.add_argument(
        'fastqs',
        help=(
            'FASTQ files. If `-x smartseq`, this is a single manifest CSV file where '
            'the first column contains cell IDs and the next two columns contain '
            'paths to FASTQs (the third column may contain a dash `-` for single-end reads).'
        ),
        nargs='+'
    )

    return parser_count


def parse_ref(parser, args, temp_dir=None):
    """Parser for the `ref` command.
    :param args: Command-line arguments dictionary, as parsed by argparse
    :type args: dict
    """
    if os.path.exists(args.i):
        parser.error(
            f'STAR index directory {args.i} already exists. '
            'Please provide a different directory or remove the existing one.'
        )
    ref(args.fasta, args.gtf, args.i, n_threads=args.t, memory=args.m * (1024**3), temp_dir=temp_dir)


def parse_count(parser, args, temp_dir=None):
    """Parser for the `count` command.
    :param args: Command-line arguments dictionary, as parsed by argparse
    :type args: dict
    """
    # Whitelist must be provided if we want to use corrected barcodes
    if args.use_corrected_barcodes and args.w is None:
        parser.error('Whitelist must be provided for `--use-corrected-barcodes`')
    if args.w is not None and not args.use_corrected_barcodes:
        logger.warning('`--use-corrected-barcodes` is highly recommended when providing a whitelist')

    # Check quality
    if args.quality < 0 or args.quality > 41:
        parser.error('`--quality` must be in [0, 42)')

    # Check group by
    if args.p_group_by.lower() == 'none':
        args.p_group_by = None
    else:
        args.p_group_by = args.p_group_by.split(',')
    if args.pi_group_by.lower() == 'none':
        args.pi_group_by = None
    else:
        args.pi_group_by = args.pi_group_by.split(',')
    if args.pi_group_by is None and args.p_group_by is not None:
        parser.error('`--p-group-by` must be `None` if `--pi-group-by` is `None`')
    if args.p_group_by is not None and args.pi_group_by is not None and not set(args.p_group_by).issubset(set(
            args.pi_group_by)):
        parser.error('`--p-group-by` must be a subset of `--pi-group-by`')

    technology = TECHNOLOGIES_MAP[args.x]
    if technology.name == 'smartseq':
        if len(args.fastqs) != 1:
            parser.error('A single manifest TSV must be provided for technology `smartseq`')
        with open(args.fastqs[0], 'r') as f:
            cell_ids = []
            for line in f:
                cell_id, fastq_1, fastq_2 = line.strip().split(',')
                if cell_id in cell_ids:
                    parser.error(f'Found duplicate cell ID in manifest CSV: {cell_id}')
                if not os.path.exists(fastq_1):
                    parser.error(
                        f'{fastq_1} does not exist. '
                        f'All paths in {args.fastqs[0]} must be relative to the '
                        'current working directory or absolute.'
                    )
                if fastq_2 != '-' and not os.path.exists(fastq_2):
                    parser.error(
                        f'{fastq_2} does not exist. '
                        f'All paths in {args.fastqs[0]} must be relative to the '
                        'current working directory or absolute.'
                    )
                cell_ids.append(cell_id)

    count(
        args.fastqs,
        args.i,
        args.o,
        TECHNOLOGIES_MAP[args.x],
        use_corrected=args.use_corrected_barcodes,
        quality=args.quality,
        conversion=args.conversion,
        p_group_by=args.p_group_by,
        pi_group_by=args.pi_group_by,
        whitelist_path=args.w,
        n_threads=args.t,
        re=args.re,
        temp_dir=temp_dir,
    )


COMMAND_TO_FUNCTION = {
    'ref': parse_ref,
    'count': parse_count,
}


def main():
    parser = argparse.ArgumentParser(description=f'{__version__}')
    parser._actions[0].help = parser._actions[0].help.capitalize()
    subparsers = parser.add_subparsers(
        dest='command',
        metavar='<CMD>',
    )

    # Add common options to this parent parser
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument('--tmp', metavar='TMP', help='Override default temporary directory', type=str, default='tmp')
    parent.add_argument('--keep-tmp', help='Do not delete the tmp directory', action='store_true')
    parent.add_argument('--verbose', help='Print debugging information', action='store_true')
    parent.add_argument('-t', metavar='THREADS', help='Number of threads to use (default: 8)', type=int, default=8)

    # Command parsers
    parser_ref = setup_ref_args(subparsers, parent)
    parser_count = setup_count_args(subparsers, parent)
    command_to_parser = {
        'ref': parser_ref,
        'count': parser_count,
    }

    # Show help when no arguments are given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if len(sys.argv) == 2:
        if sys.argv[1] in command_to_parser:
            command_to_parser[sys.argv[1]].print_help(sys.stderr)
        else:
            parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    logging.basicConfig(
        format='[%(asctime)s] %(levelname)7s %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO,
        force=True,
    )
    logging.getLogger('numba').setLevel(logging.CRITICAL)
    logging.getLogger('pystan').setLevel(logging.CRITICAL)
    logging.getLogger('anndata').setLevel(logging.CRITICAL)
    logger.debug('Printing verbose output')
    logger.debug(f'Input args: {args}')
    logger.debug(f'STAR binary located at {get_STAR_binary_path()}')
    logger.debug(f'Creating {args.tmp} directory')
    if os.path.exists(args.tmp):
        parser.error(
            f'Temporary directory {args.tmp} already exists. '
            'Is another process running? Please specify a different temporary '
            'directory with the `--tmp` option, or remove the one that already '
            'exists.'
        )
    os.makedirs(args.tmp)
    os.environ['NUMEXPR_MAX_THREADS'] = str(args.t)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            COMMAND_TO_FUNCTION[args.command](parser, args, temp_dir=args.tmp)
    except Exception:
        logger.exception('An exception occurred')
    finally:
        if not args.keep_tmp:
            logger.debug(f'Removing {args.tmp} directory')
            shutil.rmtree(args.tmp, ignore_errors=True)
