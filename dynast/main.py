import argparse
import logging
import os
import shutil
import sys

from . import __version__
from .config import get_STAR_binary_path, RE_CHOICES
from .count import count
from .ref import ref


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
        description='',
        help='',
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
            'corrected barcodes instead of raw barcodes.'
        ),
        action='store_true'
    )
    parser_count.add_argument(
        '--group-by',
        metavar='GROUPBY',
        help=(
            'Column name to group by when calculating p_e and p_c estimates. '
            'Available choices are: `barcode`, `GX`. `barcode` corresponds to '
            'either raw or corrected barcodes, depending on whether '
            '--use-corrected-barcodes was specified. `GX` corresponds to genes.'
        ),
        type=str,
        choices=['barcode', 'GX'],
        default='barcode',
    )
    parser_count.add_argument('fastqs', help='FASTQ files', nargs='+')

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
    ref()


def parse_count(parser, args, temp_dir=None):
    """Parser for the `count` command.
    :param args: Command-line arguments dictionary, as parsed by argparse
    :type args: dict
    """
    # Check quality
    if args.quality < 0 or args.quality > 41:
        parser.error('`--quality` must be in [0, 42)')

    count(
        args.fastqs,
        args.i,
        args.o,
        use_corrected=args.use_corrected_barcodes,
        quality=args.quality,
        group_by=args.group_by,
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
    )
    logger = logging.getLogger(__name__)
    logger.debug('Printing verbose output')
    logger.debug(args)
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
        COMMAND_TO_FUNCTION[args.command](parser, args, temp_dir=args.tmp)
    except Exception:
        logger.exception('An exception occurred')
    finally:
        if not args.keep_tmp:
            logger.debug(f'Removing {args.tmp} directory')
            shutil.rmtree(args.tmp, ignore_errors=True)
