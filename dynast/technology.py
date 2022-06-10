import re
from collections import namedtuple

import ngs_tools as ngs
import pysam

Technology = namedtuple('Technology', ['name', 'chemistry', 'additional_args'])

BARCODE_UMI_TECHNOLOGIES = [
    Technology('dropseq',
               ngs.chemistry.get_chemistry('dropseq').reorder([1, 0]), None),
    Technology('scifate',
               ngs.chemistry.get_chemistry('scifate').reorder([1, 0]), None),
    Technology('10xv2',
               ngs.chemistry.get_chemistry('10xv2').reorder([1, 0]), None),
    Technology('10xv3',
               ngs.chemistry.get_chemistry('10xv3').reorder([1, 0]), None),
]

PLATE_TECHNOLOGIES = [
    Technology(
        'smartseq', ngs.chemistry.get_chemistry('smartseq2'), {
            '--soloUMIdedup': 'Exact',
            '--outSAMattributes': ['RG'],
        }
    ),
]

TECHNOLOGIES = BARCODE_UMI_TECHNOLOGIES + PLATE_TECHNOLOGIES
TECHNOLOGIES_MAP = {t.name: t for t in TECHNOLOGIES}
STRAND_MAP = {
    ngs.chemistry.SequencingStrand.UNSTRANDED: 'unstranded',
    ngs.chemistry.SequencingStrand.FORWARD: 'forward',
    ngs.chemistry.SequencingStrand.REVERSE: 'reverse',
}
BAM_STRAND_PARSER = re.compile('(--readStrand|--soloStrand)(=| +)(?P<strand>Forward|Reverse|Unstranded)')


def detect_strand(bam_path):
    """Attempt to detect strandness by parsing the BAM header.

    :param bam_path: Path to BAM
    :type bam_path: str

    :return: 'unstranded', 'forward', or 'reverse if the strand was successfully detected.
        `None` otherwise.
    :rtype: str, optional
    """
    with pysam.AlignmentFile(bam_path, 'rb') as f:
        pg = f.header.get('PG')
        if pg:
            for entry in pg:
                if entry.get('PN') == 'STAR':
                    search = BAM_STRAND_PARSER.search(entry.get('CL', ''))
                    if search:
                        return search.groupdict()['strand'].lower()
    return None
