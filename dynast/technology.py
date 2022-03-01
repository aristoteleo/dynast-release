from collections import namedtuple

import ngs_tools as ngs

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
