import os
from collections import namedtuple

from .config import STAR_SOLO_BAM_TAGS, WHITELIST_DIR

Technology = namedtuple('Technology', ['name', 'STAR_args', 'bam_tags', 'whitelist_path'])

BARCODE_UMI_TECHNOLOGIES = [
    Technology(
        'dropseq', [
            '--soloType',
            'CB_UMI_Simple',
            '--soloCBstart',
            1,
            '--soloCBlen',
            12,
            '--soloUMIstart',
            13,
            '--soloUMIlen',
            8,
            '--soloFeatures',
            'Gene',
            'Velocyto',
        ], STAR_SOLO_BAM_TAGS, None
    ),
    Technology(
        'scifate', [
            '--soloType',
            'CB_UMI_Simple',
            '--soloCBstart',
            9,
            '--soloCBlen',
            10,
            '--soloUMIstart',
            1,
            '--soloUMIlen',
            8,
            '--soloFeatures',
            'Gene',
            'Velocyto',
        ], STAR_SOLO_BAM_TAGS, os.path.join(WHITELIST_DIR, 'scifate.txt.gz')
    ),
]

PLATE_TECHNOLOGIES = [
    Technology(
        'smartseq', [
            '--soloType',
            'SmartSeq',
            '--soloUMIdedup',
            'Exact',
            '--soloStrand',
            'Unstranded',
            '--soloFeatures',
            'Gene',
        ], ['RG'], None
    )
]

TECHNOLOGIES = BARCODE_UMI_TECHNOLOGIES + PLATE_TECHNOLOGIES
TECHNOLOGIES_MAP = {t.name: t for t in TECHNOLOGIES}
