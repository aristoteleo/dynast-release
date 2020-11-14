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
            '--outSAMmultNmax',
            1,
            '--outFilterScoreMinOverLread',
            0.3,
            '--outFilterMatchNminOverLread',
            0.3,
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
            '--outSAMmultNmax',
            1,
            '--outFilterScoreMinOverLread',
            0.3,
            '--outFilterMatchNminOverLread',
            0.3,
        ], STAR_SOLO_BAM_TAGS, os.path.join(WHITELIST_DIR, 'scifate.txt.gz')
    ),
    Technology(
        '10xv2', [
            '--soloType',
            'CB_UMI_Simple',
            '--soloCBstart',
            1,
            '--soloCBlen',
            16,
            '--soloUMIstart',
            17,
            '--soloUMIlen',
            10,
            '--soloFeatures',
            'Gene',
            'Velocyto',
            '--outSAMmultNmax',
            1,
            '--outFilterScoreMinOverLread',
            0.3,
            '--outFilterMatchNminOverLread',
            0.3,
        ], STAR_SOLO_BAM_TAGS, os.path.join(WHITELIST_DIR, '10xv2.txt.gz')
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
            '--outSAMmultNmax',
            1,
            '--outFilterScoreMinOverLread',
            0.3,
            '--outFilterMatchNminOverLread',
            0.3,
        ], ['RG'], None
    )
]

TECHNOLOGIES = BARCODE_UMI_TECHNOLOGIES + PLATE_TECHNOLOGIES
TECHNOLOGIES_MAP = {t.name: t for t in TECHNOLOGIES}
