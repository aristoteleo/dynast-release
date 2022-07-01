import os
from unittest import TestCase, mock

import dynast.align as align
import dynast.config as config
import dynast.utils as utils
from tests import mixins


class TestAlign(mixins.TestMixin, TestCase):

    def test_STAR_solo(self):
        with mock.patch('dynast.align.utils.get_STAR_binary_path', return_value='path/to/STAR'), \
            mock.patch('dynast.align.utils.arguments_to_list') as arguments_to_list, \
            mock.patch('dynast.align.utils.run_executable') as run_executable:
            result = align.STAR_solo(
                ['fastq1', 'fastq2'],
                'path/to/index',
                'output/path',
                self.umi_technology,
                whitelist_path=None,
                n_threads=8,
                temp_dir=self.temp_dir,
                nasc=False,
                overrides=None,
            )

            self.assertEqual({
                'bam': 'output/path/Aligned.sortedByCoord.out.bam',
                'gene': {
                    'raw': {
                        'barcodes': 'output/path/Solo.out/Gene/raw/barcodes.tsv',
                        'features': 'output/path/Solo.out/Gene/raw/features.tsv',
                        'matrix': 'output/path/Solo.out/Gene/raw/matrix.mtx'
                    },
                    'filtered': {
                        'barcodes': 'output/path/Solo.out/Gene/filtered/barcodes.tsv',
                        'features': 'output/path/Solo.out/Gene/filtered/features.tsv',
                        'matrix': 'output/path/Solo.out/Gene/filtered/matrix.mtx'
                    }
                }
            }, result)
            arguments = utils.combine_arguments(
                utils.combine_arguments(config.STAR_ARGUMENTS, self.umi_technology.chemistry.to_starsolo_arguments()),
                config.STAR_SOLO_ARGUMENTS
            )
            arguments.update({
                '--soloStrand': 'Forward',
                '--soloCBwhitelist': 'None',
                '--genomeDir': 'path/to/index',
                '--runThreadN': 8,
                '--outFileNamePrefix': 'output/path/',
                '--outTmpDir': mock.ANY,
                '--outBAMsortingBinsN': mock.ANY,
                '--readFilesIn': ['fastq1', 'fastq2']
            })
            arguments_to_list.assert_called_once_with(arguments)
            run_executable.assert_called_once()

    def test_STAR_solo_whitelist(self):
        with mock.patch('dynast.align.utils.get_STAR_binary_path', return_value='path/to/STAR'), \
            mock.patch('dynast.align.utils.arguments_to_list') as arguments_to_list, \
            mock.patch('dynast.align.utils.run_executable') as run_executable:
            result = align.STAR_solo(
                ['fastq1', 'fastq2'],
                'path/to/index',
                'output/path',
                self.umi_technology,
                whitelist_path='path/to/whitelist',
                n_threads=8,
                temp_dir=self.temp_dir,
                nasc=False,
                overrides=None,
            )

            self.assertEqual({
                'bam': 'output/path/Aligned.sortedByCoord.out.bam',
                'gene': {
                    'raw': {
                        'barcodes': 'output/path/Solo.out/Gene/raw/barcodes.tsv',
                        'features': 'output/path/Solo.out/Gene/raw/features.tsv',
                        'matrix': 'output/path/Solo.out/Gene/raw/matrix.mtx'
                    },
                    'filtered': {
                        'barcodes': 'output/path/Solo.out/Gene/filtered/barcodes.tsv',
                        'features': 'output/path/Solo.out/Gene/filtered/features.tsv',
                        'matrix': 'output/path/Solo.out/Gene/filtered/matrix.mtx'
                    }
                }
            }, result)
            arguments = utils.combine_arguments(
                utils.combine_arguments(config.STAR_ARGUMENTS, self.umi_technology.chemistry.to_starsolo_arguments()),
                config.STAR_SOLO_ARGUMENTS
            )
            arguments.update({
                '--soloStrand': 'Forward',
                '--soloCBwhitelist': 'path/to/whitelist',
                '--genomeDir': 'path/to/index',
                '--runThreadN': 8,
                '--outFileNamePrefix': 'output/path/',
                '--outTmpDir': mock.ANY,
                '--outBAMsortingBinsN': mock.ANY,
                '--readFilesIn': ['fastq1', 'fastq2']
            })
            arguments_to_list.assert_called_once_with(arguments)
            run_executable.assert_called_once()

    def test_STAR_solo_smartseq(self):
        with mock.patch('dynast.align.utils.get_STAR_binary_path', return_value='path/to/STAR'), \
            mock.patch('dynast.align.utils.arguments_to_list') as arguments_to_list, \
            mock.patch('dynast.align.utils.run_executable') as run_executable:
            # Write fake manifest
            manifest_path = os.path.join(self.temp_dir, 'manifest.csv')
            with open(manifest_path, 'w') as f:
                f.write('cell_0,fastq1,fastq2\ncell_1,fastq3,fastq4\n')

            result = align.STAR_solo(
                [manifest_path],
                'path/to/index',
                'output/path',
                self.smartseq_technology,
                whitelist_path=None,
                n_threads=8,
                temp_dir=self.temp_dir,
                nasc=False,
                overrides=None,
            )

            self.assertEqual({
                'bam': 'output/path/Aligned.sortedByCoord.out.bam',
                'gene': {
                    'raw': {
                        'barcodes': 'output/path/Solo.out/Gene/raw/barcodes.tsv',
                        'features': 'output/path/Solo.out/Gene/raw/features.tsv',
                        'matrix': 'output/path/Solo.out/Gene/raw/matrix.mtx'
                    },
                    'filtered': {
                        'barcodes': 'output/path/Solo.out/Gene/filtered/barcodes.tsv',
                        'features': 'output/path/Solo.out/Gene/filtered/features.tsv',
                        'matrix': 'output/path/Solo.out/Gene/filtered/matrix.mtx'
                    }
                }
            }, result)
            arguments = utils.combine_arguments(
                config.STAR_ARGUMENTS, self.smartseq_technology.chemistry.to_starsolo_arguments()
            )
            arguments = utils.combine_arguments(arguments, self.smartseq_technology.additional_args)
            arguments.update({
                '--soloStrand': 'Forward',
                '--genomeDir': 'path/to/index',
                '--runThreadN': 8,
                '--outFileNamePrefix': 'output/path/',
                '--outTmpDir': mock.ANY,
                '--outBAMsortingBinsN': mock.ANY,
                '--readFilesManifest': mock.ANY
            })
            arguments_to_list.assert_called_once_with(arguments)
            run_executable.assert_called_once()

    def test_STAR_solo_smartseq_nasc(self):
        with mock.patch('dynast.align.utils.get_STAR_binary_path', return_value='path/to/STAR'), \
            mock.patch('dynast.align.utils.arguments_to_list') as arguments_to_list, \
            mock.patch('dynast.align.utils.run_executable') as run_executable:
            # Write fake manifest
            manifest_path = os.path.join(self.temp_dir, 'manifest.csv')
            with open(manifest_path, 'w') as f:
                f.write('cell_0,fastq1,fastq2\ncell_1,fastq3,fastq4\n')

            result = align.STAR_solo(
                [manifest_path],
                'path/to/index',
                'output/path',
                self.smartseq_technology,
                whitelist_path=None,
                n_threads=8,
                temp_dir=self.temp_dir,
                nasc=True,
                overrides=None,
            )

            self.assertEqual({
                'bam': 'output/path/Aligned.sortedByCoord.out.bam',
                'gene': {
                    'raw': {
                        'barcodes': 'output/path/Solo.out/Gene/raw/barcodes.tsv',
                        'features': 'output/path/Solo.out/Gene/raw/features.tsv',
                        'matrix': 'output/path/Solo.out/Gene/raw/matrix.mtx'
                    },
                    'filtered': {
                        'barcodes': 'output/path/Solo.out/Gene/filtered/barcodes.tsv',
                        'features': 'output/path/Solo.out/Gene/filtered/features.tsv',
                        'matrix': 'output/path/Solo.out/Gene/filtered/matrix.mtx'
                    }
                }
            }, result)
            arguments = utils.combine_arguments(
                utils.combine_arguments(
                    config.STAR_ARGUMENTS, self.smartseq_technology.chemistry.to_starsolo_arguments()
                ), config.NASC_ARGUMENTS
            )
            arguments = utils.combine_arguments(arguments, self.smartseq_technology.additional_args)
            arguments.update({
                '--soloStrand': 'Forward',
                '--genomeDir': 'path/to/index',
                '--runThreadN': 8,
                '--outFileNamePrefix': 'output/path/',
                '--outTmpDir': mock.ANY,
                '--outBAMsortingBinsN': mock.ANY,
                '--readFilesManifest': mock.ANY
            })
            arguments_to_list.assert_called_once_with(arguments)
            run_executable.assert_called_once()

    def test_STAR_solo_overrides(self):
        with mock.patch('dynast.align.utils.get_STAR_binary_path', return_value='path/to/STAR'), \
            mock.patch('dynast.align.utils.arguments_to_list') as arguments_to_list, \
            mock.patch('dynast.align.utils.run_executable') as run_executable:
            overrides = {'--override1': 'override_value1', '--outFilterScoreMinOverLread': 1}
            result = align.STAR_solo(
                ['fastq1', 'fastq2'],
                'path/to/index',
                'output/path',
                self.umi_technology,
                whitelist_path=None,
                n_threads=8,
                temp_dir=self.temp_dir,
                nasc=False,
                overrides=overrides,
            )

            self.assertEqual({
                'bam': 'output/path/Aligned.sortedByCoord.out.bam',
                'gene': {
                    'raw': {
                        'barcodes': 'output/path/Solo.out/Gene/raw/barcodes.tsv',
                        'features': 'output/path/Solo.out/Gene/raw/features.tsv',
                        'matrix': 'output/path/Solo.out/Gene/raw/matrix.mtx'
                    },
                    'filtered': {
                        'barcodes': 'output/path/Solo.out/Gene/filtered/barcodes.tsv',
                        'features': 'output/path/Solo.out/Gene/filtered/features.tsv',
                        'matrix': 'output/path/Solo.out/Gene/filtered/matrix.mtx'
                    }
                }
            }, result)
            arguments = utils.combine_arguments(
                utils.combine_arguments(config.STAR_ARGUMENTS, self.umi_technology.chemistry.to_starsolo_arguments()),
                config.STAR_SOLO_ARGUMENTS
            )
            arguments.update({
                '--soloStrand': 'Forward',
                '--soloCBwhitelist': 'None',
                '--genomeDir': 'path/to/index',
                '--runThreadN': 8,
                '--outFileNamePrefix': 'output/path/',
                '--outTmpDir': mock.ANY,
                '--outBAMsortingBinsN': mock.ANY,
                '--readFilesIn': ['fastq1', 'fastq2']
            })
            arguments.update(overrides)
            arguments_to_list.assert_called_once_with(arguments)
            run_executable.assert_called_once()

    def test_align(self):
        pass
