import os
from unittest import TestCase, mock

import dynast.ref as ref
from tests import mixins


class TestRef(mixins.TestMixin, TestCase):

    def test_STAR_genomeGenerate(self):
        with mock.patch('dynast.ref.utils.get_STAR_binary_path', return_value='path/to/STAR'), \
            mock.patch('dynast.ref.utils.run_executable') as run_executable:
            # Write fake FASTA
            fasta_path = os.path.join(self.temp_dir, 'fasta.fa')
            with open(fasta_path, 'w') as f:
                f.write('>entry1\nACGT\n')

            result = ref.STAR_genomeGenerate(
                fasta_path, 'path/to/gtf', 'path/to/index', n_threads=8, memory=16 * 1024**3, temp_dir=self.temp_dir
            )
            self.assertEqual({'index': 'path/to/index'}, result)
            run_executable.assert_called_once_with([
                'path/to/STAR',
                '--runMode',
                'genomeGenerate',
                '--genomeDir',
                'path/to/index',
                '--genomeFastaFiles',
                fasta_path,
                '--sjdbGTFfile',
                'path/to/gtf',
                '--runThreadN',
                8,
                '--limitGenomeGenerateRAM',
                16 * 1024**3,
                '--genomeSAsparseD',
                mock.ANY,
                '--genomeSAindexNbases',
                mock.ANY,
                '--genomeChrBinNbits',
                mock.ANY,
                '--outTmpDir',
                mock.ANY,
            ])
