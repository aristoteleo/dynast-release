import logging
import math
import os
import tempfile

from .config import get_STAR_binary_path
from .utils import open_as_text, run_executable

logger = logging.getLogger(__name__)


def STAR_genomeGenerate(
    fasta_path,
    gtf_path,
    index_dir,
    n_threads=8,
    memory=16 * 1024**3,
    temp_dir=None,
):
    """Generate a STAR index from a reference.

    :param fasta_path:
    :type fasta_path:
    :param gtf_path:
    :type gtf_path:
    :param index_dir:
    :type index_dir:
    :param n_threads:
    :type n_threads:
    :param memory:
    :type memory:
    :param temp_dir:
    :type temp_dir:

    :return:
    :rtype:
    """
    print(f'Calculating optimal STAR index parameters for {memory} bytes of memory')
    # Calculate genome length and number of FASTA entries
    genome_length = 0
    n_entries = 0
    with open_as_text(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                n_entries += 1
            else:
                genome_length += len(line.strip())

    # Taken from STAR manual (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
    # genomeSAindexNbases = min(14, log2(GenomeLength)/2 - 1)
    sa_index_n_bases = min(14, int(math.log2(genome_length) / 2 - 1))
    # genomeChrBinNbits = min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)])
    chr_bin_n_bits = min(18, int(math.log2(genome_length / n_entries)))

    # Adapted from Cellranger
    sa_index_memory = (4**sa_index_n_bases) * 8
    genome_memory = genome_length
    sa_sparse_d = max(
        1, math.ceil((8 * genome_length) / (max(1, memory - 2 * 1024**3) - genome_memory - sa_index_memory))
    )
    min_memory = genome_memory + sa_index_memory + 3 * 1024**3
    if memory < min_memory:
        raise Exception(f'STAR requires at least {min_memory} to index this reference.')

    print(
        'Calculated parameters: '
        f'genome_length={genome_length} '
        f'n_entries={n_entries} '
        f'sa_index_n_bases={sa_index_n_bases} '
        f'chr_bin_n_bits={chr_bin_n_bits} '
        f'sa_sparse_d={sa_sparse_d}'
    )

    print(f'Generating STAR index to {index_dir}')
    command = [get_STAR_binary_path(), '--runMode', 'genomeGenerate']
    command += ['--genomeDir', index_dir]
    command += ['--genomeFastaFiles', fasta_path]
    command += ['--sjdbGTFfile', gtf_path]
    command += ['--runThreadN', n_threads]
    command += ['--limitGenomeGenerateRAM', memory]
    command += ['--genomeSAsparseD', sa_sparse_d]
    command += ['--genomeSAindexNbases', sa_index_n_bases]
    command += ['--genomeChrBinNbits', chr_bin_n_bits]
    command += [
        '--outTmpDir',
        os.path.join(temp_dir, f'{tempfile.gettempprefix()}{next(tempfile._get_candidate_names())}')
    ]

    run_executable(command)

    return {'index': index_dir}


def ref():
    pass
