import array
from unittest import TestCase

import pysam

import dynast.preprocessing.consensus as consensus

from .. import mixins


class TestConsensus(mixins.TestMixin, TestCase):

    def test_call_consensus_from_reads(self):
        header = pysam.AlignmentHeader.from_dict({
            'HD': {
                'VN': pysam.version.__samtools_version__,
                'SO': 'unsorted'
            },
        })
        # Desired final consensus sequence is: ACTGCCAGCTGC
        # R1 and R2 are simple overlaps
        # R3 should override R2 due to base quality, and one base has a tie for which
        # the reference should be taken
        # R4 has an insertion that should be ignored (not include in the consensus)
        # R5 has a deletion that should be skipped in the consensus
        # R6 has a low quality base, so the reference base should be taken
        read1 = pysam.AlignedSegment(header)
        read1.query_name = 'read1'
        read1.query_sequence = 'ACTG'
        read1.query_qualities = array.array('B', [30] * 4)
        read1.reference_start = 1
        read1.cigarstring = '4M'
        read1.set_tag('MD', '4')

        read2 = pysam.AlignedSegment(header)
        read2.query_name = 'read2'
        read2.query_sequence = 'TGCT'
        read2.query_qualities = array.array('B', [30] * 4)
        read2.reference_start = 3
        read2.cigarstring = '4M'
        read2.set_tag('MD', '4')

        read3 = pysam.AlignedSegment(header)
        read3.query_name = 'read3'
        read3.query_sequence = 'ACAG'
        read3.query_qualities = array.array('B', [30, 40, 30, 30])
        read3.reference_start = 5
        read3.cigarstring = '4M'
        read3.set_tag('MD', '0C0T2')

        read4 = pysam.AlignedSegment(header)
        read4.query_name = 'read4'
        read4.query_sequence = 'AGCTA'
        read4.query_qualities = array.array('B', [30] * 5)
        read4.reference_start = 7
        read4.cigarstring = '4M1I'
        read4.set_tag('MD', '4')

        read5 = pysam.AlignedSegment(header)
        read5.query_name = 'read5'
        read5.query_sequence = 'CTG'
        read5.query_qualities = array.array('B', [30] * 3)
        read5.reference_start = 9
        read5.cigarstring = '2M1D1M'
        read5.set_tag('MD', '2^A1')

        read6 = pysam.AlignedSegment(header)
        read6.query_name = 'read6'
        read6.query_sequence = 'GG'
        read6.query_qualities = array.array('B', [20] * 2)
        read6.reference_start = 12
        read6.cigarstring = '2M'
        read6.set_tag('MD', '1C')

        result = consensus.call_consensus_from_reads([read6, read1, read3, read2, read5, read4], header, quality=27)
        self.assertEqual('ACTGCCAGCTGC', result.query_sequence)
        self.assertEqual(1, result.reference_start)
        self.assertEqual(14, result.reference_end)
        self.assertEqual('10M1D2M', result.cigarstring)
        self.assertEqual('5T4^A2', result.get_tag('MD'))
        print(result.get_aligned_pairs(matches_only=True, with_seq=True))
