from unittest import TestCase

import dynast.preprocessing.gtf as gtf

from .. import mixins


class TestGtf(mixins.TestMixin, TestCase):

    def test_parser(self):
        line = '2\thavana\tgene\t2\t3\t.\t+\t.\tgene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC";'  # noqa
        match = gtf.GTF.PARSER.match(line)
        self.assertIsNotNone(match)
        self.assertEqual(
            {
                'seqname':
                    '2',
                'feature':
                    'gene',
                'start':
                    '2',
                'end':
                    '3',
                'strand':
                    '+',
                'group':
                    'gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC";'  # noqa
            },
            match.groupdict()
        )

    def test_group_parser(self):
        group = 'gene_id "ENSMUSG00000102693";'
        match = gtf.GTF.GROUP_PARSER.match(group)
        self.assertIsNotNone(match)
        self.assertEqual({'key': 'gene_id', 'value': 'ENSMUSG00000102693'}, match.groupdict())

    def test_parse_entry(self):
        line = '2\thavana\tgene\t2\t3\t.\t+\t.\tgene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC";'  # noqa
        self.assertEqual({
            'seqname': '2',
            'feature': 'gene',
            'start': 2,
            'end': 3,
            'strand': '+',
            'group': {
                'gene_id': 'ENSMUSG00000102693',
                'gene_version': '1',
                'gene_name': '4933401J01Rik',
                'gene_source': 'havana',
                'gene_biotype': 'TEC'
            }
        }, gtf.GTF.parse_entry(line))

    def test_entries(self):
        g = gtf.GTF(self.gtf_path)
        self.assertEqual(8, len(list(g.entries())))

    def test_parse_gtf(self):
        gene_infos, transcript_infos = gtf.parse_gtf(self.gtf_path)
        self.assertEqual({
            'GENE_C': {
                'transcripts': ['TRANSCRIPT_C'],
                'segment': gtf.Segment(1, 14),
                'chr': '2',
                'strand': '-',
                'gene_name': 'GENE_C_NAME'
            },
            'GENE_B': {
                'transcripts': ['TRANSCRIPT_B'],
                'segment': gtf.Segment(0, 10),
                'chr': '2',
                'strand': '+',
                'gene_name': 'GENE_B_NAME'
            },
            'GENE_A': {
                'transcripts': ['TRANSCRIPT_A'],
                'segment': gtf.Segment(0, 2),
                'chr': '1',
                'strand': '+',
                'gene_name': 'GENE_A_NAME'
            }
        }, gene_infos)
        self.assertEqual({
            'TRANSCRIPT_C': {
                'gene_id': 'GENE_C',
                'chr': '2',
                'segment': gtf.Segment(1, 14),
                'strand': '-',
                'exons': gtf.SegmentCollection([gtf.Segment(1, 3), gtf.Segment(9, 14)]),
                'introns': gtf.SegmentCollection([gtf.Segment(3, 9)])
            },
            'TRANSCRIPT_B': {
                'gene_id': 'GENE_B',
                'chr': '2',
                'segment': gtf.Segment(0, 10),
                'strand': '+',
                'exons': gtf.SegmentCollection([gtf.Segment(1, 2), gtf.Segment(4, 8)]),
                'introns': gtf.SegmentCollection([gtf.Segment(0, 1),
                                                  gtf.Segment(2, 4),
                                                  gtf.Segment(8, 10)])
            },
            'TRANSCRIPT_A': {
                'gene_id': 'GENE_A',
                'chr': '1',
                'segment': gtf.Segment(0, 2),
                'strand': '+',
                'exons': gtf.SegmentCollection([gtf.Segment(0, 2)]),
                'introns': gtf.SegmentCollection(),
            }
        }, transcript_infos)


class TestSegment(mixins.TestMixin, TestCase):

    def test_init(self):
        segment = gtf.Segment(3, 4)
        self.assertEqual(3, segment.start)
        self.assertEqual(4, segment.end)
        self.assertEqual(1, segment.width)
        with self.assertRaises(AssertionError):
            gtf.Segment(1, 0)

    def test_is_in(self):
        segment = gtf.Segment(0, 10)
        self.assertTrue(segment.is_in(0))
        self.assertFalse(segment.is_in(10))

    def test_is_exclusive(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(10, 20)
        self.assertTrue(segment1.is_exclusive(segment3))
        self.assertFalse(segment1.is_exclusive(segment2))

    def test_is_overlapping(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(10, 20)
        self.assertFalse(segment1.is_overlapping(segment3))
        self.assertTrue(segment1.is_overlapping(segment2))

    def test_is_subset(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(0, 5)
        segment4 = gtf.Segment(5, 15)
        self.assertTrue(segment2.is_subset(segment1))
        self.assertTrue(segment3.is_subset(segment1))
        self.assertFalse(segment4.is_subset(segment1))

    def test_is_superset(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(0, 5)
        segment4 = gtf.Segment(5, 15)
        self.assertTrue(segment1.is_superset(segment2))
        self.assertTrue(segment1.is_superset(segment3))
        self.assertFalse(segment1.is_superset(segment4))

    def test_comparison(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(0, 10)
        segment3 = gtf.Segment(5, 10)
        segment4 = gtf.Segment(0, 5)
        self.assertTrue(segment1 == segment2)
        self.assertTrue(segment1 < segment3)
        self.assertTrue(segment1 > segment4)


class TestSegmentCollection(mixins.TestMixin, TestCase):

    def test_init(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(5, 15)
        collection = gtf.SegmentCollection(segments=[segment1, segment2])
        self.assertEqual(0, collection.start)
        self.assertEqual(15, collection.end)
        self.assertEqual(1, len(collection))
        self.assertEqual(gtf.Segment(0, 15), collection.segments[0])

    def test_add_segment(self):
        segment1 = gtf.Segment(0, 5)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(10, 15)
        collection = gtf.SegmentCollection(segments=[segment1, segment3])
        collection.add_segment(segment2)
        self.assertEqual(3, len(collection))
        self.assertEqual([segment1, segment2, segment3], collection.segments)

    def test_add_collection(self):
        segment1 = gtf.Segment(0, 5)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(10, 15)
        collection1 = gtf.SegmentCollection(segments=[segment1, segment3])
        collection2 = gtf.SegmentCollection(segments=[segment2])
        collection1.add_collection(collection2)
        self.assertEqual(3, len(collection1))
        self.assertEqual([segment1, segment2, segment3], collection1.segments)

    def test_collapse(self):
        collection = gtf.SegmentCollection()
        collection.segments = [gtf.Segment(0, 10), gtf.Segment(5, 15)]
        collection.collapse()
        self.assertEqual(1, len(collection))
        self.assertEqual(gtf.Segment(0, 15), collection.segments[0])

    def test_span_is_exclusive(self):
        collection1 = gtf.SegmentCollection(segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)])
        collection2 = gtf.SegmentCollection(segments=[gtf.Segment(5, 10)])
        collection3 = gtf.SegmentCollection(segments=[gtf.Segment(15, 20)])
        self.assertFalse(collection1.span_is_exclusive(collection2))
        self.assertTrue(collection1.span_is_exclusive(collection3))

    def test_is_overlapping(self):
        collection1 = gtf.SegmentCollection(segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)])
        collection2 = gtf.SegmentCollection(segments=[gtf.Segment(5, 10)])
        collection3 = gtf.SegmentCollection(segments=[gtf.Segment(0, 10)])
        self.assertFalse(collection1.is_overlapping(collection2))
        self.assertTrue(collection1.is_overlapping(collection3))

    def test_is_subset(self):
        collection1 = gtf.SegmentCollection(segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)])
        collection2 = gtf.SegmentCollection(segments=[gtf.Segment(3, 13)])
        collection3 = gtf.SegmentCollection(segments=[gtf.Segment(1, 2)])
        self.assertFalse(collection2.is_subset(collection1))
        self.assertTrue(collection3.is_subset(collection1))

    def test_is_superset(self):
        collection1 = gtf.SegmentCollection(segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)])
        collection2 = gtf.SegmentCollection(segments=[gtf.Segment(3, 13)])
        collection3 = gtf.SegmentCollection(segments=[gtf.Segment(1, 2)])
        self.assertFalse(collection1.is_superset(collection2))
        self.assertTrue(collection1.is_superset(collection3))

    def test_from_positions(self):
        positions = [0, 1, 2, 3, 3, 5, 6, 7, 8]
        collection = gtf.SegmentCollection.from_positions(positions)
        self.assertEqual([gtf.Segment(0, 4), gtf.Segment(5, 9)], collection.segments)

    def test_from_collections(self):
        segment1 = gtf.Segment(0, 5)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(9, 15)
        collection1 = gtf.SegmentCollection(segments=[segment1, segment3])
        collection2 = gtf.SegmentCollection(segments=[segment2])
        collection = gtf.SegmentCollection.from_collections(collection1, collection2)
        self.assertEqual([gtf.Segment(0, 5), gtf.Segment(5, 15)], collection.segments)

    def test_comparison(self):
        collection1 = gtf.SegmentCollection(segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)])
        collection2 = gtf.SegmentCollection(segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)])
        collection3 = gtf.SegmentCollection(segments=[gtf.Segment(1, 2)])
        self.assertTrue(collection1 == collection2)
        self.assertFalse(collection1 == collection3)
