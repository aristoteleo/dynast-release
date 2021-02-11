import bisect
import re
from itertools import product

from .. import utils
from ..logging import logger


class GTF:
    """Utility class to easily read and parse GTF files.

    :param gtf_path: path to GTF file
    :type gtf_path: str
    """
    PARSER = re.compile(
        r'''
        ^(?P<seqname>.+?)\s+    # chromosome
        .*?\t                   # source
        (?P<feature>.+?)\s+     # feature: transcript, exon, etc.
        (?P<start>[0-9]+?)\s+   # start position (1-indexed)
        (?P<end>[0-9]+?)\s+     # end position (1-indexed, inclusive)
        .*?\s+                  # score
        (?P<strand>\+|-|\.)\s+  # +, -, . indicating strand
        .*?\s+                  # frame
        (?P<group>.*)           # groups
    ''', re.VERBOSE
    )
    GROUP_PARSER = re.compile(r'(?P<key>\S+?)\s*"(?P<value>.+?)"')

    def __init__(self, gtf_path):
        self.gtf_path = gtf_path

    @staticmethod
    def parse_entry(line):
        """Parse a single GTF entry.

        :param line: a line in the GTF file
        :type line: str

        :return: parsed GTF information
        :rtype: dict
        """
        match = GTF.PARSER.match(line)
        if match:
            groupdict = match.groupdict()
            groupdict['start'] = int(groupdict['start'])
            groupdict['end'] = int(groupdict['end'])
            groupdict['group'] = dict(GTF.GROUP_PARSER.findall(groupdict.get('group', '')))
            if not groupdict['group']:
                logger.warning(f'Failed to parse GTF attributes of entry: {line}')

            return groupdict
        logger.warning(f'Failed to parse GTF entry: {line}')
        return None

    def entries(self):
        """Generator that yields one GTF entry at a time.

        :return: a generator that yields a dict of the GTF entry
        :rtype: generator
        """
        with utils.open_as_text(self.gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#') or line.isspace():
                    continue

                yield GTF.parse_entry(line)


class Segment:
    """Class to represent an integer interval segment, zero-indexed.

    :param start: start position
    :type start: int
    :param end: end position
    :type end: int
    """

    def __init__(self, start, end):
        assert end > start
        self._start = start
        self._end = end

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def width(self):
        return self.end - self.start

    def is_in(self, i):
        return i >= self.start and i < self.end

    def is_exclusive(self, segment):
        return self.end <= segment.start or self.start >= segment.end

    def is_overlapping(self, segment):
        return not self.is_exclusive(segment)

    def is_subset(self, segment):
        return self.start >= segment.start and self.end <= segment.end

    def is_superset(self, segment):
        return self.start <= segment.start and self.end >= segment.end

    def __iter__(self):
        return iter((self.start, self.end))

    def __eq__(self, other):
        return (self.start, self.end) == (other.start, other.end)

    def __lt__(self, other):
        return (self.start, self.end) < (other.start, other.end)

    def __gt__(self, other):
        return (self.start, self.end) > (other.start, other.end)

    def __str__(self):
        return str((self.start, self.end))

    def __repr__(self):
        return str(self)


class SegmentCollection:
    """Class to represent a collection of integer interval segments, zero-indexed.

    :param segments: list of initial segments, defaults to `None`
    :type segments: list, optional
    """

    def __init__(self, segments=None):
        self.segments = sorted(segments) if segments else []
        if self.segments:
            self.collapse()

    @property
    def start(self):
        return self.segments[0].start if self.segments else -1

    @property
    def end(self):
        return self.segments[-1].end if self.segments else -1

    def add_segment(self, segment):
        bisect.insort_left(self.segments, segment)
        self.collapse()

    def add_collection(self, collection):
        self.segments = sorted(self.segments + collection.segments)
        self.collapse()

    def __iter__(self):
        return iter(self.segments)

    def __getitem__(self, i):
        return self.segments[i]

    def __len__(self):
        return len(self.segments)

    def collapse(self):
        segments = []
        combined = None
        for segment in self.segments:
            if combined is None:
                combined = segment
                continue
            if combined.is_overlapping(segment):
                combined = Segment(min(combined.start, segment.start), max(combined.end, segment.end))
            else:
                segments.append(combined)
                combined = segment
        if combined is not None:
            segments.append(combined)

        self.segments = segments

    def span_is_exclusive(self, collection):
        return self.end <= collection.start or self.start >= collection.end

    def is_overlapping(self, collection):
        if self.span_is_exclusive(collection):
            return False
        return any(segment_1.is_overlapping(segment_2) for segment_1, segment_2 in product(self, collection))

    def is_subset(self, collection):
        if self.span_is_exclusive(collection):
            return False
        # Assume segments are sorted
        collection_i = 0
        for cls_segment in self.segments:
            found = False
            for collection_segment in collection.segments[collection_i:]:
                if cls_segment.is_subset(collection_segment):
                    found = True
                    break
                collection_i += 1

            if not found:
                return False
        return True

    def is_superset(self, collection):
        return collection.is_subset(self)

    @classmethod
    def from_positions(cls, positions):
        collection = cls()
        positions = sorted(positions)
        start = positions[0]
        end = positions[0] + 1
        for i in range(1, len(positions)):
            left = positions[i - 1]
            right = positions[i]

            if right - left > 1:
                collection.add_segment(Segment(start, end))
                start = right
            end = right + 1

        collection.add_segment(Segment(start, end))
        return collection

    @classmethod
    def from_collections(cls, *collections):
        segments = []
        for collection in collections:
            segments.extend(collection.segments)
        return cls(segments=segments)

    def __str__(self):
        return f'SegmentCollection {str(self.segments)}'

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if len(self) != len(other):
            return False
        for segment1, segment2 in zip(self.segments, other.segments):
            if segment1 != segment2:
                return False
        return True


def parse_gtf(gtf_path):
    """Parse GTF for gene and transcript informations.

    :param gtf_path: path to GTF
    :type gtf_path: str

    :return: (gene information, transcript information)
    :rtype: (dictionary, dictionary)
    """
    gtf = GTF(gtf_path)

    gene_infos = {}
    transcript_exons = {}
    transcript_infos = {}
    count = 0
    for gtf_entry in gtf.entries():
        count += 1

        start = gtf_entry['start'] - 1
        end = gtf_entry['end']
        strand = gtf_entry['strand']
        chrom = gtf_entry['seqname']
        if gtf_entry['feature'] == 'exon':
            transcript_id = gtf_entry['group']['transcript_id']
            transcript_exons.setdefault(transcript_id, SegmentCollection()).add_segment(Segment(start, end))
        elif gtf_entry['feature'] == 'transcript':
            transcript_id = gtf_entry['group']['transcript_id']
            gene_id = gtf_entry['group']['gene_id']
            gene_name = gtf_entry['group'].get('gene_name')
            gene_infos.setdefault(gene_id, {}).setdefault('transcripts', set()).add(transcript_id)
            segment = gene_infos[gene_id].setdefault('segment', Segment(start, end))
            gene_infos[gene_id].update({
                'segment': Segment(min(segment.start, start), max(segment.end, end)),
                'chr': chrom,
                'strand': strand,
                'gene_name': gene_name,
            })

            if transcript_id not in transcript_infos:
                attributes = {
                    'gene_id': gene_id,
                    'chr': chrom,
                    'segment': Segment(start, end),
                    'strand': strand,
                }
                transcript_infos[transcript_id] = attributes

    # Calculate introns
    for transcript_id in transcript_infos:
        attributes = transcript_infos[transcript_id]
        transcript_interval = attributes['segment']

        introns = SegmentCollection()
        exons = transcript_exons.get(transcript_id, SegmentCollection())
        if exons:
            if exons[0].start > transcript_interval.start:
                introns.add_segment(Segment(transcript_interval.start, exons[0].start))

            for i in range(len(exons) - 1):
                introns.add_segment(Segment(exons[i].end, exons[i + 1].start))

            if exons[-1].end < transcript_interval.end:
                introns.add_segment(Segment(exons[-1].end, transcript_interval.end))
        else:
            introns.add_segment(transcript_interval)

        attributes['exons'] = exons
        attributes['introns'] = introns

    # Cast all sets of transcript IDs to lists (this makes it faster to iterate through later)
    for info in gene_infos.values():
        info['transcripts'] = list(info.get('transcripts', set()))

    logger.info(f'Parsed {len(gene_infos)} genes and {len(transcript_infos)} ' f'transcripts from {count} GTF entries')
    return gene_infos, transcript_infos
