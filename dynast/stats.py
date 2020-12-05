import datetime as dt
import json
import sys

from . import __version__, config


class Stats:
    """Class used to collect run statistics.
    """

    def __init__(self):
        self.call = None

        self.start_time = None
        self.end_time = None

        self.STAR = {
            'version': config.get_STAR_version(),
            'command': None,
        }

        self.elapsed = None
        self.version = __version__

    def start(self):
        """Start collecting statistics.

        Sets start time, the command line call.
        """
        self.start_time = dt.datetime.now()
        self.call = ' '.join(sys.argv)

    def end(self):
        """End collecting statistics.
        """
        self.end_time = dt.datetime.now()
        self.elapsed = (self.end_time - self.start_time).total_seconds()

    def save(self, path):
        """Save statistics as JSON to path.

        :param path: path to JSON
        :type path: str

        :return: path to saved JSON
        :rtype: str
        """
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)
        return path

    def to_dict(self):
        """Convert statistics to dictionary, so that it is easily parsed
        by the report-rendering functions.
        """
        return {
            'version': self.version,
            'call': self.call,
            'start_time': self.start_time.isoformat(),
            'end_time': self.end_time.isoformat(),
            'elapsed': self.elapsed,
            'STAR': self.STAR,
        }


STATS = Stats()
