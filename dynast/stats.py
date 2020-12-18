import datetime as dt
import json
import sys
from contextlib import contextmanager

from . import __version__, config


class Step:

    def __init__(self, skipped=False, **kwargs):
        self.start_time = None
        self.end_time = None
        self.elapsed = None
        self.skipped = skipped
        self.extra = kwargs

    def start(self):
        self.start_time = dt.datetime.now()

    def end(self):
        self.end_time = dt.datetime.now()
        self.elapsed = (self.end_time - self.start_time).total_seconds()

    def to_dict(self):
        return {
            'start_time': None if self.skipped else self.start_time.isoformat(),
            'end_time': None if self.skipped else self.end_time.isoformat(),
            'elapsed': None if self.skipped else self.elapsed,
            'skipped': self.skipped,
            **self.extra
        }


class Stats:
    """Class used to collect run statistics.
    """

    def __init__(self):
        self.call = None

        self.start_time = None
        self.end_time = None
        self.elapsed = None

        self.steps = {}
        self.step_order = []

        self.STAR = {
            'version': config.get_STAR_version(),
            'command': None,
        }

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

    @contextmanager
    def step(self, key, skipped=False, **kwargs):
        step = Step(skipped=skipped, **kwargs)
        self.steps[key] = step
        self.step_order.append(key)
        if not skipped:
            step.start()
        yield
        if not skipped:
            step.end()

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
            'step_order': self.step_order,
            'steps': {key: step.to_dict()
                      for key, step in self.steps.items()}
        }


STATS = Stats()
