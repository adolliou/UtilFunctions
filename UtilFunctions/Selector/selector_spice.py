from selector import Selector
import re

class SelectorSpice(Selector):

    def __init__(self):
        super().__init__()

    def get_regex(self):
        self.re_filename = re.compile(
            r'''
            solo
            _(?P<level>L[123])
            _spice
                (?P<concat>-concat)?
                -(?P<slit>[wn])
                -(?P<type>(ras|sit|exp))
                (?P<db>-db)?
                (?P<int>-int)?
            _(?P<time>\d{8}T\d{6})
            _(?P<version>V\d{2})
            _(?P<SPIOBSID>\d+)-(?P<RASTERNO>\d+)
            (?P<ext>\..*)
            ''',
            re.VERBOSE
            )