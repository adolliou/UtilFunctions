import os.path
from .selector import Selector
from urllib.parse import urljoin
import re
from astropy.time import Time



class SelectorAspiicsMPS(Selector):

    default_base_url = "/scratch/proba3/"

    def __init__(self, release="public", level=2, base_url=None, release_dict=None, level_dict=None, 
                 year_suffix="", month_suffix="", day_suffix="", 
                 verbose=1):
        """

        :param release: Release type (e.g. public)
        :param level: level of the data (e.g. 2)
        :param base_url: url or path where the release folders are located. If none, then set to the default_base_url
        :param release_dict: dict containing the names of each release folders. Change if needed
        :param level_dict: dict containing the names of each level folders. Change if needed.
        :param year_suffix: Suffix in the year index of the path
        :param month_suffix: Suffix in the month index of the path
        :param day_suffix: Suffix in the day index of the path
        :param verbose: Level of printing on the terminal. 

        """
        if release_dict is None:
            self.release_dict = {
                "public": "public",
            }
        else:
            self.release_dict = release_dict
        if level_dict is None:
            self.level_dict = {
                2: "l2",
            }
        else:
            self.level_dict = level_dict

        if base_url is None:
            base_url = SelectorAspiicsMPS.default_base_url
        url = base_url + '/' + self.release_dict[str(release)] + '/' + self.level_dict[level]
        super().__init__(release_url_basis=url, year_suffix=year_suffix, month_suffix=month_suffix, day_suffix=day_suffix, verbose=verbose,)

    
    
    def _find_url_from_file(self, fits_file_name):
        time = self._find_time_from_file(fits_file_name)
        return (self.release_url_basis + '/' + f'{time.ymdhms[0]:04d}' + self.year_suffix + '-' + f"{time.ymdhms[1]:02d}" + self.month_suffix + '-' +
                f"{time.ymdhms[2]:02d}") + self.day_suffix
    
    def _find_url_from_time(self, time: Time):
        url = self.release_url_basis + '/' + f"{time.ymdhms[0]:04d}" + self.year_suffix + '-' + f"{time.ymdhms[1]:02d}" + self.month_suffix \
        + '-' + f"{time.ymdhms[2]:02d}" + self.day_suffix
        return url

    def get_regex(self):
        self.re_filename = re.compile(
            r'''
        aspiics
        _(?P<parameter>(wb|p3|p2|p1|he|fe))?
        _(?P<level>l[123])
        _(?P<version>V\w{14})
        _(?P<time>\d{8}T\d{6})
        (?P<miliseconds>\d+)?
       .fits
       (?P<compression>(.gz))?
        ''',
        re.VERBOSE
        )

# "aspiics_fe_l2_1700404A000121_20251002T125214.fits"