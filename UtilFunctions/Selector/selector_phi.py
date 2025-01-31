import os.path
from .selector import Selector
from urllib.parse import urljoin
import re



class SelectorPhiMPS(Selector):

    default_base_url = "/data/solo/phi/data/fmdb"

    def __init__(self, release=6.0, level=2, base_url=None, release_dict=None, level_dict=None, 
                 year_suffix="", month_suffix="", day_suffix=""):
        """

        :param release: Release number (e.g. 6.0)
        :param level: level of the data (e.g. 2)
        :param base_url: url or path where the release folders are located. If none, then set to the default_base_url
        :param release_dict: dict containing the names of each release folders. Change if needed
        :param level_dict: dict containing the names of each level folders. Change if needed.
        :param year_suffix: Suffix in the year index of the path
        :param month_suffix: Suffix in the month index of the path
        :param day_suffix: Suffix in the day index of the path

        """
        if release_dict is None:
            self.release_dict = {
                "public": "public",
            }
        else:
            self.release_dict = release_dict
        if level_dict is None:
            self.level_dict = {
                "2": "l2",
            }
        else:
            self.level_dict = level_dict

        if base_url is None:
            base_url = SelectorPhiMPS.default_base_url
        url = base_url + '/' + self.release_dict[str(release)] + '/' + self.level_dict[str(level)]
        super().__init__(release_url_basis=url, year_suffix=year_suffix, month_suffix=month_suffix, day_suffix=day_suffix)

    
    def get_regex(self):
        self.re_filename = re.compile(
            r'''
        solo
        _(?P<level>L[123])
        _phi
        -(?P<instrument>(hrt|fdt))
        (?P<parameter>(vlos|stockes|icnt|FullModel|chi2|bmag|blos|binc|bazi|))?
        (?P<exposure>-short)?
        _(?P<time>\d{8}T\d{6})
        (?P<miliseconds>\d+)?
       _(?P<version>V\d{12})
       _(?P<reference_number>\d{10})
       .fits
       (?P<compression>(.gz))?
        ''',
        re.VERBOSE
        )

# "solo_L2_phi-hrt-vlos_20230410T214409_V202408301102_0344100916.fits.gz"