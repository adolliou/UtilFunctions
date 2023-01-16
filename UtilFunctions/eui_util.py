from . import common_util
import numpy as np
from astropy.wcs import WCS


class EUIUtil:
    @staticmethod
    def extract_EUI_coordinates(hdr):
        w = WCS(hdr)
        idx_lon = np.where(np.array(w.wcs.ctype, dtype="str") == "HPLN-TAN")[0][0]
        idx_lat = np.where(np.array(w.wcs.ctype, dtype="str") == "HPLT-TAN")[0][0]
        x, y = np.meshgrid(np.arange(w.pixel_shape[idx_lon]),
                           np.arange(w.pixel_shape[idx_lat]), )  # t dépend de x,
        # should reproject on a new coordinate grid first : suppose slits at the same time :
        longitude, latitude = w.pixel_to_world(x, y)
        dsun_obs_large = w.to_header()["DSUN_OBS"]
        return common_util.UtilFunctions.ang2pipi(longitude),\
            common_util.UtilFunctions.ang2pipi(latitude), dsun_obs_large
