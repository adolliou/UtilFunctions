import astropy.io.fits as fits
import astropy.units as u
from astropy.time import Time
import astropy.constants
import numpy as np
from create_dict_file import *
from tqdm import tqdm

class UtilFunctions:

    @staticmethod
    def find_closest_time(list_to_find: list, list_ref: list, window_to_find=-1, window_ref=-1, time_delay=True,
                          maximal_threshold=15 * u.s):
        """
        Returns a list (of size length(list_to_find) with the index of the closest Fits files.
        Returns
        """
        list_index = np.arr
        for ii, fits_path_to_find in enumerate(tqdm(list_to_find, desc="find_closest_time")):
            with fits.open(fits_path_to_find) as hdul:
                hdu = hdul[window_to_find]
                date_obs_to_find = Time(hdu.header["DATE_OBS"])
                if time_delay:
                    dsun_obs_to_find = hdu.header["DSUN_OBS"]
                time_diff = []
                for jj, fits_path_ref in enumerate(list_ref):
                    with fits.open(fits_path_ref) as hdul_tmp:
                        hdu_tmp = hdul_tmp[window_ref]
                        date_obs_ref = Time(hdu_tmp.header["DATE_OBS"])

                        if time_delay:
                            dsun_obs_ref = hdu_tmp.header["DSUN_OBS"]
                            dt = (np.array(dsun_obs_to_find) - np.array(dsun_obs_ref)) / astropy.constants.c.value
                            date_obs_ref = date_obs_ref + dt*u.s
                        time_diff.append(np.abs((date_obs_to_find - date_obs_ref).to(u.s).value))
                        hdul_tmp.close()
                list_index.append(np.array(time_diff).argmin())
                hdul.close()
            list_index = np.array(list_index, dtype='int')
            selection_error = (list_index > maximal_threshold.to(u.s).value)
            if selection_error.sum() > 0:
                raise ValueError("Threshold delta time of %i s attained" % (maximal_threshold.to(u.s).value))








