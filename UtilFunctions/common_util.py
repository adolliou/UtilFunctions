import astropy.io.fits as fits
import astropy.units as u
from astropy.time import Time
import astropy.constants
import numpy as np
from tqdm import tqdm
from astropy.wcs import WCS
import cv2
import scipy


class CommonUtil:

    @staticmethod
    def find_closest_dict_index(utc_eui, dict_file_reference, threshold_time):

        delta_time = np.array([np.abs((utc_eui - n).to(u.s).value) for n in dict_file_reference["date-avg"]])

        closest_index = delta_time.argmin()
        delta_time_min = delta_time[closest_index]
        if delta_time_min > threshold_time:
            raise ValueError("Delta time between EUI and SPICE file "
                             "equal to %2f s > %.2f" % (delta_time_min, threshold_time))
        return closest_index, delta_time_min
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
                date_obs_to_find = Time(hdu.header["DATE-AVG"])
                if time_delay:
                    dsun_obs_to_find = hdu.header["DSUN_OBS"]
                time_diff = []
                for jj, fits_path_ref in enumerate(list_ref):
                    with fits.open(fits_path_ref) as hdul_tmp:
                        hdu_tmp = hdul_tmp[window_ref]
                        date_obs_ref = Time(hdu_tmp.header["DATE-AVG"])

                        if time_delay:
                            dsun_obs_ref = hdu_tmp.header["DSUN_OBS"]
                            dt = (np.array(dsun_obs_to_find) - np.array(dsun_obs_ref)) / astropy.constants.c.value
                            date_obs_ref = date_obs_ref + dt * u.s
                        time_diff.append(np.abs((date_obs_to_find - date_obs_ref).to(u.s).value))
                        hdul_tmp.close()
                list_index.append(np.array(time_diff).argmin())
                hdul.close()
            list_index = np.array(list_index, dtype='int')
            selection_error = (list_index > maximal_threshold.to(u.s).value)
            if selection_error.sum() > 0:
                raise ValueError("Threshold delta time of %i s attained" % (maximal_threshold.to(u.s).value))

    @staticmethod
    def ang2pipi(ang):
        """ put angle between ]-180, +180] deg """
        pi = u.Quantity(180, 'deg')
        return - ((- ang + pi) % (2 * pi) - pi)

    @staticmethod
    def interpol2d(image, x, y, order=1, fill=0, opencv=False, dst=None):
        """
        Interpolates in 2D image using either map_coordinates or opencv

        data: image to interpolate
        x, y: coordinates (in pixels) at which to interpolate the image
        order: if opencv is True:  0=nearest neighbor, 1=linear, 2=cubic
               if opencv is False: the order of the spline interpolation used by
                                   map_coordinates (see scipy documentation)
        opencv: If True, uses opencv
                If False, uses scipy.ndimage.map_coordinates
                opencv can use only 32 bits floating point coordinates input
        fill: constant value usesd to fill in the edges
        dst: if present, ndarray in which to place the result
        """

        bad = np.logical_or(x == np.nan, y == np.nan)
        x = np.where(bad, -1, x)
        y = np.where(bad, -1, y)

        if dst is None:
            dst = np.empty(x.shape, dtype=image.dtype)

        if opencv:
            if order == 0:
                inter = cv2.INTER_NEAREST
            elif order == 1:
                inter = cv2.INTER_LINEAR
            elif order == 2:
                inter = cv2.INTER_CUBIC
            cv2.remap(image,
                      x.astype(np.float32),  # converts to float 32 for opencv
                      y.astype(np.float32),  # does nothing with default dtype
                      inter,  # interpolation method
                      dst,  # destination array
                      cv2.BORDER_CONSTANT,  # fills in with constant value
                      fill)  # constant value
        else:
            coords = np.stack((y.ravel(), x.ravel()), axis=0)

            scipy.ndimage.map_coordinates(image,  # input array
                                          coords,  # array of coordinates
                                          order=order,  # spline order
                                          mode='constant',  # fills in with constant value
                                          cval=fill,  # constant value
                                          output=dst.ravel(),
                                          prefilter=False)

        return dst
