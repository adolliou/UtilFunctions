from . import common_util
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord

class EUIUtil:
    @staticmethod
    def extract_EUI_coordinates(hdr, dsun=True):
        w = WCS(hdr)
        idx_lon = np.where(np.array(w.wcs.ctype, dtype="str") == "HPLN-TAN")[0][0]
        idx_lat = np.where(np.array(w.wcs.ctype, dtype="str") == "HPLT-TAN")[0][0]
        x, y = np.meshgrid(np.arange(w.pixel_shape[idx_lon]),
                           np.arange(w.pixel_shape[idx_lat]), )  # t dépend de x,
        # should reproject on a new coordinate grid first : suppose slits at the same time :
        longitude, latitude = w.pixel_to_world(x, y)
        if dsun:
            dsun_obs_large = hdr["DSUN_OBS"]
            return common_util.CommonUtil.ang2pipi(longitude), \
                common_util.CommonUtil.ang2pipi(latitude), dsun_obs_large
        else:
            return common_util.CommonUtil.ang2pipi(longitude), common_util.CommonUtil.ang2pipi(latitude)

    @staticmethod
    def diff_rot(lat, wvl='default'):
        ''' Return the angular velocity difference between differential and
        Carrington rotation.
        Parameters
        ==========
        lat : float
            The latitude, in radians
        wvl : str (default: 'default'
            The wavelength, or the band to return the rotation from.
        Returns
        =======
        corr : float
            The difference in angular velocities between the differential and
            Carrington rotations, in radians per second:
                Δω(θ) = ω_Car - ω_diff(θ)
                with ω_Car = 360° / (25.38 days)
                and  ω_diff(θ) = A + B sin² θ + C sin⁴ θ
        '''
        p = {
            # ° day⁻¹; Hortin (2003):
            'EIT 171': (14.56, -2.65, 0.96),
            'EIT 195': (14.50, -2.14, 0.66),
            'EIT 284': (14.60, -0.71, -1.18),
            'EIT 304': (14.51, -3.12, 0.34),
        }
        p['default'] = p['EIT 195']
        A, B, C = p[wvl]
        A_car = 360 / 25.38  # ° day⁻¹
        corr = A - A_car + B * np.sin(lat) ** 2 + C * np.sin(lat) ** 4  # ° day⁻¹
        corr = np.deg2rad(corr / 86400)  # rad s⁻¹
        return corr

    @staticmethod
    def recenter_crpix_in_header(hdr):
        w = WCS(hdr)
        if "ZNAXIS1" in hdr:
            naxis1 = hdr["ZNAXIS1"]
            naxis2 = hdr["ZNAXIS2"]
        else:
            naxis1 = hdr["NAXIS1"]
            naxis2 = hdr["NAXIS2"]
        x_mid = (naxis1 - 1) / 2
        y_mid = (naxis2 - 1) / 2
        try:
            lon_mid, lat_mid = w.pixel_to_world(np.array([x_mid]), np.array([y_mid]))
        except:
            coords = w.pixel_to_world(np.array([x_mid]), np.array([y_mid]))
            lon_mid = coords.Tx
            lat_mid = coords.Ty
        lon_mid = lon_mid[0].to(hdr["CUNIT1"]).value
        lat_mid = lat_mid[0].to(hdr["CUNIT2"]).value
        hdr["CRVAL1"] = lon_mid
        hdr["CRVAL2"] = lat_mid
        hdr["CRPIX1"] = (naxis1 + 1) / 2
        hdr["CRPIX2"] = (naxis2 + 1) / 2

    @staticmethod
    def extend_fov_in_header(hdr, extension_naxis1_pixels: int, extension_naxis2_pixels: int, sunpy=False):

        # first, recenter crpix.
        w = WCS(hdr)
        if "ZNAXIS1" in hdr:
            naxis1 = hdr["ZNAXIS1"]
            naxis2 = hdr["ZNAXIS2"]
        else:
            naxis1 = hdr["NAXIS1"]
            naxis2 = hdr["NAXIS2"]
        x_mid = (naxis1 - 1) / 2
        y_mid = (naxis2 - 1) / 2
        if sunpy:
            coords = w.pixel_to_world(np.array([x_mid]), np.array([y_mid]))
            lon_mid = coords.Tx.to(hdr["CUNIT1"]).value[0]
            lat_mid = coords.Ty.to(hdr["CUNIT1"]).value[0]

        else:
            lon_mid, lat_mid = w.pixel_to_world(np.array([x_mid]), np.array([y_mid]))
            lon_mid = lon_mid[0].to(hdr["CUNIT1"]).value
            lat_mid = lat_mid[0].to(hdr["CUNIT2"]).value
        hdr["CRVAL1"] = lon_mid
        hdr["CRVAL2"] = lat_mid
        hdr["CRPIX1"] = (naxis1 + 1) / 2
        hdr["CRPIX2"] = (naxis2 + 1) / 2
        #second extend naxis

        hdr["NAXIS1"] += extension_naxis1_pixels
        hdr["NAXIS2"] += extension_naxis2_pixels

        hdr["CRPIX1"] += extension_naxis1_pixels/2
        hdr["CRPIX2"] += extension_naxis2_pixels/2






