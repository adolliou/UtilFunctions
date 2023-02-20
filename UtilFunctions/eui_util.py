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