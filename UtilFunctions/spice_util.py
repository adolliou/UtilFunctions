import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


class SpiceUtil:

    @staticmethod
    def slit_pxl(header):

        """ Compute the first and last pixel of the slit from a FITS header """
        ybin = header['NBIN2']
        h_detector = 1024 / ybin
        if header['DETECTOR'] == 'SW':
            h_slit = 600 / ybin
        elif header['DETECTOR'] == 'LW':
            h_slit = 626 / ybin
        else:
            raise ValueError(f"unknown detector: {header['DETECTOR']}")
        slit_beg = (h_detector - h_slit) / 2
        slit_end = h_detector - slit_beg
        slit_beg = slit_beg - header['PXBEG2'] / ybin + 1
        slit_end = slit_end - header['PXBEG2'] / ybin + 1
        slit_beg = int(np.ceil(slit_beg))
        slit_end = int(np.floor(slit_end))
        return slit_beg, slit_end

    @staticmethod
    def vertical_edges_limits(header):
        iymin, iymax = SpiceUtil.slit_pxl(header)
        iymin += int(20 / header['NBIN2'])
        iymax -= int(20 / header['NBIN2'])
        return iymin, iymax

    @staticmethod
    def create_intensity_map(path_to_l3, index_amplitude=0, index_width=2):
        hdul = fits.open(path_to_l3)
        hdu_results = hdul[0]
        data = hdu_results.data.copy()
        hdr = hdu_results.header.copy()
        hdul.close()

        missing = hdr["ANA_MISS"]
        condition_missing = np.array((data[:, :, index_amplitude] == missing) | (data[:, :, index_width] == missing),
                                     dtype="bool")

        intensity = data[:, :, index_amplitude] * data[:, :, index_width] * np.sqrt(2 * np.pi)
        intensity[condition_missing] = np.nan

        return hdr, intensity

    @staticmethod
    def _data_treatment_spice(self, hdr, data):
        w_small = WCS(hdr)
        w2 = w_small.deepcopy()
        w2.wcs.pc[3, 0] = 0
        w2.wcs.pc[3, 1] = 0
        w_xyt = w2.dropaxis(2)
        # TODO 1) select slit positions 2) correct for solar rotation
        w_xy_small = w_xyt.dropaxis(2)
        return w_xy_small

    @staticmethod
    def extract_spice_coordinates_l3(hdr, return_type='xy'):
        w_small = WCS(hdr)
        w2 = w_small.deepcopy()

        w2.wcs.pc[3, 0] = 0

        w_xyt = w2.dropaxis(0)
        if return_type == 'xy':
            w2.wcs.pc[3, 1] = 0
            w_xy = w_xyt.dropaxis(2)
            idx_lon = np.where(np.array(w_xy.wcs.ctype, dtype="str") == "HPLN-TAN")[0][0]
            idx_lat = np.where(np.array(w_xy.wcs.ctype, dtype="str") == "HPLT-TAN")[0][0]
            x_small, y_small = np.meshgrid(np.arange(w_xy.pixel_shape[idx_lon]),
                                           np.arange(w_xy.pixel_shape[idx_lat]),
                                           indexing='ij')  # t dépend de x,
            # should reproject on a new coordinate grid first : suppose slits at the same time :
            longitude_small, latitude_small = w_xy.pixel_to_world(x_small, y_small)
            return longitude_small, latitude_small
        elif return_type == 'xyt':
            idx_lon = np.where(np.array(w_xyt.wcs.ctype, dtype="str") == "HPLN-TAN")[0][0]
            idx_lat = np.where(np.array(w_xyt.wcs.ctype, dtype="str") == "HPLT-TAN")[0][0]
            idx_utc = np.where(np.array(w_xyt.wcs.ctype, dtype="str") == "UTC")[0][0]
            x_small, y_small, z_small = np.meshgrid(np.arange(w_xyt.pixel_shape[idx_lon]),
                                                    np.arange(w_xyt.pixel_shape[idx_lat]),
                                                    np.arange(w_xyt.pixel_shape[idx_utc]),
                                                    indexing='ij')  # t dépend de x,
            longitude_small, latitude_small, utc_small = w_xyt.pixel_to_world(x_small, y_small, z_small)
            return longitude_small, latitude_small, utc_small
