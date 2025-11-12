from matplotlib import pyplot as plt
import numpy as np
from . import eui_util
from UtilFunctions.common_util import CommonUtil
from astropy.wcs import WCS
import astropy.units as u
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches
import matplotlib.colors as colors

from astropy.visualization import ImageNormalize, AsymmetricPercentileInterval, SqrtStretch, LinearStretch, LogStretch


class CmapUtil:

    @staticmethod
    def create_cdict(r, g, b):
        """
        Create the color tuples in the correct format.
        """
        i = np.linspace(0, 1, r.size)
        cdict = {name: list(zip(i, el / 255.0, el / 255.0))
                 for el, name in [(r, 'red'), (g, 'green'), (b, 'blue')]}
        return cdict

    @staticmethod
    def get_idl3(path_idl3: str):
        # The following values describe color table 3 for IDL (Red Temperature)
        return np.loadtxt(path_idl3, delimiter=',')

    @staticmethod
    def _cmap_from_rgb(r, g, b, name):
        cdict = CmapUtil.create_cdict(r, g, b)
        return colors.LinearSegmentedColormap(name, cdict)

    @staticmethod
    def create_aia_wave_dict(path_idl3: str):
        idl_3 = CmapUtil.get_idl3(path_idl3)
        r0, g0, b0 = idl_3[:, 0], idl_3[:, 1], idl_3[:, 2]

        c0 = np.arange(256, dtype='f')
        c1 = (np.sqrt(c0) * np.sqrt(255.0)).astype('f')
        c2 = (np.arange(256) ** 2 / 255.0).astype('f')
        c3 = ((c1 + c2 / 2.0) * 255.0 / (c1.max() + c2.max() / 2.0)).astype('f')

        aia_wave_dict = {
            1600 * u.angstrom: (c3, c3, c2),
            1700 * u.angstrom: (c1, c0, c0),
            4500 * u.angstrom: (c0, c0, b0 / 2.0),
            94 * u.angstrom: (c2, c3, c0),
            131 * u.angstrom: (g0, r0, r0),
            171 * u.angstrom: (r0, c0, b0),
            193 * u.angstrom: (c1, c0, c2),
            211 * u.angstrom: (c1, c0, c3),
            304 * u.angstrom: (r0, g0, b0),
            335 * u.angstrom: (c2, c0, c1)
        }
        return aia_wave_dict

    @staticmethod
    def solohri_lya1216_color_table(path_idl3: str):
        solohri_lya1216 = CmapUtil.get_idl3(path_idl3)
        solohri_lya1216[:, 2] = solohri_lya1216[:, 0] * np.linspace(0, 1, 256)
        return CmapUtil._cmap_from_rgb(*solohri_lya1216.T, 'SolO EUI HRI Lyman Alpha')

    @staticmethod
    def aia_color_table(wavelength: u.angstrom, path_idl3: str):
        """
        Returns one of the fundamental color tables for SDO AIA images.

        Based on aia_lct.pro part of SDO/AIA on SSWIDL written by Karel
        Schrijver (2010/04/12).

        Parameters
        ----------
        wavelength : `~astropy.units.quantity`
            Wavelength for the desired AIA color table.
        """
        aia_wave_dict = CmapUtil.create_aia_wave_dict(path_idl3)
        try:
            r, g, b = aia_wave_dict[wavelength]
        except KeyError:
            raise ValueError("Invalid AIA wavelength. Valid values are "
                             "1600,1700,4500,94,131,171,193,211,304,335.")

        return CmapUtil._cmap_from_rgb(r, g, b, 'SDO AIA {:s}'.format(str(wavelength)))


class PlotFits:
    @staticmethod
    def get_range(data, stre='linear', imax=99.5, imin=2, a=1000):
        """
        :param data:
        :param stretch: 'sqrt', 'log', or 'linear' (default)
        :return: norm
        """
        if np.isnan(data).sum() == data.size:
            return None

        isnan = np.isnan(data)
        data = data[~isnan]
        do = False
        if imax > 100:
            vmin, vmax = AsymmetricPercentileInterval(imin, 100).get_limits(data)
            vmax = vmax * imax/100
        else:
            vmin, vmax = AsymmetricPercentileInterval(imin, imax).get_limits(data)

        #    print('Vmin:', vmin, 'Vmax', vmax)
        if stre == 'linear':
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
        elif stre == 'sqrt':
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())
        elif stre == 'log':
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch(a))
        else:
            raise ValueError('Bad stre value: either linear, sqrt or log')
        return norm
    @staticmethod
    def plot_fov_rectangle(data, slc=None, path_save=None, show=True, plot_colorbar=True, norm=None, angle=0):
        fig = plt.figure()
        ax = fig.add_subplot()
        if norm is None:
            norm = ImageNormalize(stretch=LogStretch(5))
        PlotFits.plot_fov(data=data, show=False, fig=fig, ax=ax, norm=norm)
        rect = patches.Rectangle((slc[1].start, slc[0].start),
                                 slc[1].stop - slc[1].start, slc[0].stop - slc[0].start, linewidth=1,
                                 edgecolor='r', facecolor='none', angle=angle)
        ax.add_patch(rect)
        if show:
            fig.show()
        if path_save is not None:
            fig.savefig(path_save)

    @staticmethod
    def plot_fov(data, slc=None, path_save=None, show=True, plot_colorbar=True, fig=None, ax=None, norm=None):
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot()
        if norm is None:
            norm = ImageNormalize(stretch=LogStretch(5))
        if slc is not None:
            im = ax.imshow(data[slc[0], slc[1]], origin="lower", interpolation="none", norm=norm)
        else:
            im = ax.imshow(data, origin="lower", interpolation="None", norm=norm)
        if plot_colorbar:
            fig.colorbar(im, label="DN/s")
        if show:
            fig.show()
        if path_save is not None:
            fig.savefig(path_save)

    @staticmethod
    def simple_plot(hdr_main, data_main, path_save=None, show=True, ax=None, fig=None, norm=None,
                    show_xlabel=True, show_ylabel=True, plot_colorbar=True):

        longitude, latitude, dsun = eui_util.EUIUtil.extract_EUI_coordinates(hdr_main)
        longitude_grid, latitude_grid = PlotFits._build_regular_grid(longitude=longitude, latitude=latitude)
        w = WCS(hdr_main)
        x, y = w.world_to_pixel(longitude_grid, latitude_grid)
        image_on_regular_grid = CommonUtil.interpol2d(data_main, x=x, y=y, fill=-32762, order=1)
        image_on_regular_grid[image_on_regular_grid == -32762] = np.nan

        dlon = (longitude_grid[1, 1] - longitude_grid[0, 0]).to("arcsec").value
        dlat = (latitude_grid[1, 1] - latitude_grid[0, 0]).to("arcsec").value

        return_im = False
        if fig is None:
            fig = plt.figure()
            return_im = True
        if ax is None:
            ax = fig.add_subplot()
        if norm is None:
            norm = ImageNormalize(stretch=LogStretch(5))
        im = ax.imshow(image_on_regular_grid, origin="lower", interpolation="none", norm=norm,
                       extent=[longitude_grid[0, 0].to(u.arcsec).value - 0.5*dlon,
                               longitude_grid[-1, -1].to(u.arcsec).value + 0.5*dlon,
                               latitude_grid[0, 0].to(u.arcsec).value - 0.5*dlat,
                               latitude_grid[-1, -1].to(u.arcsec).value + 0.5*dlat])
        # im = ax.imshow(data_main, origin="lower", interpolation="none", norm=norm,)

        if show_xlabel:
            ax.set_xlabel("Solar-X [arcsec]")
        if show_ylabel:
            ax.set_ylabel("Solar-Y [arcsec]")
        if plot_colorbar:
            fig.colorbar(im, label=hdr_main["BUNIT"])
        if show:
            fig.show()
        if path_save is not None:
            fig.savefig(path_save)
        if return_im:
            return im

    @staticmethod
    def contour_plot(hdr_main, data_main, hdr_contour, data_contour, path_save=None, show=True, levels=None,
                     ax=None, fig=None, norm=None, show_xlabel=True, show_ylabel=True, plot_colorbar=True):
        longitude_main, latitude_main, dsun = eui_util.EUIUtil.extract_EUI_coordinates(hdr_contour)
        longitude_grid, latitude_grid = PlotFits._build_regular_grid(longitude=longitude_main,
                                                                     latitude=latitude_main)

        w_xy_main = WCS(hdr_main)
        x_small, y_small = w_xy_main.world_to_pixel(longitude_grid, latitude_grid)
        image_main_cut = CommonUtil.interpol2d(np.array(data_main,
                                                        dtype=np.float64), x=x_small, y=y_small,
                                               order=1, fill=-32768)
        image_main_cut[image_main_cut == -32768] = np.nan

        w_xy_contour = WCS(hdr_contour)
        x_contour, y_contour = w_xy_contour.world_to_pixel(longitude_grid, latitude_grid)
        image_contour_cut = CommonUtil.interpol2d(np.array(data_contour, dtype=np.float64),
                                                  x=x_contour, y=y_contour,
                                                  order=1, fill=-32768)
        image_contour_cut[image_contour_cut == -32768] = np.nan
        dlon = (longitude_grid[1, 1] - longitude_grid[0, 0]).to("arcsec").value
        dlat = (latitude_grid[1, 1] - latitude_grid[0, 0]).to("arcsec").value

        return_im = True
        if fig is None:
            fig = plt.figure()
            return_im = False
        if ax is None:
            ax = fig.add_subplot()
        if norm is None:
            norm = ImageNormalize(stretch=LogStretch(5))
        im = ax.imshow(image_main_cut, origin="lower", interpolation="none", norm=norm,
                       extent=[longitude_grid[0, 0].to(u.arcsec).value - 0.5*dlon,
                               longitude_grid[-1, -1].to(u.arcsec).value + 0.5*dlon,
                               latitude_grid[0, 0].to(u.arcsec).value - 0.5*dlat,
                               latitude_grid[-1, -1].to(u.arcsec).value + 0.5*dlat])
        if levels is None:
            max_small = np.nanmax(image_contour_cut)
            levels = [0.5 * max_small]
        ax.contour(image_contour_cut, levels=levels, origin='lower', linewidths=0.5, colors='w',
                   extent=[longitude_grid[0, 0].to(u.arcsec).value - 0.5 * dlon,
                           longitude_grid[-1, -1].to(u.arcsec).value + 0.5 * dlon,
                           latitude_grid[0, 0].to(u.arcsec).value - 0.5 * dlat,
                           latitude_grid[-1, -1].to(u.arcsec).value + 0.5 * dlat])
        if show_xlabel:
            ax.set_xlabel("Solar-X [arcsec]")
        if show_ylabel:
            ax.set_ylabel("Solar-Y [arcsec]")
        if plot_colorbar:
            fig.colorbar(im, label=hdr_main["BUNIT"])
        if show:
            fig.show()
        if path_save is not None:
            fig.savefig(path_save)
        if return_im:
            return im

    @staticmethod
    def compare_plot(hdr_main, data_main, hdr_contour_1, data_contour_1,
                     hdr_contour_2, data_contour_2, norm, path_save=None, show=True, levels=None, ):

        if (norm.vmin is None) or (norm.vmax is None):
            raise ValueError("Must explicit vmin and vmax in norm, so that the cbar is the same for both figures.")

        fig = plt.figure(figsize=(10, 5))
        gs = GridSpec(1, 3, width_ratios=[1, 1, 0.1], wspace=0.3)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax_cbar = fig.add_subplot(gs[2])

        im = PlotFits.contour_plot(hdr_main=hdr_main, data_main=data_main, plot_colorbar=False,
                                   hdr_contour=hdr_contour_1, data_contour=data_contour_1,
                                   path_save=None, show=False, levels=levels, fig=fig, ax=ax1, norm=norm)

        im = PlotFits.contour_plot(hdr_main=hdr_main, data_main=data_main, show_ylabel=False, plot_colorbar=False,
                                   hdr_contour=hdr_contour_2, data_contour=data_contour_2,
                                   path_save=None, show=False, levels=levels, fig=fig, ax=ax2, norm=norm)

        fig.colorbar(im, cax=ax_cbar, label=hdr_main["BUNIT"])

        if show:
            fig.show()
        if path_save is not None:
            fig.savefig(path_save)

    @staticmethod
    def build_regular_grid(longitude, latitude, lonlims=None, latlims=None, apply_ang2pipi=True):
        x = np.abs((longitude[0, 1] - longitude[0, 0]).to("deg").value)
        y = np.abs((latitude[0, 1] - latitude[0, 0]).to("deg").value)
        dlon = np.sqrt(x**2 + y**2)

        x = np.abs((longitude[1, 0] - longitude[0, 0]).to("deg").value)
        y = np.abs((latitude[1, 0] - latitude[0, 0]).to("deg").value)
        dlat = np.sqrt(x**2 + y**2)
        if apply_ang2pipi:
            ang2pipi = CommonUtil.ang2pipi
        else:
            def ang2pipi(a):
                return a

        longitude1D = np.arange(np.min(ang2pipi(longitude).to(u.deg).value),
                                np.max(ang2pipi(longitude).to(u.deg).value), dlon)
        latitude1D = np.arange(np.min(ang2pipi(latitude).to(u.deg).value),
                               np.max(ang2pipi(latitude).to(u.deg).value), dlat)
        if (lonlims is not None) or (latlims is not None):
            longitude1D = longitude1D[(longitude1D > ang2pipi(lonlims[0]).to("deg").value) &
                                      (longitude1D < ang2pipi(lonlims[1]).to("deg").value)]
            latitude1D = latitude1D[(latitude1D > ang2pipi(latlims[0]).to("deg").value) &
                                    (latitude1D < ang2pipi(latlims[1]).to("deg").value)]
        longitude_grid, latitude_grid = np.meshgrid(longitude1D, latitude1D)

        longitude_grid = longitude_grid * u.deg
        latitude_grid = latitude_grid * u.deg
        dlon = dlon * u.deg
        dlat = dlat * u.deg
        return longitude_grid, latitude_grid, dlon, dlat

    @staticmethod
    def extend_regular_grid(longitude_grid, latitude_grid, delta_longitude, delta_latitude):
        x = np.abs((longitude_grid[0, 1] - longitude_grid[0, 0]).to("deg").value)
        y = np.abs((latitude_grid[0, 1] - latitude_grid[0, 0]).to("deg").value)
        dlon = np.sqrt(x**2 + y**2)

        x = np.abs((longitude_grid[1, 0] - longitude_grid[0, 0]).to("deg").value)
        y = np.abs((latitude_grid[1, 0] - latitude_grid[0, 0]).to("deg").value)
        dlat = np.sqrt(x**2 + y**2)

        delta_longitude_deg = CommonUtil.ang2pipi(delta_longitude).to("deg").value
        delta_latitude_deg = CommonUtil.ang2pipi(delta_latitude).to("deg").value

        longitude1D = np.arange(np.min(CommonUtil.ang2pipi(longitude_grid).to(u.deg).value - 0.5 * delta_longitude_deg),
                                np.max(CommonUtil.ang2pipi(longitude_grid).to(u.deg).value) + 0.5 * delta_longitude_deg,
                                dlon)
        latitude1D = np.arange(np.min(CommonUtil.ang2pipi(latitude_grid).to(u.deg).value - 0.5 * delta_latitude_deg),
                               np.max(CommonUtil.ang2pipi(latitude_grid).to(u.deg).value) + 0.5 * delta_latitude_deg,
                               dlat)

        longitude_grid_ext, latitude_grid_ext = np.meshgrid(longitude1D, latitude1D)
        longitude_grid_ext = longitude_grid_ext * u.deg
        latitude_grid_ext = latitude_grid_ext * u.deg

        return longitude_grid_ext, latitude_grid_ext


