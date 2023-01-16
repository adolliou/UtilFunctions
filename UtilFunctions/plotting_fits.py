from matplotlib import pyplot as plt
import numpy as np
from . import eui_util
from UtilFunctions.common_util import CommonUtil

from astropy.wcs import WCS
import astropy.units as u
from matplotlib.gridspec import GridSpec
from astropy.visualization import ImageNormalize, LogStretch
import matplotlib.patches as patches


class PlotFits:

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
        return_im = False
        if fig is None:
            fig = plt.figure()
            return_im = True
        if ax is None:
            ax = fig.add_subplot()
        if norm is None:
            norm = ImageNormalize(stretch=LogStretch(5))
        im = ax.imshow(image_on_regular_grid, origin="lower", interpolation="none", norm=norm,
                       extent=[longitude_grid[0, 0].to(u.arcsec).value, longitude_grid[-1, -1].to(u.arcsec).value,
                               latitude_grid[0, 0].to(u.arcsec).value, latitude_grid[-1, -1].to(u.arcsec).value])
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

        return_im = True
        if fig is None:
            fig = plt.figure()
            return_im = False
        if ax is None:
            ax = fig.add_subplot()
        if norm is None:
            norm = ImageNormalize(stretch=LogStretch(5))
        im = ax.imshow(image_main_cut, origin="lower", interpolation="none", norm=norm,
                       extent=[longitude_main[0, 0].to(u.arcsec).value, longitude_main[-1, -1].to(u.arcsec).value,
                               latitude_main[0, 0].to(u.arcsec).value, latitude_main[-1, -1].to(u.arcsec).value])
        if levels is None:
            max_small = np.nanmax(image_contour_cut)
            levels = [0.5 * max_small]
        ax.contour(image_contour_cut, levels=levels, origin='lower', linewidths=0.5, colors='w',
                   extent=[longitude_main[0, 0].to(u.arcsec).value, longitude_main[-1, -1].to(u.arcsec).value,
                           latitude_main[0, 0].to(u.arcsec).value, latitude_main[-1, -1].to(u.arcsec).value])
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
    def _build_regular_grid(longitude, latitude):
        dlon = np.abs((longitude[1, 1] - longitude[0, 0]).to(u.deg).value)
        dlat = np.abs((latitude[1, 1] - latitude[0, 0]).to(u.deg).value)
        longitude_grid, latitude_grid = np.meshgrid(
            np.arange(np.min(CommonUtil.ang2pipi(longitude).to(u.deg).value),
                      np.max(CommonUtil.ang2pipi(longitude).to(u.deg).value), dlon),
            np.arange(np.min(CommonUtil.ang2pipi(latitude).to(u.deg).value),
                      np.max(CommonUtil.ang2pipi(latitude).to(u.deg).value), dlat))

        longitude_grid = longitude_grid * u.deg
        latitude_grid = latitude_grid * u.deg

        return longitude_grid, latitude_grid
