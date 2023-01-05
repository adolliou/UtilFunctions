import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from plot import PlottingFunctions


def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # non non pas du tout copié collé de la doc python.
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
    binwidth2 = 0.05
    binwidth1 = 0.2

    xmax = np.max(np.abs(x))
    ymax = np.max(np.abs(y))

    # lim1 = (int(xmax/binwidth1) + 1) * binwidth1
    # lim2 = (int(ymax/binwidth2) + 1) * binwidth2
    lim1 = 2 + binwidth1 * 0.5
    lim2 = 1 + binwidth2
    ax.scatter(x, y)

    bins1 = np.arange(-lim1, lim1, binwidth1)
    bins2 = np.arange(0.2, lim2 + binwidth2, binwidth2)

    ax_histx.hist(x, bins=bins1, )
    ax_histy.hist(y, bins=bins2, orientation='horizontal', )

    return bins1, bins2


def scatter_hist_2(x, y, ax, ax_histx, ax_histy, channel: int, limit_corr: float, binstyle='surface_max_corr',
                   bins1=None, bins2=None):
    # non non pas du tout copié collé de la doc python.
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # print('voilà : %i' %(len(x)))
    # print(len(y))
    C1 = np.isfinite(x)
    C2 = np.isfinite(y)
    C = C1 & C2
    x = x[C]
    y = y[C]

    # the scatter plot:
    # ax.scatter(x, y)

    # now determine nice limits by hand:

    xmax = np.max(np.abs(x))
    ymax = np.max(np.abs(y))
    ymin = np.min(np.abs(y))
    if binstyle == 'surface_max_corr':

        if bins1 is None:
            binwidth1 = 2
            lim1 = xmax + binwidth1
            bins1 = np.arange(0, lim1 + binwidth1, binwidth1)

        if bins2 is None:
            ymax = np.max(np.abs(y))
            binwidth2 = ymax / 20
            lim2 = 1 + binwidth2
            bins2 = np.arange(limit_corr, lim2 + binwidth2, binwidth2)

    if binstyle == 'surface_tau':

        if bins1 is None:
            binwidth1 = 2
            lim1 = xmax + binwidth1
            bins1 = np.arange(0, lim1 + binwidth1, binwidth1)

        if bins2 is None:
            binwidth2 = 0.2
            lim2 = 2 - binwidth2 * 0.5
            bins2 = np.arange(-lim2, lim2, binwidth2)

    if binstyle == 'duration_maxcorr':
        if bins1 is None:
            binwidth1 = 5
            lim1 = xmax + binwidth1
            bins1 = np.arange(0, lim1 + binwidth1, binwidth1)
        if bins2 is None:
            binwidth2 = ymax / 20
            lim2 = 1 + binwidth2
            bins2 = np.arange(limit_corr, lim2 + binwidth2, binwidth2)

    if binstyle == 'duration_tau':

        if bins1 is None:
            binwidth1 = 5
            lim1 = xmax + binwidth1
            bins1 = np.arange(0, lim1 + binwidth1, binwidth1)
        if bins2 is None:
            binwidth2 = 0.2
            lim2 = 2 - binwidth2 * 0.5
            bins2 = np.arange(-lim2, lim2, binwidth2)

    if binstyle == 'height_tau':
        binwidth1 = xmax / 20.0
        binwidth2 = 0.2

        lim1 = xmax + binwidth1
        lim2 = 2 - binwidth2 * 0.5

        bins1 = np.arange(0, lim1 + binwidth1, binwidth1)
        bins2 = np.arange(-lim2, lim2, binwidth2)

    if binstyle == 'height_duration':
        binwidth1 = xmax / 20.0
        binwidth2 = 5

        lim1 = xmax + binwidth1
        lim2 = ymax + binwidth2

        bins1 = np.arange(0, lim1 + binwidth1, binwidth1)
        bins2 = np.arange(0, lim2 + binwidth2, binwidth2)

    if binstyle == 'surface_signal':

        binwidth1 = 2
        binwidth2 = ymax / 30.0

        lim1 = xmax + binwidth1
        lim2 = ymax + binwidth2

        bins1 = np.arange(0, lim1 + binwidth1, binwidth1)
        if channel == 335:
            bins2 = np.arange(ymin - binwidth2, lim2, binwidth2)

        else:
            bins2 = np.arange(0, lim2, binwidth2)

    ax.hist2d(x, y, bins=(bins1, bins2), cmap=plt.cm.Reds, norm=LogNorm())

    H, xedge, yedge = np.histogram2d(x, y, bins=(bins1, bins2), )
    H = np.swapaxes(H, axis1=0, axis2=1)
    # ax.hist2d(x, y, bins = 30, cmap=plt.cm.jet)
    xcoord = np.zeros(len(xedge) - 1)
    ycoord = np.zeros(len(yedge) - 1)

    for xx in range(len(xcoord)):
        xcoord[xx] = 0.5 * (xedge[xx] + xedge[xx + 1])
    for yy in range(len(ycoord)):
        ycoord[yy] = 0.5 * (yedge[yy] + yedge[yy + 1])

    levels = np.arange(int(0.15 * np.max(H)), int(0.7 * np.max(H)), int(0.15 * np.max(H)))
    cs = ax.contour(xcoord, ycoord, H, levels=levels, linewidths=1, origin='lower', colors='k')

    ax.clabel(cs, inline=True, fontsize=10)

    ax_histx.hist(x, bins=bins1, histtype='step', color='r')
    ax_histy.hist(y, bins=bins2, orientation='horizontal', histtype='step', color='r')

    return bins1, bins2


def scatter_hist_3(x1, y1, x2, y2, ax, channel: int, limit_corr=0, binstyle='maxcorr_signal',
                   ax_histx=None, ax_histy=None, vmin=1, vmax=100, limit_xcorr=None, limit_hist1d_y=None,
                   color='Reds', bins1=None, bins2=None, limit_hist1d_x=None, show_xhist_label=True):
    color1D = {'Reds': 'red',
               'Greens': 'green',
               'Oranges': 'orange'}
    # no labels
    if ax_histx is not None:
        ax_histx.tick_params(axis="x", labelbottom=False)
    if ax_histy is not None:
        ax_histy.tick_params(axis="y", labelleft=False)

    # print('voilà : %i' %(len(x)))
    # print(len(y))
    C1 = np.isfinite(x1)
    C2 = np.isfinite(y1)
    C = C1 & C2
    x1 = x1[C]
    y1 = y1[C]

    C1 = np.isfinite(x2)
    C2 = np.isfinite(y2)
    # remove dark pixels for HRI174

    C = C1 & C2
    x2 = x2[C]
    y2 = y2[C]

    # plt.colorbar(im)

    # ax.hist2d(x2, y2, bins=15,)

    # ax.hist2d(x, y, bins = 30, cmap=plt.cm.jet)

    # the scatter plot:
    # ax.scatter(x, y)

    # now determine nice limits by hand:

    if binstyle == 'maxcorr_signal':

        if bins1 is None:
            xmax1 = np.max(np.abs(x1))
            binwidth1 = xmax1 / 30
            lim1 = 1.2
            bins1 = np.arange(limit_corr, lim1 + binwidth1, binwidth1)
        if bins2 is None:
            ymax1 = np.max(np.abs(y1))
            ymin1 = np.min(np.abs(y1))
            binwidth2 = ymax1 / 35
            lim2 = ymax1 - 0.5 * binwidth2
            print(binwidth2)
            print(lim2)
            if channel == 335:
                bins2 = np.arange(ymin1 - binwidth2, lim2, binwidth2)
            else:
                bins2 = np.arange(binwidth2 * 0.5, lim2, binwidth2)

        H1, xedge1, yedge1 = np.histogram2d(x1, y1, bins=(bins1, bins2), )
        H1 = np.swapaxes(H1, axis1=0, axis2=1)

        H, xedge, yedge = np.histogram2d(x2, y2, bins=(bins1, bins2), )
        H = np.swapaxes(H, axis1=0, axis2=1)
        # ax.hist2d(x, y, bins = 30, cmap=plt.cm.jet)

        xcoord1 = np.zeros(len(xedge1) - 1)
        ycoord1 = np.zeros(len(yedge1) - 1)

        xcoord = np.zeros(len(xedge) - 1)
        ycoord = np.zeros(len(yedge) - 1)

        for xx in range(len(xcoord)):
            xcoord[xx] = 0.5 * (xedge[xx] + xedge[xx + 1])
        for yy in range(len(ycoord)):
            ycoord[yy] = 0.5 * (yedge[yy] + yedge[yy + 1])

        for xx in range(len(xcoord1)):
            xcoord1[xx] = 0.5 * (xedge1[xx] + xedge1[xx + 1])
        for yy in range(len(ycoord1)):
            ycoord1[yy] = 0.5 * (yedge1[yy] + yedge1[yy + 1])
        # ax.plot(x2, y2, '.', markersize = 0.2, color='gray')

        h = ax.hist2d(x1, y1, bins=(bins1, bins2), cmap=color, norm=LogNorm(vmin=vmin, vmax=vmax), label='campf')
        # plt.legend()

        # levels1 = np.arange(0.2 * np.max(H1), 0.8 * np.max(H1), 0.2 * np.max(H1))

        # levels = np.arange(0.2 * np.max(H), 0.8 * np.max(H), 0.2 * np.max(H))
        levels = [0.2 * np.max(H), 0.4 * np.max(H), 0.6 * np.max(H), 0.8 * np.max(H)]
        levels1 = [0.2 * np.max(H1), 0.4 * np.max(H1), 0.6 * np.max(H1), 0.8 * np.max(H1)]

        cs = ax.contour(xcoord, ycoord, H, levels=levels, linewidths=1,
                        origin='lower', colors=('#b0aaff', '#8880ff', '#6055ff', '#392bff',))

        # ax.clabel(cs, inline=True, fontsize=10)
        cs = ax.contour(xcoord1, ycoord1, H1, levels=levels1, linewidths=1,
                        origin='lower', colors=('#ff0000', '#db0000', '#b80000', '#960000',))
        # cs = ax.contour(xcoord1, ycoord1, H1, levels=levels1, linewidths=1, origin='lower', colors='k')
        # ax.clabel(cs, inline=True, fontsize=10)
    elif binstyle == 'maxcorr_tau':
        ymax2 = np.max(np.abs(y2))

        if bins1 is None:
            binwidth1 = 0.2
            lim1 = 2 - binwidth1 * 0.5
            bins1 = np.arange(-lim1, lim1 + binwidth1, binwidth1)
        if bins2 is None:
            binwidth2 = ymax2 / 20
            lim2 = 1 + binwidth2
            print(lim2)
            print(binwidth2)
            bins2 = np.arange(0, lim2 + binwidth2, binwidth2)

        H1, xedge1, yedge1 = np.histogram2d(x1, y1, bins=(bins1, bins2), )
        H1 = np.swapaxes(H1, axis1=0, axis2=1)

        H, xedge, yedge = np.histogram2d(x2, y2, bins=(bins1, bins2))
        H = np.swapaxes(H, axis1=0, axis2=1)
        # ax.hist2d(x, y, bins = 30, cmap=plt.cm.jet)

        xcoord1 = np.zeros(len(xedge1) - 1)
        ycoord1 = np.zeros(len(yedge1) - 1)

        xcoord = np.zeros(len(xedge) - 1)
        ycoord = np.zeros(len(yedge) - 1)

        for xx in range(len(xcoord)):
            xcoord[xx] = 0.5 * (xedge[xx] + xedge[xx + 1])
        for yy in range(len(ycoord)):
            ycoord[yy] = 0.5 * (yedge[yy] + yedge[yy + 1])

        for xx in range(len(xcoord1)):
            xcoord1[xx] = 0.5 * (xedge1[xx] + xedge1[xx + 1])
        for yy in range(len(ycoord1)):
            ycoord1[yy] = 0.5 * (yedge1[yy] + yedge1[yy + 1])
        # ax.plot(x2, y2, '.', markersize = 0.2, color='gray')

        # ax.plot(x2, y2, '.', markersize = 0.2, color='gray')

        levels1 = np.arange(int(0.13 * np.max(H1)), int(0.7 * np.max(H1)), int(0.15 * np.max(H1)))

        levels = np.arange(int(0.13 * np.max(H)), int(0.7 * np.max(H)), int(0.15 * np.max(H)))

        cs = ax.contour(xcoord, ycoord, H, levels=levels, linewidths=1, origin='lower', colors='blue')

        ax.clabel(cs, inline=True, fontsize=10)

        cs = ax.contour(xcoord1, ycoord1, H1, levels=levels1, linewidths=1, origin='lower', colors='k')
        ax.clabel(cs, inline=True, fontsize=10)

        h = ax.hist2d(x1, y1, bins=(bins1, bins2), cmap=color, norm=LogNorm(vmin=vmin, vmax=vmax))

    elif binstyle == 'tau_signal':

        if bins1 is None:
            binwidth1 = 0.2
            lim1 = 2 - binwidth1 * 0.5
            bins1 = np.arange(-lim1, lim1, binwidth1)
        if bins2 is None:
            ymin1 = np.min(np.abs(y1))
            ymax1 = np.max(np.abs(y1))
            binwidth2 = ymax1 / 20
            lim2 = ymax1
            if channel == 335:
                bins2 = np.arange(ymin1 - binwidth2, lim2, binwidth2)
            else:
                bins2 = np.arange(binwidth2 * 0.5, lim2, binwidth2)

        H1, xedge1, yedge1 = np.histogram2d(x1, y1, bins=(bins1, bins2), )
        H1 = np.swapaxes(H1, axis1=0, axis2=1)

        H, xedge, yedge = np.histogram2d(x2, y2, bins=(bins1, bins2))
        H = np.swapaxes(H, axis1=0, axis2=1)
        # ax.hist2d(x, y, bins = 30, cmap=plt.cm.jet)

        xcoord1 = np.zeros(len(xedge1) - 1)
        ycoord1 = np.zeros(len(yedge1) - 1)

        xcoord = np.zeros(len(xedge) - 1)
        ycoord = np.zeros(len(yedge) - 1)

        for xx in range(len(xcoord)):
            xcoord[xx] = 0.5 * (xedge[xx] + xedge[xx + 1])
        for yy in range(len(ycoord)):
            ycoord[yy] = 0.5 * (yedge[yy] + yedge[yy + 1])

        for xx in range(len(xcoord1)):
            xcoord1[xx] = 0.5 * (xedge1[xx] + xedge1[xx + 1])
        for yy in range(len(ycoord1)):
            ycoord1[yy] = 0.5 * (yedge1[yy] + yedge1[yy + 1])
        # ax.plot(x2, y2, '.', markersize = 0.2, color='gray')

        h = ax.hist2d(x1, y1, bins=(bins1, bins2), cmap='Reds', norm=LogNorm(vmin=vmin, vmax=vmax), label='campf')

        levels1 = np.arange(0.2 * np.max(H1), 0.8 * np.max(H1), 0.2 * np.max(H1))

        levels = np.arange(0.2 * np.max(H), 0.8 * np.max(H), 0.2 * np.max(H))

        cs = ax.contour(xcoord, ycoord, H, levels=levels, linewidths=1, origin='lower', colors='blue')

        ax.clabel(cs, inline=True, fontsize=10)

        cs = ax.contour(xcoord1, ycoord1, H1, levels=levels1, linewidths=1, origin='lower', colors='k')
        ax.clabel(cs, inline=True, fontsize=10)
    if ax_histx is not None:
        # ax_histx.hist(x2, bins=bins1,histtype='step', color = 'b', label='not cmpf', density = True)
        PlottingFunctions.normalised_histogram(x=x2, ax=ax_histx, bins=bins1, color='b')
        # ax_histx.legend(loc='center left', bbox_to_anchor=(1.035, 0.5))
    if ax_histy is not None:
        # ax_histy.hist(y2, bins=bins2, orientation='horizontal', color='b', histtype='step', label='not cmpf',
        #              density=True)
        PlottingFunctions.normalised_histogram(x=y2, ax=ax_histy, bins=bins2, color='b', orientation='horizontal')

    if binstyle == 'maxcorr_tau':
        if ax_histx is not None:
            ax_histx.hist(x1, bins=bins1, color=color1D[color], label='cmpf', histtype='step', density=True)
        if ax_histy is not None:
            ax_histy.hist(y1, bins=bins2, orientation='horizontal', color=color1D[color], label='cmpf', histtype='step',
                          density=True)
    else:
        if ax_histx is not None:
            # ax_histx.hist(x1, bins=bins1, color='r', label='cmpf', histtype='step', density=True)
            PlottingFunctions.normalised_histogram(x=x1, ax=ax_histx, bins=bins1, color='r')

        if ax_histy is not None:
            # ax_histy.hist(y1, bins=bins2, orientation='horizontal', color = 'r', label='cmpf', histtype='step',
            #          density = True)
            PlottingFunctions.normalised_histogram(x=y1, ax=ax_histy, bins=bins2, color='r', orientation='horizontal')

    # ax_histy.legend()
    if limit_hist1d_x is not None:
        ax_histx.set_ylim([0, limit_hist1d_x])
        if not (show_xhist_label):
            ax_histx.set_yticklabels([])
    if limit_hist1d_y is not None:
        ax_histy.set_xlim([0, limit_hist1d_y])

    return bins1, bins2, h


def scatter_hist_4(x1, y1, x2, y2, x3, y3, ax, ax_histx, ax_histy, channel: int, limit_corr: float,
                   binstyle='maxcorr_signal', color='Reds', bins1=None, bins2=None,
                   show_hist_x=True, show_hist_y=True, limit_hist1d_x=None, limit_hist1d_y=None,
                   show_background=True, vmin=1, vmax=100):
    color1D = {'Reds': 'red',
               'Greens': 'green',
               'Oranges': 'orange'}
    # no labels
    if ax_histx is not None:
        ax_histx.tick_params(axis="x", labelbottom=False)
    if ax_histy is not None:
        ax_histy.tick_params(axis="y", labelleft=False)

    # print('voilà : %i' %(len(x)))
    # print(len(y))

    for ii, (xx, yy) in enumerate(zip([x1, x2, x3], [y1, y2, y3])):
        C1 = np.isfinite(xx)
        C2 = np.isfinite(yy)
        C = C1 & C2
        if ii == 0:
            x1 = xx[C]
            y1 = yy[C]
        if ii == 1:
            x2 = xx[C]
            y2 = yy[C]
        if ii == 2:
            x3 = xx[C]
            y3 = yy[C]

    if binstyle == 'maxcorr_tau':

        ymax2 = np.max(np.abs(y2))
        if bins1 is None:
            binwidth1 = 0.2
            lim1 = 2 - binwidth1 * 0.5
            bins1 = np.arange(-lim1, lim1 + binwidth1, binwidth1)
        if bins2 is None:
            binwidth2 = ymax2 / 20
            lim2 = 1 + binwidth2
            bins2 = np.arange(limit_corr, lim2 + binwidth2, binwidth2)

        H1, xedge1, yedge1 = np.histogram2d(x1, y1, bins=(bins1, bins2), )
        H1 = np.swapaxes(H1, axis1=0, axis2=1)

        H, xedge, yedge = np.histogram2d(x2, y2, bins=(bins1, bins2))
        H = np.swapaxes(H, axis1=0, axis2=1)
        # ax.hist2d(x, y, bins = 30, cmap=plt.cm.jet)

        xcoord1 = np.zeros(len(xedge1) - 1)
        ycoord1 = np.zeros(len(yedge1) - 1)

        xcoord = np.zeros(len(xedge) - 1)
        ycoord = np.zeros(len(yedge) - 1)

        for xx in range(len(xcoord)):
            xcoord[xx] = 0.5 * (xedge[xx] + xedge[xx + 1])
        for yy in range(len(ycoord)):
            ycoord[yy] = 0.5 * (yedge[yy] + yedge[yy + 1])

        for xx in range(len(xcoord1)):
            xcoord1[xx] = 0.5 * (xedge1[xx] + xedge1[xx + 1])
        for yy in range(len(ycoord1)):
            ycoord1[yy] = 0.5 * (yedge1[yy] + yedge1[yy + 1])
        # ax.plot(x2, y2, '.', markersize = 0.2, color='gray')

        # ax.plot(x2, y2, '.', markersize = 0.2, color='gray')

        # levels1 = np.arange(int(0.13 * np.max(H1)), int(0.7 * np.max(H1)), int(0.15 * np.max(H1)))

        # levels = np.arange(int(0.13 * np.max(H)), int(0.7 * np.max(H)), int(0.15 * np.max(H)))
        levels = [0.2 * np.max(H), 0.4 * np.max(H), 0.6 * np.max(H), 0.8 * np.max(H)]
        if show_background:
            cs = ax.contour(xcoord, ycoord, H, levels=levels,
                            linewidths=1, origin='lower', colors=('#b0aaff', '#8880ff', '#6055ff', '#392bff',))

        # ax.clabel(cs, inline=True, fontsize=10)

        # cs = ax.contour(xcoord1, ycoord1, H1, levels=levels1, linewidths=1, origin='lower', colors='k')
        # ax.clabel(cs, inline=True, fontsize=10)
        h = ax.hist2d(x1, y1, bins=(bins1, bins2), cmap=color, norm=LogNorm(vmin=vmin, vmax=vmax))

    if show_hist_x:
        if show_background & (ax_histx is not None):
            # ax_histx.hist(x2, bins=bins1,histtype='step', color = 'b', label='not cmpf', density = True)
            PlottingFunctions.normalised_histogram(x=x2, ax=ax_histx, bins=bins1, color='b')
        PlottingFunctions.normalised_histogram(x=x1, ax=ax_histx, bins=bins1, color='r')
        PlottingFunctions.normalised_histogram(x=x3, ax=ax_histx, bins=bins1, color='green')

        # ax_histx.hist(x1, bins=bins1, color='r', label='cmpf', histtype='step', density=True)
        # ax_histx.hist(x3, bins=bins1, color='green', label='SLC', histtype='step', density=True)
    if show_hist_y & (ax_histy is not None):
        PlottingFunctions.normalised_histogram(x=y1, ax=ax_histy, bins=bins2, color='r', orientation='horizontal')

        # ax_histy.hist(y1, bins=bins2, orientation='horizontal', color='r', label='cmpf', histtype='step',
        #              density=True)
        if show_background:
            # ax_histy.hist(y2, bins=bins2, orientation='horizontal', color='b', histtype='step', label='not cmpf',
            #          density=True)
            PlottingFunctions.normalised_histogram(x=y2, ax=ax_histy, bins=bins2, color='b', orientation='horizontal')

        # ax_histy.hist(y3, bins=bins2, orientation='horizontal', color='green', histtype='step', label='SLC',
        #              density=True)
        PlottingFunctions.normalised_histogram(x=y3, ax=ax_histy, bins=bins2, color='green', orientation='horizontal')

    # ax_histx.legend(loc='center left', bbox_to_anchor=(1.035, 0.5))
    if (limit_hist1d_x is not None) & (ax_histx is not None):
        ax_histx.set_ylim([0, limit_hist1d_x])
    if (limit_hist1d_y is not None) & (ax_histy is not None):
        ax_histy.set_xlim([0, limit_hist1d_y])
    # ax_histy.legend()
    return bins1, bins2, h
