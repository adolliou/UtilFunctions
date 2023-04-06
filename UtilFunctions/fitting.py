import numpy as np
from matplotlib import pyplot as plt
import itertools


class FittingUtil:

    @staticmethod
    def gaussian(x: np.array, I: np.float64, mu: np.float64, sigma: np.float64, back: np.float64):
        return I * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2)) + back

    @staticmethod
    def multiple_gaussian(x: np.array, I: list, mu: list, sigma: list, back: np.float64):
        s = np.array([I_ * np.exp(-((x - mu_) ** 2) / (2 * sigma_ ** 2)) for I_, mu_, sigma_ in zip(I, mu, sigma)],
                     dtype=np.float64)
        return s.sum(axis=0) + back


class PlotSpectrum:

    @staticmethod
    def _get_edges(x: np.array):

        edges = [x[n] + (x[n + 1] - x[n]) / 2 for n in range(len(x) - 1)]
        edges.append(x[-1] + (x[-1] - x[-2]) / 2)
        edges.insert(0, x[0] - (x[1] - x[0]) / 2)
        return edges

    @staticmethod
    def plot_spectrum(lam: np.array, spectrum: np.array, error_spectrum=None, label="rd", color="r",
                      linewidth_stair=0.5, linewidth_fit=0.5,
                      fig=None, ax=None, show_legend=True, save_path=None, fitting_function=None, kwargs_fitting=None,
                      kwargs_sigma_fitting=None, sigma=3
                      ):

        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot()
        sort = lam.argsort()

        lam = lam[sort]
        spectrum = spectrum[sort]
        isnan = np.isnan(spectrum)
        spectrum = spectrum[~isnan]
        lam = lam[~isnan]
        if error_spectrum is not None:
            error_spectrum = error_spectrum[~isnan]

        edges_lam = PlotSpectrum._get_edges(lam)
        ax.stairs(spectrum, edges=edges_lam, label=label, color=color, linewidth=linewidth_stair, )
        if error_spectrum is not None:
            ax.errorbar(x=lam, y=spectrum, yerr=error_spectrum, linestyle="", marker="", capsize=1,
                        elinewidth=linewidth_stair, color=color)
        if fitting_function is not None:
            if fitting_function == "gaussian":
                fit = FittingUtil.gaussian(lam, **kwargs_fitting)
                if (kwargs_sigma_fitting is not None) & (kwargs_fitting is not None):
                    fits_sigma = PlotSpectrum._extract_all_possible_fits(kwargs_fitting,
                                                            kwargs_sigma_fitting, lam,
                                                            sigma)



            elif fitting_function == "multiple_gaussian":
                fit = FittingUtil.multiple_gaussian(lam, **kwargs_fitting)
                fits_sigma = PlotSpectrum._extract_all_possible_fits(kwargs_fitting,
                                                                     kwargs_sigma_fitting, lam,
                                                                     sigma)
        ax.plot(lam, fit, color=color, linewidth=linewidth_fit, label="_nolegend_")
        if show_legend:
            ax.legend()
        if save_path is not None:
            fig.savefig(save_path)
        return edges_lam, spectrum, lam, fit, fits_sigma

    @staticmethod
    def _extract_all_possible_fits(kwargs_fitting, kwargs_sigma_fitting, lam, sigma):
        fits_sigma = []
        keys = kwargs_sigma_fitting.keys()
        for ii in range(len(keys)):
            for subset in itertools.combinations(keys, ii):
                kwargs_fit_tmp = kwargs_fitting
                for key in subset:
                    kwargs_fit_tmp[key] += 0.5 * sigma * kwargs_sigma_fitting[key]
                    fits_sigma.append(FittingUtil.gaussian(lam, **kwargs_fit_tmp))
                    kwargs_fit_tmp[key] -= 0.5 * sigma * kwargs_sigma_fitting[key]
                    fits_sigma.append(FittingUtil.gaussian(lam, **kwargs_fit_tmp))
        return fits_sigma
