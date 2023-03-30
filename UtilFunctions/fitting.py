import numpy as np
from matplotlib import pyplot as plt


class FittingUtil:

    @staticmethod
    def gaussian(x: np.array, I: np.float64, mu: np.float64, sigma: np.float64, back: np.float64):
        return I * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2)) + back

    @staticmethod
    def multiple_gaussian(x, I: list, mu: list, sigma: list, back: np.float64):
        s = back

        s = np.array([I_ * np.exp(-((x - mu_) ** 2) / (2 * sigma_ ** 2)) for I_, mu_, sigma_ in zip(I, mu, sigma)],
                     dtype=np.float64)
        return s.sum()


class PlotSpectrum:

    @staticmethod
    def _get_edges(x: np.array):

        edges = [x[n + 1] - x[n] for n in range(len(x) - 1)]
        edges.insert(0, x[0] - edges[0])
        edges.append(x[-1] + edges[-1])
        return edges

    @staticmethod
    def plot_spectrum(lam: np.array, spectrum: np.array, error_spectrum=None, label="rd", color="r",
                      linewidth_stair=0.5, linewidth_fit=0.5,
                      fig=None, ax=None, show_legend=True, save_path=None, fitting_function=None, **kwargs_fitting, ):

        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot()
        sort = lam.argsort()
        lam = lam[sort]
        spectrum = spectrum[sort]
        edges_lam = PlotSpectrum._get_edges(lam)
        ax.stairs(spectrum, edges=edges_lam, label=label, color=color, linewidth=linewidth_stair, )
        if error_spectrum is not None:
            ax.errorbar(x=lam, y=spectrum, yerr=error_spectrum, linestyle="", marker="", capsize=1,
                        elinewidth=linewidth_stair, color=color)
        if fitting_function is not None:
            if fitting_function == "gaussian":
                fit = FittingUtil.gaussian(lam, **kwargs_fitting)
            elif fitting_function == "multiple_gaussian":
                fit = FittingUtil.multiple_gaussian(lam, **kwargs_fitting)
        ax.plot(lam, fit, color=color, linewidth=linewidth_fit, label="_nolegend_")
        if show_legend:
            ax.legend()
        if save_path is not None:
            fig.savefig(save_path)
