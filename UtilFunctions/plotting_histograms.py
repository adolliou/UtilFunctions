import numpy as np


class PlottingFunctions:
    @staticmethod
    def normalised_histogram(x, bins, ax, color: str, orientation='vertical', legend=None):
        hist, edges = np.histogram(x, bins=bins)
        hist_normed = hist / hist.sum()
        # centers = (edges[1:] + edges[:-1])/2
        if legend is None:
            ax.stairs(values=hist_normed, edges=edges, color=color, orientation=orientation)
        else:
            ax.stairs(values=hist_normed, edges=edges, color=color, orientation=orientation, label=legend)
