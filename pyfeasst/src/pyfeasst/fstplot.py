"""
This module provides some utility plotting functions for use with the FEASST simulation program.
"""

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

def val2map(values, cmname="viridis"):
    """ Given list of values, return scalar_map for colors """
    jet = plt.get_cmap(cmname)
    c_norm = colors.Normalize(vmin=values[0], vmax=values[-1])
    return cmx.ScalarMappable(norm=c_norm, cmap=jet)

def display(vals, label='', axes=None, fig=None, cmname="viridis"):
    """ Display colorbar """
    empty_plot = plt.contourf([[100000, 100000], [100000, 100000]], np.arange(vals[0], vals[-1] + 0.01, 0.01), cmap=cmname)
    #empty_plot = plt.contourf([[0, 0], [0, 0]], np.arange(vals[0], vals[-1] + 0.01, 0.01), cmap=cmname)
    bound = vals
    if (axes is None):
        cbar = plt.colorbar(empty_plot, ticks=bound)
    else:
        assert(not fig is None)
        cbar_ax = fig.add_axes(axes)
        cbar = plt.colorbar(empty_plot, ticks=bound, cax=cbar_ax)
    cbar.set_label(label, fontsize=14)
    return cbar

if __name__ == "__main__":
    import doctest
    doctest.testmod()
