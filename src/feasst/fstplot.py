"""
This module provides some utility plotting functions for use with the FEASST simulation program.
"""

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

def val2map(values, cmname="viridis"):
    """
    Given list of values, return scalar_map for colors.

    Usage:
    color=fstplot.val2map(range)
    for i in range:
        plt.plot(..., color=color.to_rgba(range[i]))
    fstplot.display(range, label='label')
    """
    jet = plt.get_cmap(cmname)
    mn=values[0]
    mx=values[-1]
    if mn > mx:
        mn=values[-1]
        mx=values[0]
    c_norm = colors.Normalize(vmin=mn, vmax=mx)
    return cmx.ScalarMappable(norm=c_norm, cmap=jet)

def display(vals, label='', axes=None, fig=None, cmname="viridis", fontsize=16, location='right', shrink=1):
    """ Display colorbar """
    mn=vals[0]
    mx=vals[-1]
    if mn > mx:
        mn=vals[-1]
        mx=vals[0]
    empty_plot = plt.contourf([[100000, 100000], [100000, 100000]], np.arange(mn, mx + 0.01, 0.01), cmap=cmname)
#    empty_plot = plt.contourf([[0, 0], [0, 0]], np.arange(vals[0], vals[-1] + 0.01, 0.01), cmap=cmname)
    bound = vals
    if axes is None:
        cbar = plt.colorbar(empty_plot, ticks=bound, location=location, shrink=shrink)
    else:
        assert(not fig is None)
        cbar_ax = fig.add_axes(axes)
        cbar = plt.colorbar(empty_plot, ticks=bound, location=location, shrink=shrink, cax=cbar_ax)
    cbar.set_label(label, fontsize=fontsize)
    return cbar

if __name__ == "__main__":
    import doctest
    doctest.testmod()
