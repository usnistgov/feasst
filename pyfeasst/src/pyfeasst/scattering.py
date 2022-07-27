"""
This module computes quantities to compare with scattering experiments.
"""

import math
import pandas as pd

def structure_factor(distances, radial_distribution, frequency, number_density):
    """
    Return the structure factor at a given frequency by Fourier transform.
    Assume that the distances correspond with the radial distribution and are evenly spaced.
    Each distance is the center of the bin, and the width is twice the first distance.

    :param float number_density: the number of particles divided by the volume.

    >>> from pyfeasst import scattering
    >>> gr = pd.read_csv('../../tests/gr.csv')
    >>> sq = scattering.structure_factor(
    ...     distances=gr['r'],
    ...     radial_distribution=gr['g_HH'],
    ...     frequency=0.4,
    ...     number_density=530/90**3)
    >>> round(sq, 8)
    0.59292662
    """
    struct_fac = 0.
    assert len(distances) == len(radial_distribution)
    bin_width = 2*distances[0]
    for index, distance in enumerate(distances):
        qr = distance*frequency
        struct_fac += bin_width*distance*distance*math.sin(qr)*(radial_distribution[index] - 1)/qr
    return 1 + 4*math.pi*number_density*struct_fac

if __name__ == "__main__":
    import doctest
    doctest.testmod()
