"""
This module accumulates statistical values to analyze its mean, standard deviation and moments.
"""

import numpy as np

class Accumulator:
    """Accumulate values."""
    def __init__(self, num_moments=3, values=None):
        self._moments = np.zeros(num_moments)
        if values is not None:
            self.add(values)

    def num_moments(self):
        """Return the number of moments."""
        return len(self._moments)

    def add(self, value):
        """Add a value or values to the moments."""
        if isinstance(value, list):
            for val in value:
                self.add(val)
        else:
            mom = 1
            for index, _ in enumerate(self._moments):
                self._moments[index] += mom
                mom *= value

    def num_values(self):
        """Return the number of values."""
        return self._moments[0]

    def mean(self):
        """
        Return the mean.

        >>> from pyfeasst import accumulator
        >>> acc = accumulator.Accumulator()
        >>> acc.add(1.)
        >>> acc.add(2.)
        >>> acc.add(4.)
        >>> acc.add([3, 5, 7])
        >>> round(float(acc.mean()), 6)
        3.666667
        """
        mean = float('nan')
        if self.num_values() > 0 and self.num_moments() > 1:
            mean = self.sum_moment(1)/self.num_values()
        return mean

    def stdev(self):
        """
        Return the standard deviation.

        >>> from pyfeasst import accumulator
        >>> acc = accumulator.Accumulator(values=[1, 2, 3, 4, 5, 7])
        >>> round(float(acc.stdev()), 6)
        2.160247
        """
        std = float('nan')
        if self.num_values() > 2 and self.num_moments() > 2:
            num = self.num_values()
            std = np.sqrt((self.sum_moment(2)/num - (self.sum_moment(1)/num)**2)*num/(num-1))
        return std

    def sum_moment(self, power):
        """
        Return the sum of the moment, x**power.

        >>> from pyfeasst import accumulator
        >>> acc = accumulator.Accumulator(values=[1, 2, 3, 4, 5, 7])
        >>> round(float(acc.sum_moment(0)), 6)
        6.0
        >>> round(float(acc.sum_moment(1)), 6)
        22.0
        >>> round(float(acc.sum_moment(2)), 6)
        104.0
        """
        mom = float('nan')
        if self.num_values() > 0 and self.num_moments() > power:
            mom = self._moments[power]
        return mom

if __name__ == "__main__":
    import doctest
    doctest.testmod()
