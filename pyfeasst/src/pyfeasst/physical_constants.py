"""
This module provides some physical constants and conversions for convenience.
These constants are taken from CODATA 2018 [https://doi.org/10.1103/RevModPhys.93.025010].
The main purpose of this module is to facilitate comparison of FEASST results to machine percision.
For example, the molar gas constant is derived from the Boltzmann and Avogadro constants.
However, their product to machine percision does not match
scipy.constants.physical_constants["molar gas constant"].
And the conversion of elementary charge per distance to energy units can be even more problematic.
"""

import math

class ElementaryCharge:
    """The elementary charge."""
    def __init__(self):
        self._value = 1.602176634e-19
        self._units = "C"

    def units(self):
        """
        Return the units.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.ElementaryCharge().units()
        'C'
        """
        return self._units

    def value(self):
        """
        Return the value.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.ElementaryCharge().value()
        1.602176634e-19
        """
        return self._value

class BoltzmannConstant:
    """The Boltzmann constant."""
    def __init__(self):
        self._value = 1.380649e-23
        self._units = "J K^-1"

    def units(self):
        """
        Return the units.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.BoltzmannConstant().units()
        'J K^-1'
        """
        return self._units

    def value(self):
        """
        Return the value.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.BoltzmannConstant().value()
        1.380649e-23
        """
        return self._value

class AvogadroConstant:
    """The Avogadro constant."""
    def __init__(self):
        self._value = 6.02214076e+23
        self._units = "mol^-1"

    def units(self):
        """
        Return the units.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.AvogadroConstant().units()
        'mol^-1'
        """
        return self._units

    def value(self):
        """
        Return the value.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.AvogadroConstant().value()
        6.02214076e+23
        """
        return self._value

class VacuumElectricPermittivity:
    """The vacuum electric permittivity."""
    def __init__(self):
        self._value = 8.8541878128e-12
        self._units = "C^2 J^-1 m^-1"

    def units(self):
        """
        Return the units.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.VacuumElectricPermittivity().units()
        'C^2 J^-1 m^-1'
        """
        return self._units

    def value(self):
        """
        Return the value.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.VacuumElectricPermittivity().value()
        8.8541878128e-12
        """
        return self._value

class MolarGasConstant:
    """
    The molar gas constant or ideal gas constant.
    This is equivalent to the product of the Boltzmann and Avogadro constants.
    """
    def __init__(self):
        self._value = BoltzmannConstant().value()*AvogadroConstant().value()
        self._units = "J mol^-1 K^-1"

    def units(self):
        """
        Return the units.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.MolarGasConstant().units()
        'J mol^-1 K^-1'
        """
        return self._units

    def value(self):
        """
        Return the value.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.MolarGasConstant().value()
        8.31446261815324
        """
        return self._value

class ESqPerAngstromToKJPerMol:
    """
    Convert e^2/Angstrom to kJ/mol by a factor of units (kJ*A/e^2/mol).
    This is ElementaryCharge**2/(4*PI*VacuumElectricPermittivity*1e3/1e10/AvogadroConstant).
    """
    def __init__(self):
        self._value = ElementaryCharge().value()**2/ \
            (4*math.pi*VacuumElectricPermittivity().value()*1e3/1e10/AvogadroConstant().value())
        self._units = "kJ A e^-2 mol^-1"

    def units(self):
        """
        Return the units.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.ESqPerAngstromToKJPerMol().units()
        'kJ A e^-2 mol^-1'
        """
        return self._units

    def value(self):
        """
        Return the value.

        >>> from pyfeasst import physical_constants
        >>> physical_constants.ESqPerAngstromToKJPerMol().value()
        1389.35457644382
        """
        return self._value

if __name__ == "__main__":
    import doctest
    doctest.testmod()
