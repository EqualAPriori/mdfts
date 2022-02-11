"""
potential_gaussian.py: Contains the Gaussian class that specifies a Gaussian
interaction between two specified BeadTypes.
"""
from __future__ import absolute_import, division, print_function

import math

from mdfts.forcefield.potential import _Parameter, _PairPotential

__all__ = ['Gaussian']


###############################################################################
# PARAMETERS
###############################################################################


class excl_vol(_Parameter):
    """Excluded volume parameter (nm^3)"""

    def __init__(self, potential, value=0.0, fixed=False):
        """Constructor for the excluded volume parameter"""
        super(excl_vol, self).__init__(potential, float, value, fixed)


class B(_Parameter):
    """Exponential prefactor parameter

    This parameter is directly proportional to the excluded volume.
    """

    def __init__(self, potential):
        """Constructor for the exponential prefactor parameter"""
        super(B, self).__init__(potential, float, None, False)

    @property
    def value(self):
        """Returns the value of the B parameter, which is computed using the
        values of the excluded volume and Kappa parameters"""
        excl_vol = self.potential.excl_vol.value
        Kappa = self.potential.Kappa.value
        return excl_vol * (Kappa / math.pi)**1.5

    @value.setter
    def value(self, val):
        """Set the value of the B parameter. The value of the excluded volume
        parameter is adjusted proportionally"""
        Kappa = self.potential.Kappa.value
        excl_vol = val / (Kappa / math.pi)**1.5
        self.potential.excl_vol.value = excl_vol


class Kappa(_Parameter):
    """Length scale parameter (nm^-2)

    This parameter is used in both the exponent and the prefactor."""

    def __init__(self, potential):
        """Constructor for the Kappa parameter"""
        super(Kappa, self).__init__(potential, float, None, True)

    @property
    def value(self):
        """Value of the Kappa parameter, which is computed using the smearing
        lengths of the BeadTypes in the potential"""
        smear_length_sqd_sum = sum([bt.smear_length**2 for bt in self.potential.bead_types])
        return 0.5 / smear_length_sqd_sum

    @value.setter
    def value(self, val):
        """The value of the Kappa parameter can't be set, as it depends on the
        smearing lengths of the BeadTypes in the potential"""
        raise NotImplementedError("cannot set value of Kappa parameter")


###############################################################################
# GAUSSIAN POTENTIAL
###############################################################################


class Gaussian(_PairPotential):
    """Container object for a Gaussian potential"""

    _SERIALIZED_PARAMETERS = [excl_vol]
    _NON_SERIALIZED_PARAMETERS = [B, Kappa]
