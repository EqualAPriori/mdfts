"""
potential_harmonic_bond.py: Contains the HarmonicBond class that specifies a
harmonic bond interaction between two specified BeadTypes.
"""
from __future__ import absolute_import, division, print_function

import math
import warnings

from mdfts.forcefield.potential import _Parameter, _PairPotential

__all__ = ['HarmonicBondNoOffset']


###############################################################################
# PARAMETERS
###############################################################################


class b(_Parameter):
    """Root-mean-square bond length parameter (nm)"""

    def __init__(self, potential, value=1.0, fixed=False):
        """Constructor for b"""
        super(b, self).__init__(potential, float, value, fixed)


class K(_Parameter):
    """Force constant parameter (kT/nm^2)"""

    def __init__(self, potential):
        """Constructor for K"""
        super(K, self).__init__(potential, float, None, True)

    @property
    def value(self):
        """Returns the value of the K parameter, which is computed using the
        root-mean-square bond length"""
        return 3. / (2. * self.potential.b.value**2)

    @value.setter
    def value(self, val):
        """Set the value of the K parameter. The value of the root-mean-square
        bond length parameter is adjusted proportionally"""
        self.potential.b.value = math.sqrt(3. / (2. * val))


###############################################################################
# HARMONIC BOND NO OFFSET POTENTIAL
###############################################################################


class HarmonicBondNoOffset(_PairPotential):
    """Container object for a harmonic bond potential"""

    _SERIALIZED_PARAMETERS = [b]
    _NON_SERIALIZED_PARAMETERS = [K]

    def from_sim_specification(self, sim_spec, kT=1.0):
        super(HarmonicBondNoOffset, self).from_sim_specification(sim_spec, kT=kT)
        if sim_spec.parameters["Dist0"] != 0.0:
            warnings.warn("Dist0 parameters is not zero. Ignoring.")
        self.K.value = sim_spec.parameters["FConst"] / kT
