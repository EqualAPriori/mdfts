"""
potential_harmonic_bond.py: Contains the HarmonicBond class that specifies a
harmonic bond interaction between two specified BeadTypes.
"""
from __future__ import absolute_import, division, print_function

from mdfts.forcefield.potential import _Parameter, _PairPotential

__all__ = ['HarmonicBond']


###############################################################################
# PARAMETERS
###############################################################################


class K(_Parameter):
    """Force constant parameter"""

    def __init__(self, potential, value=1.0, fixed=False):
        """Constructor for K"""
        super(K, self).__init__(potential, float, value, fixed)


class r0(_Parameter):
    """Equilibrium bond length parameter"""

    def __init__(self, potential, value=0.0, fixed=False):
        """Constructor for r0"""
        super(r0, self).__init__(potential, float, value, fixed)


###############################################################################
# HARMONIC BOND POTENTIAL
###############################################################################


class HarmonicBond(_PairPotential):
    """Container object for a harmonic bond potential"""

    _SERIALIZED_PARAMETERS = [K, r0]
