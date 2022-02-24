"""
potential_gaussian.py: Contains the Gaussian class that specifies a Gaussian
interaction between two specified BeadTypes.
"""
from __future__ import absolute_import, division, print_function

import math
import warnings

from mdfts.forcefield.potential import _Parameter, _PairPotential

__all__ = ["Gaussian"]


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


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
        return excl_vol / (math.pi / Kappa) ** 1.5

    @value.setter
    def value(self, val):
        """Set the value of the B parameter. The value of the excluded volume
        parameter is adjusted proportionally"""
        Kappa = self.potential.Kappa.value
        excl_vol = val * (math.pi / Kappa) ** 1.5
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
        # Todo: HARD-CODED ASSUMPTION of bead type structure, update
        smear_length_sqd_sum = sum(
            [bt[0].smear_length ** 2 for bt in self.potential.bead_types]
        )
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

    def from_sim_specification(self, sim_spec, kT=1.0):
        super(Gaussian, self).from_sim_specification(sim_spec, kT=kT)
        if sim_spec.bead_names[0] == sim_spec.bead_names[1]:
            self.bead_types[0][0].smear_length = 1.0 / (
                2 * math.sqrt(sim_spec.parameters["Kappa"])
            )
        else:
            if not isclose(self.Kappa.value, sim_spec.parameters["Kappa"]):
                warnings.warn(
                    "Kappas don't match ({} vs {})".format(
                        self.Kappa.value, sim_spec.parameters["Kappa"]
                    )
                )
        self.B.value = sim_spec.parameters["B"] / kT
