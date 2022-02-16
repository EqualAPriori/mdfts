"""
potential.py: Contains base classes used to specify interactions between
BeadTypes and parameters for those interactions.

All parameters will inherit from the _Parameter class while potentials will
inherit from either _Potential or _PairPotential (might add more base classes
later).
"""
from __future__ import absolute_import, division, print_function

from .beadtype import BeadType
from mdfts.utils import serial

__all__ = ['_Parameter', '_Potential', '_PairPotential']


###############################################################################
# PARAMETERS
###############################################################################


@serial.serialize(['_value', '_fixed'])
class _Parameter(object):
    """Container object for _Potential parameters

    Each parameter in a potential will inherit from this base class.
    """

    def __init__(self, potential, parameter_type, value, fixed):
        """Constructor for the _Parameter class"""
        self._potential = potential
        self._type = parameter_type
        self._value = value
        self._fixed = fixed

    @property
    def potential(self):
        """_Potential instance the parameter is associated with"""
        return self._potential

    @property
    def type(self):
        """Type of the parameter"""
        return self._type

    @property
    def value(self):
        """Value of the parameter"""
        return self._value

    @value.setter
    def value(self, val):
        """Set value of the parameter"""
        if not isinstance(val, self.type):
            raise TypeError("value of {} parameter must of type {}".format(self.name, self.type.__name__))
        self._value = val

    @property
    def val(self):
        """Alias for value attribute"""
        return self.value

    @val.setter
    def val(self, value):
        """Alias for value attribute"""
        self.value = value

    @property
    def fixed(self):
        """Flag that indicates whether this parameter is fixed during relative
        entropy minimization"""
        return self._fixed

    @fixed.setter
    def fixed(self, value):
        """Set flag that indicates whether this parameter is fixed during
        relative entropy minimization"""
        if not isinstance(value, bool):
            raise TypeError("'fixed' attribute must be of type bool")
        self._fixed = value


###############################################################################
# POTENTIAL BASE CLASSES
###############################################################################


@serial.serialize(['_bead_types'])
class _Potential(object):
    """Container object for a potential

    A _Potential represents an interaction between a specified number of
    BeadTypes. This class is a base class which all other potentials should
    inherit from.
    """
    
    _NUM_BEAD_TYPES = 1
    _SERIALIZED_PARAMETERS = []
    _NON_SERIALIZED_PARAMETERS = []

    def __init__(self, *args, **kwargs):
        """Constructor for the _Potential base class"""
        # set bead types of the potential
        # TODO: have ability to have multiple bead types
        self.bead_types = args
        self._parameters = serial.SerializableTypedDict()
        for p in self._SERIALIZED_PARAMETERS:
            self._parameters[p.__name__] = p(self)
            self._serial_vars.append(p.__name__)
        for p in self._NON_SERIALIZED_PARAMETERS:
            self._parameters[p.__name__] = p(self)

    def __getattr__(self, item):
        return self._parameters[item]

    @property
    def bead_types(self):
        """List of BeadTypes in the Potential"""
        return list(self._bead_types)

    @bead_types.setter
    def bead_types(self, value):
        """Set BeadTypes of the Potential"""
        if len(value) != self._NUM_BEAD_TYPES:
            raise ValueError("attempted to set incorrect number of BeadTypes")
        self._bead_types = serial.SerializableTypedList(BeadType, *list(value))

    @property
    def bead_names(self):
        """List of names of BeadTypes in the Potential"""
        return [bt.name for bt in self.bead_types]

    def compare_bead_types(self, other_potential):
        """Checks if another Potential shares the same BeadTypes as this
        Potential"""
        # check that other potential is of type _Potential
        if not isinstance(other_potential, _Potential):
            raise TypeError("other potential must be type _Potential or inherit from it")

        # return False if potentials specify different number of bead types
        if self._NUM_BEAD_TYPES != other_potential._NUM_BEAD_TYPES:
            return False

        # return True if bead types of other potential match, False if else
        # TODO: condition only works for _NUM_BEAD_TYPES<=2. Need to generalize
        return set(self.bead_types) == set(other_potential.bead_types)

    def from_sim_specification(self, sim_spec, kT=1.0):
        """Reads in a _SimPotentialSpecification instance and sets values for
        and sets values for the parameters"""
        if sim_spec.bead_names == self.bead_names or sim_spec.bead_names == list(reversed(self.bead_names)):
            pass
        else:
            raise ValueError("bead names in sim specification don't match bead names of the Potential")


class _PairPotential(_Potential):
    """Container object for pair potentials

    The only difference between this class and _Potential is the number of
    BeadTypes"""

    _NUM_BEAD_TYPES = 2
