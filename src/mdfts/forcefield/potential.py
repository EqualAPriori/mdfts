"""
potential.py: Contains a base class used to specify interactions between
BeadTypes.
"""
from __future__ import absolute_import, division, print_function

from collections import OrderedDict

from .beadtype import BeadType
from mdfts.utils import serial

__all__ = ['_Parameter', '_Potential']


@serial.serialize(['value', 'fixed'])
class _Parameter(object):
    """Container object for _Potential parameters

    Each parameter in a potential will inherit from this base class.
    """

    def __init__(self, name, parameter_type, value, fixed):
        """Constructor for the _Parameter class"""
        self.name = name
        self._type = parameter_type
        self._value = value
        self._fixed = fixed

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


@serial.serialize(['bead_types', 'parameters'])
class _Potential(object):
    """Container object for a potential

    A _Potential represents an interaction between a specified number of
    BeadTypes. This class is a base class which all other potentials should
    inherit from.
    """
    
    _NUM_BEAD_TYPES = 1
    _PARAMETER_NAMES = []

    def __init__(self, *args, **kwargs):
        """Constructor for the _Potential base class"""
        self.bead_types = args
        self._parameters = serial.SerializableTypedList(_Parameter)

    def __getattr__(self, item):
        if item in self._PARAMETER_NAMES:
            for p in self._parameters:
                if p.name == item:
                    return p
        raise AttributeError("'{}' object has no attribute '{}'".format(self.__class__.__name__, item))

    @property
    def bead_types(self):
        """List of BeadTypes in the Potential"""
        return self._bead_types

    @bead_types.setter
    def bead_types(self, value):
        """Set BeadTypes of the Potential"""
        if len(value) != self._NUM_BEAD_TYPES:
            raise ValueError("attempted to set incorrect number of BeadTypes")
        self._bead_types = serial.SerializableTypedList(BeadType, *list(value))

    def compare_bead_types(self, other_potential):
        """Checks if another Potential shares the same BeadTypes as this
        Potential."""

        # check that other potential is of type _Potential
        if not isinstance(other_potential, _Potential):
            raise TypeError("other potential must be type _Potential or inherit from it")

        # return False if potentials specify different number of bead types
        if self._NUM_BEAD_TYPES != other_potential._NUM_BEAD_TYPES:
            return False

        # return True if bead types of other potential match, False if else
        # TODO: condition only works for _NUM_BEAD_TYPES<=2. Need to generalize
        return set(self.bead_types) == set(other_potential.bead_types)

    @property
    def parameters(self):
        """List of _Parameters in the Potential"""
        return self._parameters

    @parameters.setter
    def parameters(self, value):
        """Set _Parameters of the Potential"""

        # check that correct number of parameters are being set
        if len(value) != len(self._PARAMETER_NAMES):
            raise ValueError("attempted to set incorrect number of Parameters")

        # add parameters to temporary list
        _tmp = [None] * len(self._PARAMETER_NAMES)
        for p in value:
            if p.name not in self._PARAMETER_NAMES:
                raise ValueError("{} is not a valid parameter for the {} potential".format(p.name,
                                                                                           self.__class__.__name__))
            _tmp[self._PARAMETER_NAMES.index(p.name)] = p

        # set parameters equal to that of the temporary variable
        self._parameters = serial.SerializableTypedList(_Parameter, *_tmp)
