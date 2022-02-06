"""
forcefield.py: Constructs a force field for a system. Contains bead types and
potentials between those bead types.
"""
from __future__ import absolute_import, division, print_function

from .potential import _Potential
from .beadtype import BeadType
from mdfts.utils import serial

__all__ = ['ForceField']


@serial.serialize(['kT', 'bead_types', 'potentials'])
class ForceField(object):
    """Container object for a force field

    A ForceField represents the interactions in a system. The ForceField stores
    BeadTypes and Potentials between those bead types. BeadTypes and
    Potentials can be added to the force field through various methods.
    """

    def __init__(self, kT=1.0):
        """Constructor for the ForceField class"""
        self.kT = kT
        self._bead_types = serial.SerializableTypedList(BeadType)
        self._potentials = serial.SerializableTypedList(_Potential)

    @property
    def kT(self):
        """Thermal energy of the system. Used when outputting parameters of
        potentials in terms of real energy units."""
        return self._kT

    @kT.setter
    def kT(self, value):
        """Set the thermal energy of the system"""
        try:
            self._kT = float(value)
        except ValueError:
            raise ValueError("must be able to convert new value of kT to float")

    @property
    def bead_types(self):
        """List of BeadTypes in the ForceField"""
        return self._bead_types

    @bead_types.setter
    def bead_types(self, value):
        """Set BeadTypes of the ForceField"""
        self._bead_types = serial.SerializableTypedList(BeadType, value)

    @property
    def bead_names(self):
        """List of bead names of the ForceField"""
        return [bt.name for bt in self.bead_types]

    def add_bead_type(self, bead_type):
        """Add a BeadType to the ForceField"""

        # check that ForceField doesn't already contain a BeadType of the same name
        if bead_type.name in self.bead_names:
            raise ValueError("ForceField instance already contains BeadType '{}'".format(bead_type.name))

        self._bead_types.append(bead_type)

    def get_bead_type(self, bead_name):
        """Returns a BeadType of the specified name"""

        # return bead type of specified name if it exists in ForceField
        for bt in self.bead_types:
            if bead_name == bt.name:
                return bt

        # raise error if BeadType with specified name doesn't exist in ForceField
        raise ValueError("ForceField instance does not contain BeadType '{}'".format(bead_name))

    def reorder_bead_types(self, bead_names):
        """Rearranges the order of BeadTypes in the ForceField. Useful when
        converting to PolyFTS format """
        if set(bead_names) != set(self.bead_names):
            raise ValueError("provided bead names don't match bead names in ForceField instance")
        self._bead_types = serial.SerializableTypedList(BeadType, *[self.get_bead_type(bn) for bn in bead_names])

    @property
    def potentials(self):
        """List of Potentials in the ForceField"""
        return self._potentials

    def add_potential(self, potential):
        """Add a BeadType to the ForceField"""
        # TODO: add warning if Potential using same BeadTypes already exists in ForceField
        self._potentials.append(potential)

    def __str__(self):
        """Default method is overridden to allow for more transparent printing
        for debugging purposes"""
        ret = "\n"
        ret += "BeadTypes: {}\n".format(
            [bt.__str__() for bt in self._bead_types])
        ret += "Potentials: {}".format([p.__str__() for p in self._potentials])
        return ret
