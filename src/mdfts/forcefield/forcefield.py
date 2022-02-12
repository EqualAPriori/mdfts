"""
forcefield.py: Constructs a force field for a system. Contains bead types and
potentials between those bead types.
"""
from __future__ import absolute_import, division, print_function

from .potential import _Potential
from .beadtype import BeadType
from mdfts.utils import serial

__all__ = ["ForceField"]


@serial.serialize(["_kT", "_bead_types", "_potentials"])
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
        self._potentials = serial.SerializableTypedDict()

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
        return list(self._bead_types)

    @property
    def bead_names(self):
        """List of bead names of the ForceField"""
        return [bt.name for bt in self.bead_types]

    def add_bead_type(self, bead_type):
        """Add a BeadType to the ForceField"""
        # check that ForceField doesn't already contain a BeadType of the same name
        if bead_type.name in self.bead_names:
            raise ValueError(
                "ForceField instance already contains BeadType '{}'".format(
                    bead_type.name
                )
            )

        self._bead_types.append(bead_type)

    def get_bead_type(self, bead_name):
        """Returns a BeadType of the specified name"""
        # return bead type of specified name if it exists in ForceField
        for bt in self.bead_types:
            if bead_name == bt.name:
                return bt

        # raise error if BeadType with specified name doesn't exist in ForceField
        raise ValueError(
            "ForceField instance does not contain BeadType '{}'".format(bead_name)
        )

    def reorder_bead_types(self, reordered_bead_names):
        """Rearranges the order of BeadTypes in the ForceField. Useful when
        converting to PolyFTS format."""
        if set(reordered_bead_names) != set(self.bead_names):
            raise ValueError(
                "provided bead names don't match bead names in ForceField instance"
            )
        self._bead_types = serial.SerializableTypedList(
            BeadType, *[self.get_bead_type(bn) for bn in reordered_bead_names]
        )

    @property
    def potentials(self):
        """List of _Potentials in the ForceField"""
        _potential_list = []
        for v in self._potentials.values():
            _potential_list.extend(v)
        return _potential_list

    def add_potential(self, potential):
        """Add a _Potential to the ForceField"""
        # check that the potential inherits from base _Potential class
        if not isinstance(potential, _Potential):
            raise TypeError("potential must be of type _Potential or inherit from it")

        # create entry for potential type if it doesn't exist already
        potential_class = potential.__class__
        potential_class_name = potential_class.__name__
        if potential_class_name not in self._potentials.keys():
            self._potentials.add_entry_type(potential_class_name, potential_class, has_many=True)

        # add potential
        # TODO: add warning if _Potential using same BeadTypes already exists in ForceField
        self._potentials[potential_class_name].append(potential)

    def __str__(self):
        """Default method is overridden to allow for more transparent printing
        for debugging purposes"""
        s = "\n"
        s += "BeadTypes: {}\n".format([bt.__str__() for bt in self.bead_types])
        s += "Potentials: {}".format([p.__str__() for p in self.potentials])
        return s
