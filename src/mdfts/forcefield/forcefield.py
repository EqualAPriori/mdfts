from __future__ import absolute_import, division, print_function

from .potential import _Potential
from mdfts.utils import serial

__all__ = ['BeadType', 'ForceField']


@serial.serialize(['name', 'smear_length', 'charge'])
class BeadType(object):

    def __init__(self, name, smear_length=1.0, charge=0.0):
        self.name = name
        self.smear_length = smear_length
        self.charge = charge

    def __str__(self):
        return self.to_dict()

    def __repr__(self):
        return 'BeadType: {}'.format(self.to_dict())


@serial.serialize(['kT', 'bead_types', 'potentials'])
class ForceField(object):

    def __init__(self, kT=1.0):
        self.kT = kT
        self._bead_types = serial.SerializableTypedList(BeadType)
        self._potentials = serial.SerializableTypedList(_Potential)

    @property
    def kT(self):
        return self._kT

    @kT.setter
    def kT(self, value):
        try:
            self._kT = float(value)
        except ValueError:
            raise ValueError("must be able to convert new value of kT to float")

    @property
    def bead_types(self):
        return self._bead_types

    @bead_types.setter
    def bead_types(self, value):
        self._bead_types = serial.SerializableTypedList(BeadType, value)

    @property
    def bead_names(self):
        return [bt.name for bt in self.bead_types]

    def add_bead_type(self, bead_type):
        if bead_type.name in self.bead_names:
            raise ValueError("ForceField instance already contains BeadType '{}'".format(bead_type.name))
        self._bead_types.append(bead_type)

    def get_bead_type(self, bead_name):
        for bt in self.bead_types:
            if bead_name == bt.name:
                return bt
        raise ValueError("ForceField instance does not contain BeadType '{}'".format(bead_name))

    def reorder_bead_types(self, bead_names):
        if set(bead_names) != set(self.bead_names):
            raise ValueError("provided bead names don't match bead names in ForceField instance")
        self._bead_types = serial.SerializableTypedList(BeadType, *[self.get_bead_type(bn) for bn in bead_names])

    @property
    def potentials(self):
        return self._potentials

    # More transparent printing for debugging. I.e. make sure bead types are getting instantiated
    def __str__(self):
        ret = "\n"
        ret += "BeadTypes: {}\n".format(
            [bt.__str__() for bt in self._bead_types])
        ret += "Potentials: {}".format([p.__str__() for p in self._potentials])
        return ret
