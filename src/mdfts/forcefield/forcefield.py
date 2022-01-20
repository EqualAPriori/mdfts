from __future__ import absolute_import, division, print_function

from mdfts.utils import serialize   # this doesn't work

__all__ = ['BeadType', 'ForceField']


class BeadType(object):

    def __init__(self, name, smear_length=1.0, charge=0.0):
        self.name = name
        self.smear_length = smear_length
        self.charge = charge
        self._sim_AtomType = None


class ForceField(object):

    def __init__(self, kT=1.0):
        self.kT = kT
        self._bead_types = []
        self._potentials = []

    @property
    def kT(self):
        return self._kT

    @kT.setter
    def kT(self, value):
        try:
            self._kT = float(value)
        except ValueError:
            raise ValueError("must be able to convert new value of kT to float")

    def add_bead_type(self, bead_type):
        if bead_type.name in self.bead_names:
            raise ValueError("ForceField already contains BeadType '{}'".format(bead_type.name))
        self._bead_types.append(bead_type)

    def get_bead_type(self, bead_name):
        for bt in self.bead_types:
            if bead_name == bt.name:
                return bt
        raise ValueError("ForceField does not contain BeadType '{}'".format(bead_name))

    @property
    def bead_types(self):
        return iter(self._bead_types)

    @property
    def bead_names(self):
        for bt in self.bead_types:
            yield bt.name

