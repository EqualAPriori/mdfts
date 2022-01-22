from __future__ import absolute_import, division, print_function

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
        self._bead_types = _BeadTypes()
        self._potentials = _Potentials()

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

    @property
    def bead_names(self):
        return self._bead_types.bead_names()

    def add_bead_type(self, bead_type):
        self._bead_types.add_bead_type(bead_type)

    def get_bead_type(self, bead_name):
        return self._bead_types.get_bead_type(bead_name)

    def reorder_bead_types(self, bead_names):
        if set(bead_names) != set(self.bead_names):
            raise ValueError("provided bead names don't match bead names in ForceField instance")
        self._bead_types = [self.get_bead_type(bn) for bn in bead_names]

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
        # return self.to_dict()

    # Need customized-handling for _bead_types and _potentials to iterate through
    # alternative is to do what `sim` does: define yet another `bead_types(list)`
    # object that extends list, and put the custom get/set in that new object
    # in the future, can further modify these functions to allow support for alternative representations
    def custom_get(self, k):
        if k == '_bead_types':
            return [p.to_dict() for p in self._bead_types]
        elif k == '_potentials':
            return [p.to_dict() for p in self._potentials]
        else:
            return serial.Serializable.custom_get(self, k)

    def custom_set(self, k, v):
        """right now completely over-writes _bead_types and _potentials.
             
        maybe want to raise some kind of error if bead types don't match?
        """
        if k == '_bead_types':
            self._bead_types = []
            for bead_def in v:  # assume `v` is a `list` of `bead_type` dictionaries
                self._bead_types.append(BeadType(**bead_def))
        elif k == '_potentials':
            self._potentials = []
            for potential_def in v:     # assume `v` is a `list` of `` dictionaries
                pass # Implement once potentials are implemented
                # perhaps looks like _POTENTIALTYPES[potential_name](**otherargs)
        else:
            self.custom_set(k, v)
        

@serial.serialize([])
class _BeadTypes(list):

    def __init__(self):
        super(_BeadTypes, self).__init__()

    def add_bead_type(self, bead_type):
        if not isinstance(bead_type, BeadType):
            raise TypeError("bead_type argument must be of type mdfts.forcefield.BeadType")
        if bead_type.name in self.bead_names:
            raise ValueError("BeadType '{}' already exists in ForceField instance".format(bead_type.name))
        self.append(bead_type)

    def get_bead_type(self, bead_name):
        for bt in self:
            if bead_name == bt.name:
                return bt
        raise ValueError("ForceField instance does not contain BeadType '{}'".format(bead_name))

    def bead_names(self):
        return [bt.name for bt in self]


@serial.serialize([])
class _Potentials(list):

    def __init__(self):
        super(_Potentials, self).__init__()
