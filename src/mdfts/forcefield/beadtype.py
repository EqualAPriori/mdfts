"""
beadtype.py: BeadType class that stores information about a specific bead type
in a system.
"""
from __future__ import absolute_import, division, print_function

from mdfts.utils import serial


@serial.serialize(["name", "smear_length", "charge"])
class BeadType(object):
    """Container object for a bead type

    A BeadType stores characteristics about a single bead type in a system,
    such as its name, smear length, and charge. The smear length is in units
    of nanometer while the charge is in terms of the elementary charge
    (1.6e-19 C).
    """

    def __init__(self, name, smear_length=1.0, charge=0.0):
        """Constructor for the BeadType class"""
        self.name = name
        self.smear_length = smear_length
        self.charge = charge

    def __eq__(self, other):
        """Overrides default behavior and compares serialized attributes with
        other BeadType"""

        # check if other object is a BeadType
        if not isinstance(other, BeadType):
            return False

        # compare serialized attributes; return False if any don't match
        for attr in self._serial_vars:
            if self.__getattribute__(attr) != other.__getattribute__(attr):
                return False
        return True

    def __str__(self):
        return self.to_dict()

    def __repr__(self):
        return "BeadType: {}".format(self.to_dict())

    @classmethod
    def init_from_dict(cls, d, *args, **kwargs):
        try:
            obj = cls("dummy")
        except:
            print("Could not initialize {}".format(cls))

        # final update
        obj.from_dict(d)
        return obj


class BeadFilter(serial.SerializableFilterSet):
    """Order-insensitive BeadFilter
    (the potential itself should know whether the filter should be ordered or not).
    """

    def __init__(self, *args, ordered=False):
        super(serial.SerializableFilterSet, self).__init__(*args, ordered=ordered)
        if self._oktype != BeadType and len(self._pattern) > 0:
            raise TypeError("BeadFilter must have contain type `BeadType`")
        self._oktype = BeadType
        self._serial_vars = [
            "_pattern",
            "_ordered",
        ]

    def to_dict(self, shorthand=False):
        res = []
        # res.append(self._ordered)
        for subpattern in self._pattern:
            subres = []
            for sp in subpattern:
                if isinstance(sp, serial.Serializable) and not shorthand:
                    subres.append(sp.to_dict())
                elif isinstance(sp, serial.Serializable) and shorthand:
                    subres.append(sp.name)  # should be a BeadType
                else:
                    subres.append(sp)
            res.append(subres)
        return res

    def from_dict(self, d):
        # self._ordered = d[0]
        res = []
        for subpattern in d[:]:  # d[1:] if including self._ordered
            subres = []
            for sp in subpattern:
                if not isinstance(sp, str):
                    subres.append(BeadType.init_from_dict(sp))
                else:  # treat the entry as just a beadname
                    subres.append(BeadType(sp))
            res.append(tuple(subres))

        self._pattern, self._hashable = serial.process_pattern(res)

    def align_to_dict(self, d):
        """Re-set the bead type parameters in filter to those from a global definition.
        Done by matching bead names.

        Args:
            d (dict): dict of name --> BeadType
        """
        for subpattern in self._pattern:
            for sp in subpattern:
                if sp.name in d:
                    sp.from_dict(d[sp.name].to_dict())

    @property
    def bead_names(self):
        res = []
        for subpattern in self._pattern:
            tmp = []
            for bt in subpattern:
                tmp.append(bt.name)
            res.append(tuple(tmp))
        return tuple(res)  # serial.process_pattern(res)[0]

    """
    @classmethod
    def init_from_dict(cls, d, *args, **kwargs):
        try:
            obj = cls()
        except:
            print("Could not initialize {}".format(cls))

        # final update
        obj.from_dict(d)
        return obj
    """
