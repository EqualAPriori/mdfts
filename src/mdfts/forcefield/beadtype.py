"""
beadtype.py: BeadType class that stores information about a specific bead type
in a system.
"""
from __future__ import absolute_import, division, print_function

from mdfts.utils import serial


@serial.serialize(['name', 'smear_length', 'charge'])
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

    def __str__(self):
        return self.to_dict()

    def __repr__(self):
        return 'BeadType: {}'.format(self.to_dict())
