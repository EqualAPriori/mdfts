""" MD and sim-related system settings """

from __future__ import absolute_import, division, print_function
from logging import warning
import os, warnings

from mdfts.utils import serial
from mdfts.utils import yamlhelper
import mdfts.forcefield as ff
from mdfts.utils import topology
from mdfts.utils import parsify
from mdfts.system import System

__all__ = []


class SystemMD(System):
    def __init__(
        self,
        force_fields=[],
        moldefs=[],
        contents=[],
        chain_aliases={},
        bead_aliases={},
        temp=1.0,
    ):
        super().__init__(
            force_fields, moldefs, contents, chain_aliases, bead_aliases, temp
        )
        pass
        """Other wishlist:
        - box
        - barostat
        - tension
        - integrator
        - neutralize
        - sys_name (?)
        - platform?
        - initial condition
        - runtime
            e.g. equili, write_frequency
        - protocols?
        - Optimizers?
        """
