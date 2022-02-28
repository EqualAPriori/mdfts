"""
Very light system object to hold parameters

Also, hold sim optimization parameters
And fts parameters

Handle the exports here?
"""

from __future__ import absolute_import, division, print_function
import os

from mdfts.utils import serial
from mdfts.utils import yamlhelper
import mdfts.forcefield as ff
from mdfts.utils import topology
from mdfts.utils import parsify

__all__ = []


class System(serial.Serializable):
    def __init__(self, force_fields=[], moldefs=[]):
        self._serial_vars = ["paths", "ff_files", "moldef_files"]
        self.paths = ["."]
        self.ff_files = []
        self.update_force_fields(force_fields)
        self.update_moldefs(moldefs)

    def update_force_fields(self, force_fields):
        self.ffs = []
        self.ff_files = []
        if not isinstance(force_fields, (tuple, list)):
            force_fields = [force_fields]
        for ind, ffdef in enumerate(force_fields):
            if isinstance(ffdef, str) and ffdef.endswith((".yml", ".yaml")):
                ffdef_file = ffdef
                try:
                    ff_tmp = ff.ForceField().load(
                        parsify.findpath(ffdef_file, self.paths)
                    )
                except:
                    raise ValueError(
                        "Could not load force field file {}".format(ffdef_file)
                    )
            elif isinstance(ffdef, serial.collectionsABC.Mapping):
                try:
                    ff_tmp = ff.ForceField.init_from_dict(ffdef)
                    ffdef_file = ff_tmp.save("ff_{}.yaml".format(ind))
                    if os.getcwd() not in self._paths:
                        self.paths.append(os.getcwd())
                except:
                    raise ValueError("Could not parse dictionary into force field")
            elif isinstance(ffdef, ff.ForceField):
                ff_tmp = ffdef
                ffdef_file = ff_tmp.save("ff_{}.yaml".format(ind))
                if os.getcwd() not in self._paths:
                    self.paths.append(os.getcwd())
            else:
                raise TypeError("Unrecognized forcefield definition {}".format(ffdef))

            self.ff_files.append(ffdef_file)
            self.ffs.append(ff_tmp)

    def update_moldefs(self, moldefs):
        self.moldefs = []
        self.moldef_files = []
        if not isinstance(moldefs, (tuple, list)):
            moldefs = [moldefs]
        for ind, moldef in enumerate(moldefs):
            if isinstance(moldef, str) and moldef.endswith((".pdb", ".yml", ".yaml")):
                moldef_file = moldef
                try:
                    moldef_tmp = topology.load(
                        parsify.findpath(moldef_file, paths=self.paths)
                    )
                except:
                    raise ValueError(
                        "Could not load moldef file {}".format(moldef_file)
                    )
            elif isinstance(moldef, serial.collectionsABC.Mapping):
                try:
                    moldef_tmp = topology.load(moldef)
                    moldef_file = moldef_tmp.save("moldef_{}.yaml".format(ind))
                    if os.getcwd() not in self._paths:
                        self.paths.append(os.getcwd())
                except:
                    raise ValueError("Could not parse dictionary into moldef")
            elif isinstance(moldef, (topology.FTSTopology, topology.MDTopology)):
                moldef_tmp = moldef
                moldef_file = moldef_tmp.save("moldef_{}.yaml".format(ind))
                if os.getcwd() not in self._paths:
                    self.paths.append(os.getcwd())
            else:
                raise TypeError("Unrecognized moldef definition {}".format(moldef))

            self.moldef_files.append(moldef_file)
            self.moldefs.append(moldef_tmp)
