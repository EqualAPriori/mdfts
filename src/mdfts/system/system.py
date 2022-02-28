"""
Very light system object to hold parameters

Also, hold sim optimization parameters
And fts parameters

Handle the exports here?
"""

from __future__ import absolute_import, division, print_function
from logging import warning
import os, warnings

from mdfts.utils import serial
from mdfts.utils import yamlhelper
import mdfts.forcefield as ff
from mdfts.utils import topology
from mdfts.utils import parsify

__all__ = []


class System(serial.Serializable):
    def __init__(
        self,
        force_fields=[],
        moldefs=[],
        contents=[],
        chain_aliases={},
        bead_aliases={},
    ):
        """Basic system container

        Args:
            force_fields (list, optional): _description_. Defaults to [].
            moldefs (list, optional): _description_. Defaults to [].
            contents (list or tuple, optional): [(mol_name,mol_num)] pairs. Defaults to None.
            chain_aliases (dict, optional): in the form of {preferred_name:actual_name_in_def}. Defaults to None.
            bead_aliases (dict, optional): in the form of {preferred_name:actual_name_in_def}. Defaults to None.
        """
        self._serial_vars = [
            "paths",
            "ff_files",
            "moldef_files",
            "contents",
            "chain_aliases",
            "bead_aliases",
        ]
        self.paths = ["."]
        self.ff_files = []
        self.update_force_fields(force_fields)
        self.update_moldefs(moldefs)
        self.update_contents(contents)
        self.chain_aliases = chain_aliases
        self.bead_aliases = bead_aliases

    def add_force_field(self, ffdef):
        if isinstance(ffdef, str) and ffdef.endswith((".yml", ".yaml")):
            ffdef_file = ffdef
            try:
                ff_tmp = ff.ForceField().load(parsify.findpath(ffdef_file, self.paths))
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

    def update_force_fields(self, force_fields):
        self.ffs = []
        self.ff_files = []
        if not isinstance(force_fields, (tuple, list)):
            force_fields = [force_fields]
        for ind, ffdef in enumerate(force_fields):
            self.add_force_field(ffdef)

    def add_moldef(self, moldef):
        if isinstance(moldef, str) and moldef.endswith((".pdb", ".yml", ".yaml")):
            moldef_file = moldef
            try:
                moldef_tmp = topology.load(
                    parsify.findpath(moldef_file, paths=self.paths)
                )
            except:
                raise ValueError("Could not load moldef file {}".format(moldef_file))
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

    def update_moldefs(self, moldefs):
        self.moldefs = []
        self.moldef_files = []
        if not isinstance(moldefs, (tuple, list)):
            moldefs = [moldefs]
        for ind, moldef in enumerate(moldefs):
            self.add_moldef(moldef)

    def add_molecule(mol_name, mol_num):
        mol_found = False
        for moldef in self.moldefs:
            if isinstance(moldef, topology.MDTopology):
                if mol_name in moldef:
                    mol_found = True
                    break
            elif isinstance(moldef, topology.FTSTopology):
                if mol_name in moldef.chain_types:
                    mol_found = True
                    break
        if mol_found:
            self.contents.append(mol_name, mol_num)
        else:
            warnings.warn("molecule {} not defined".format(mol_name))

    def update_contents(self, contents):
        self.contents = []
        for mol_name, mol_num in contents:
            self.add_molecule(mol_name, mol_num)

    # ===== Serializable Interface =====
    def save(self, filename="system.yaml"):
        yamlhelper.save_dict(filename, self.to_dict())

    def load(self, filename="system.yaml"):
        sysdef = yamlhelper.load(filename)
        self.from_dict(sysdef)


if __name__ == "__main__":
    s = System(moldefs=["../utils/moldef.yaml", "../utils/test_Ch.pdb"])
