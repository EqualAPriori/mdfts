""" Helper functions for representing topolgies

Todo:
    1) aliases
    
    2) SMILES translator. 
    
        Think about where I'd want to use that 
        --> in residue definition? molecule definition?
        essentially, use SMILES to generate the beadlist and bondlist
    
    3) smarter pathing for loading files
"""
# standard imports
import os, copy, re
from collections import OrderedDict

# 3rd party imports
import numpy as np
import mdtraj

# local imports
from . import yamlhelper as yaml
from . import parsify

version = "list"

# ===== General Utilities =====
def create_system(tops):
    """
    tops: list
        of tuple pairs of (mdtraj topology, #replicates)

    Notes:
    ------
    in future, consider creating topology from string/list specification
    also can handle of tops are trajectories... then spit out a trajectory
    """
    n_entries = len(tops)

    sub_topologies = [None] * n_entries
    new_topology = mdtraj.Topology()
    for ic, chain in enumerate(tops):
        tmp_top = chain[0]
        if isinstance(chain[0], mdtraj.Trajectory):
            tmp_top = chain[0].topology

        sub_topologies[ic] = replicate_topology(tmp_top, chain[1])
        new_topology = new_topology.join(sub_topologies[ic])

    if all([isinstance(entry[0], mdtraj.Trajectory) for entry in tops]):
        xyzlist = []
        # then also replicate coordinates
        for ic, chain in enumerate(tops):
            xyz = chain[0].xyz[0]
            n_replicates = chain[1]
            xyzlist.append(np.tile(xyz, [1, n_replicates, 1]))
        xyzs = np.concatenate(xyzlist, axis=1)
        new_topology = mdtraj.Trajectory(xyzs, new_topology)

    return new_topology


def replicate_topology(top, n_replicates):
    """
    Try to be more efficient than sequentially adding (which I believe can cost a lot of iterations)
    """
    bits = np.binary_repr(n_replicates)
    power2_replicates = [top]
    for ia in range(len(bits) - 1):
        power2_replicates.append(power2_replicates[-1].join(power2_replicates[-1]))

    new_top = mdtraj.Topology()
    for ib, bit in enumerate(bits):
        if int(bit) == 1:
            new_top = new_top.join(power2_replicates[-(ib + 1)])

    return new_top


def symbolify(index):
    """
    Notes
    -----
    return two-character symbol.
    For index < 100, return 0-99. 
    For index > 100, return 0A, 0B, etc.
    This scheme can only handle 360 bead types...
    I think symbols are case insensitive
    """
    if index < 100:
        symbol = "{:02d}".format(index)
    else:
        prefix = str(int((index - 100) / 26))
        suffix = chr(ord("A") + (index - 100) % 26)
        symbol = prefix + suffix

    return symbol


cg_elements = [
    (e.name, e)
    for n, e in mdtraj.element.Element._elements_by_atomic_number.items()
    if n >= 200
]
cg_beadtypes = OrderedDict(cg_elements)


def add_cg_beadtype(name, mass, dummy_radius=1.0):
    """
    convention is to use numbers > 200 for my atomic numbers.
               & preprend symbols with digits to avoid conflict w/ existing atoms.
    """
    if name not in cg_beadtypes:
        print("... successfully added new element {} with mass {}".format(name, mass))
        index = len(cg_beadtypes)
        symbol = symbolify(index)
        cg_beadtypes[name] = mdtraj.element.Element(
            200 + index, name, symbol, mass, dummy_radius
        )
    else:
        print("... found existing element {} with mass {}".format(name, mass))
    return cg_beadtypes[name]


def generate_bond_list(args, style="simple"):
    if style.lower() in ["simple", "linear"]:
        "interpret args as # beads to bond"
        nbeads = args
        new_bond_list = []
        for ia in range(nbeads - 1):
            new_bond_list += [(ia, ia + 1)]
    else:
        raise ValueError("style {} not supported".format(style))

    return new_bond_list


def flatten_shorthand(mylist):
    """
    flatten shorthand of type [ (name0,#0), name1, name1, (name2,#2), ... ]
    """
    new_list = []
    for entry in mylist:
        if isinstance(entry, tuple) or isinstance(entry, list):
            if len(entry) == 2:
                name = entry[0]
                num_thing = entry[1]
                new_list += [name] * num_thing
            else:
                raise ValueError(
                    "entry {} not interpretable. Should be string or 2-element list/tuple"
                )
        elif isinstance(entry, str):
            new_list += [entry]
        else:
            raise ValueError(
                "atom entry {} is not interpretable. Should be string or 2-element list/tuple"
            )
    return new_list


def load(trajfile, top=None):
    """ To replace mdtraj.load(trajfile, topfile) shorthand, with .yaml or csv+bond dat files"""
    if isinstance(top, str):
        if top.endswith(".yaml"):
            print("detected .yaml topology file")
            new_top = Topology(top)
            processed_top = new_top.system.topology
        elif top.endswith(".pdb"):
            processed_top = top
        else:
            raise ValueError("Can not parse given top {}".format(top))
    elif isinstance(top, list) or isinstance(top, tuple):
        # if using .csv, bond.dat pair to define topology
        if len(top) != 2:
            raise ValueError(
                "Seems like given topology should be .csv, bond.dat pair, but did not get exactly two inputs"
            )
        df = pandas.read_csv(top[0])
        bonds = np.loadtxt(top[1])
        processed_top = mdtraj.Topology.from_dataframe(df, bonds)
    elif isinstance(top, Topology):
        processed_top = (
            top.system.topology
        )  # Topology object should be already processed/initialized, top.system is actually a mdtraj Trajectory
    elif isinstance(top, mdtraj.Topology):
        processed_top = top

    if top is None:
        traj = mdtraj.load(trajfile)
    else:
        traj = mdtraj.load(trajfile, top=processed_top)

    return traj


# ===== handling files =====
class Topology:
    """
    Parameters:
    -----------
    topdef: string or dict

    Notes:
    ------
    simplest specification for now:
        *) bead_types
        *) mol_types
        *) system

    In the future, build in:
        *) defining molecules in line, e.g. using residues, lists

    Todo:
    - allow for pickle, csv+bond.dat definitions
    """

    def __init__(self, topdef=None):
        self.system = None  # mdtraj trajectory object of entire system
        # objects that exporter will use, use mdtraj objects where possible:
        self.bead_types = {}  # dictionary to mdtraj elements
        self.res_types = {}  # dictionary of dictionary definitions for residues
        self.mol_types = {}  # dictionary to trajectories of each chain
        self.system_def = []  # list of (trajectory, #) pairs

        # objects having to do with yaml representation:
        self.loaded_file = None
        self.processed_file = None
        """
        self.definition = {
                'paths': './',
                'bead_types': [],
                'res_types': [],
                'mol_types': [],
                'system': None
                }
        """

        if topdef is not None:
            self.load(topdef)

    def process_system_pdb(self, filename):
        print("===== Assuming loading a .pdb topology file of the entire system =====")
        dummy_mass = 1.0
        dummy_radius = 1.0
        dummy_charge = 0.0

        self.system = mdtraj.load(filename)
        self.processed_file = {"bead_types": [], "system": filename}

        # assume every unique atom name is an atom bead type
        for atom in self.system.topology.atoms:
            if atom.name not in self.bead_types:

                # create dummy elements/bead types
                print(
                    "Assuming default mass={}, charge={} for bead {}".format(
                        dummy_mass, dummy_charge, atom.name
                    )
                )
                self.bead_types[atom.name] = add_cg_beadtype(
                    name, dummy_mass, dummy_radius=dummy_radius
                )
                # self.bead_types.append( mdtraj.element.Element(200+ib, name, symbol, dummy_mass, dummy_radius) )
                self.processed_file["bead_types"].append(
                    {"name": atom.name, "mass": dummy_mass, "charge": dummy_charge}
                )

            # reassign atom types
            atom.element = self.bead_types[atom.name]

    def process_system_dict(self):
        """
        Notes
        -----
        acts on stored self.loaded_file and self.processed_file
        system object should be list of (moltype, #mol) tuples/lists.
        if system topology defined by a full-system pdb, then might as well skip this whole process.
        """
        print("===== Building topology =====")
        # now need to parse the sections...
        if "paths" not in self.loaded_file:
            self.processed_file["paths"] = []

        ## bead types, mandatory
        print("\n=== Processing Bead Types =====")
        if "bead_types" not in self.loaded_file:
            print("no bead definitions found, skipping")
        elif self.loaded_file["bead_types"] is None:
            print("no bead definitions found, skipping")
        elif isinstance(self.loaded_file["bead_types"], list):
            num_bead_types = 0
            bead_type_fields = None
            bead_type_defaults = None
            bead_type_field_index = None
            custom_fields = False
            for ind_bead_def, bead_def in enumerate(self.loaded_file["bead_types"]):
                print("--->processing entry {}: {}".format(ind_bead_def, bead_def))
                if (
                    ind_bead_def == 0
                    and isinstance(bead_def, dict)
                    and "fields" in bead_def
                ):
                    bead_type_fields = bead_def["fields"]
                    bead_type_defaults = bead_def["defaults"]
                    bead_type_field_index = {
                        n: i for i, n in enumerate(bead_type_fields)
                    }
                    custom_fields = True
                if isinstance(bead_def, str):
                    tmp_bead_def = bead_def.split()
                    bead_def = []
                    for ie, entry in enumerate(tmp_bead_def):
                        if ie == 0:
                            bead_def.append(entry)
                        else:
                            bead_def.append(float(entry))
                if isinstance(bead_def, list):
                    tmp_bead_def = {}
                    for ie, entry in enumerate(bead_def):
                        if ie == 0:
                            tmp_bead_def["name"] = entry
                        else:
                            if custom_fields:
                                tmp_bead_def[bead_type_fields[ie - 1]] = entry
                            elif ie == 1:
                                tmp_bead_def["mass"] = entry
                            elif ie == 2:
                                tmp_bead_def["charge"] = entry
                    bead_def = tmp_bead_def
                if isinstance(bead_def, dict) and "fields" not in bead_def:
                    bead_name = bead_def["name"]
                    index = len(self.bead_types)
                    if "mass" in bead_def:
                        mass = bead_def["mass"]
                    elif custom_fields:
                        mass = bead_type_defaults[bead_type_field_index["mass"]]
                    else:
                        raise ValueError("mass not defined")

                    self.bead_types[bead_name] = add_cg_beadtype(bead_name, mass)
                    if custom_fields:
                        processed_bead_def = {}
                        processed_bead_def["name"] = bead_def["name"]
                        for ii, field in enumerate(bead_type_fields):
                            processed_bead_def[field] = bead_type_defaults[ii]
                        for field, value in bead_def.items():
                            processed_bead_def[field] = value

                        self.processed_file["bead_types"][
                            ind_bead_def
                        ] = processed_bead_def
                        print(
                            "\tfinal processed bead def: {}".format(processed_bead_def)
                        )
                    else:
                        print("\tfinal processed bead def: {}".format(bead_def))
        else:
            raise ValueError("dictionary input for bead types not supported yet")

        ## residue types
        print("\n=== Processing Residue Types =====")
        if "bead_types" not in self.loaded_file:
            print("no residue definitions found, skipping")
        elif self.loaded_file["res_types"] is None:
            print("no residue definitions found, skipping")
        elif isinstance(self.loaded_file["res_types"], list):
            for ir, res_def in enumerate(self.loaded_file["res_types"]):
                resname = res_def["name"]
                print("--->Working on residue {}: {}".format(ir, resname))

                # linearize the atom list if needed, for easier bonding enumeration
                old_bead_def = res_def["beads"]
                new_bead_list = flatten_shorthand(old_bead_def)
                self.processed_file["res_types"][ir]["beads"] = new_bead_list

                # now add bonds, intra-residue
                if "bonds" in res_def:
                    bond_def = res_def["bonds"]
                else:
                    bond_def = "simple"
                if isinstance(bond_def, str):
                    if bond_def.lower() in ["simple", "contiguous", "linear"]:
                        new_bond_list = generate_bond_list(
                            len(new_bead_list), style="simple"
                        )
                        self.processed_file["res_types"][ir]["bonds"] = new_bond_list

                # now process linker beads defaults
                if "head" not in res_def:
                    print("defaulting head link for {} to be bead 0".format(resname))
                    self.processed_file["res_types"][ir]["head"] = 0
                if "tail" not in res_def:
                    print("defaulting tail link for {} to be bead 0".format(resname))
                    self.processed_file["res_types"][ir]["tail"] = 0

                self.res_types[resname] = self.processed_file["res_types"][ir]
        else:
            raise ValueError("dictionary input for bead types not supported yet")

        ## moltypes, mandatory
        print("\n=== Processing Molecule/Chain Types =====")
        for im, mol_entry in enumerate(self.loaded_file["mol_types"]):
            print(mol_entry)
            print("--->Working on mol {}: {}".format(im, mol_entry["name"]))
            mol_name = mol_entry["name"]
            if "def" in mol_entry:
                print("\tTrying to read in {}".format(mol_entry["def"]))
                if isinstance(mol_entry["def"], str):  # a molecule definition file
                    mol_file = mol_entry["def"]
                    if mol_file.endswith("pdb"):
                        full_file = parsify.findpath(
                            mol_file, self.processed_file["paths"]
                        )
                        new_chain_top, new_chain_def = self.add_mol_type(
                            mol_name, full_file, style="pdb"
                        )
                        self.processed_file["mol_types"][im] = new_chain_def
                    else:
                        raise ValueError(
                            "unrecognized molecule definition file {}".format(
                                mol_entry["def"]
                            )
                        )
                if isinstance(mol_entry["def"], list):  # should be list of residues
                    print("\tattempting to build molecule from residues")
                    new_chain_top, new_chain_def = self.add_mol_type(
                        mol_name,
                        definition=mol_entry["def"],
                        style="residues",
                        bonds="simple",
                    )
                    self.processed_file["mol_types"][im] = new_chain_def
            else:  # build up from beads
                print("\tattempting to build molecule from beads")
                new_chain_top, new_chain_def = self.add_mol_type(
                    mol_name,
                    definition=mol_entry["beads"],
                    style="beads",
                    bonds="simple",
                )
                self.processed_file["mol_types"][im] = new_chain_def

        ## system set up
        print("\n=== Processing System =====")
        if isinstance(self.loaded_file["system"], str):
            if self.loaded_file["system"].endswith("pdb"):
                full_file = parsify.findpath(
                    self.loaded_file["system"], self.processed_file["paths"]
                )
                self.system = mdtraj.load(full_file)
            else:
                raise ValueError(
                    "uncrecognized system definition file: {}".format(
                        self.loaded_file["system"]
                    )
                )
        else:
            for im, entry in enumerate(self.loaded_file["system"]):
                print("--->Working on entry {}: {}".format(im, entry))
                if isinstance(entry, str):
                    # should be in format: name number
                    # ca use whitespace, `:`, `;` to separate
                    entry = re.split(r"[\s,:;]+", entry)
                    if len(entry) != 2:
                        raise ValueError(
                            "entry {} does not have 2 elements as required".format(
                                entry
                            )
                        )
                    entry = [entry[0], int(entry[1])]
                    self.processed_file["system"][im] = entry

                _mol_type = self.mol_types[entry[0]]
                _n_replicates = int(entry[1])
                self.system_def.append((_mol_type, _n_replicates))
            print("replicating chains and assembling system")
            self.system = create_system(self.system_def)

    def add_mol_type(self, mol_name, definition, style, bonds="simple"):
        """
        use internally stored bead_type, res_type, mol_type data structures
        """
        print("adding molecule {} with definition style {}".format(mol_name, style))
        if style.lower() in ["pdb"]:
            mol_file = definition
            tmp_traj = mdtraj.load(mol_file)
            tmp_top = tmp_traj.topology
            atoms_in_mol = [a.name for a in tmp_top.atoms]
            bonds_in_mol = [(b[0].index, b[1].index) for b in tmp_top.bonds]
            if len(bonds_in_mol) == 0:
                if isinstance(bonds, list):
                    new_bond_list = bonds
                elif bonds.lower in ["simple", "linear"]:
                    print(
                        "bonds not detected in {}, filling in using style {}".format(
                            mol_file, bonds
                        )
                    )
                    new_bond_list = generate_bond_list(
                        len(new_bead_list), style="simple"
                    )
                else:
                    new_bond_list = []
                    print(
                        "caution! no bonds detected in {}, and no supplemental bond definition given"
                    )
                for bond_pair in new_bond_list:
                    a0 = tmp_top.atom(bond_pair[0])
                    a1 = tmp_top.atom(bond_pair[1])
                    tmp_top.add_bond(a0, a1)
                bonds_in_mol = [(b[0].index, b[1].index) for b in tmp_top.bonds]

            for atom in tmp_traj.topology.atoms:
                if atom.name not in self.bead_types:
                    bead_name = atom.name
                    print(
                        "Need to add new cg bead type {}, default mass {}".format(
                            bead_name, 1.0
                        )
                    )
                    self.bead_types[atom.name] = add_cg_beadtype(bead_name, 1.0)
                atom.element = self.bead_types[atom.name]
            self.mol_types[mol_name] = tmp_traj
            mol_def = {
                "name": mol_name,
                "def": mol_file,
                "beads": atoms_in_mol,
                "bonds": bonds_in_mol,
            }

        elif style.lower() in ["residues"]:
            shorthand_res_list = definition
            flat_res_list = flatten_shorthand(shorthand_res_list)

            tmp_top = mdtraj.Topology()
            ch = tmp_top.add_chain()
            for ir, res_name in enumerate(flat_res_list):
                r = tmp_top.add_residue(res_name, ch)
                res_def = self.res_types[res_name]
                atoms_in_res = []
                for bead_name in res_def["beads"]:
                    if bead_name not in self.bead_types:  # NEED TO ADD ATOM
                        print(
                            "Need to add new cg bead type {}, default mass {}".format(
                                bead_name, 1.0
                            )
                        )
                        self.bead_types[bead_name] = add_cg_beadtype(bead_name, 1.0)
                        self.processed_file["bead_types"].append(
                            {"name": bead_name, "mass": 1.0, "charge": 0.0}
                        )
                    new_atom = tmp_top.add_atom(
                        bead_name, self.bead_types[bead_name], r
                    )

                    atoms_in_res.append(new_atom)

                if isinstance(bonds, str):
                    if bonds.lower() in ["simple", "linear"]:
                        if ir > 0:  # add bonds between residues
                            tmp_top.add_bond(prev_tail, atoms_in_res[res_def["head"]])
                        prev_tail = atoms_in_res[res_def["tail"]]
                    else:
                        raise ValueError(
                            "unknown mol bond specification: {}".format(bonds)
                        )
                elif isinstance(bonds, list):
                    print(
                        "interpreting additional bonding specification as further intra-chain-indexed bonds to add"
                    )
                    for bond_pair in bonds:
                        a0 = tmp_top.atom(bond_pair[0])
                        a1 = tmp_top.atom(bond_pair[1])
                        tmp_top.add_bond(a0, a1)
                else:
                    raise ValueError("unknown mol bond specification: {}".format(bonds))
                for bond_pair in res_def["bonds"]:
                    tmp_top.add_bond(
                        atoms_in_res[bond_pair[0]], atoms_in_res[bond_pair[1]]
                    )

            atoms_in_mol = [a.name for a in tmp_top.atoms]
            bonds_in_mol = [(b[0].index, b[1].index) for b in tmp_top.bonds]
            for atom in tmp_top.atoms:
                atom.element = self.bead_types[atom.name]
            self.mol_types[mol_name] = mdtraj.Trajectory(
                np.zeros([1, len(atoms_in_mol), 3]), tmp_top
            )
            mol_def = {
                "name": mol_name,
                "def": shorthand_res_list,
                "beads": atoms_in_mol,
                "bonds": bonds_in_mol,
            }

        elif style.lower() in ["beads"]:
            # building up from bead names, skipping residue definition!
            shorthand_bead_list = definition
            flat_bead_list = flatten_shorthand(shorthand_bead_list)

            tmp_top = mdtraj.Topology()
            ch = tmp_top.add_chain()
            for ib, bead_name in enumerate(flat_bead_list):
                r = tmp_top.add_residue(bead_name, ch)
                if bead_name not in self.bead_types:  # NEED TO ADD ATOM
                    print(
                        "Need to add new cg bead type {}, default mass {}".format(
                            bead_name, 1.0
                        )
                    )
                    self.bead_types[bead_name] = add_cg_beadtype(bead_name, 1.0)
                    self.processed_file["bead_types"].append(
                        {"name": bead_name, "mass": 1.0, "charge": 0.0}
                    )
                new_atom = tmp_top.add_atom(bead_name, self.bead_types[bead_name], r)

            if isinstance(bonds, str):
                if bonds.lower() in ["simple", "linear"]:
                    new_bond_list = generate_bond_list(tmp_top.n_atoms, style="simple")
            elif isinstance(bonds, list):
                new_bond_list = bonds
            else:
                raise ValueError("unknown mol bond specification: {}".format(bonds))
            for bond_pair in new_bond_list:
                a0 = tmp_top.atom(bond_pair[0])
                a1 = tmp_top.atom(bond_pair[1])
                tmp_top.add_bond(a0, a1)

            atoms_in_mol = [a.name for a in tmp_top.atoms]
            bonds_in_mol = [(b[0].index, b[1].index) for b in tmp_top.bonds]
            for atom in tmp_top.atoms:
                atom.element = self.bead_types[atom.name]
            self.mol_types[mol_name] = mdtraj.Trajectory(
                np.zeros([1, len(atoms_in_mol), 3]), tmp_top
            )
            mol_def = {"name": mol_name, "beads": atoms_in_mol, "bonds": bonds_in_mol}

        return tmp_top, mol_def

    def load(self, topdef):
        if isinstance(topdef, str):
            if topdef.endswith(
                "pdb"
            ):  # assume atom types, bonds, etc. defined via pdb for the entire system
                self.process_system_pdb(topdef)
            else:  # assume topdef is a .yaml file
                self.loaded_file = yaml.load(topdef)
                self.processed_file = yaml.load(topdef)
                self.process_system_dict()
        else:
            # Allow for simplified, in-line definition of molecules, i.e. a dictionary
            print("===== Feeding in a full system-definition dictionary =====")
            self.loaded_file = topdef
            self.processed_file = copy.deepcopy(topdef)
            self.process_system_dict()

    def save(self, prefix):
        # prefix = os.path.splitext(filename)[0]
        filename = "{}_processed.yaml".format(prefix)
        yaml.save_dict(
            filename,
            self.processed_file,
            header="#Processed system topology definition.\n",
        )

        # save individual chain pdb files
        for mol_name, mol in self.mol_types.items():
            mol.save("{}_chain_{}.pdb".format(prefix, mol_name))

        # save system
        if self.system.n_atoms < 10000:
            self.system.save("{}_system.pdb".format(prefix))
        else:
            self.system.save("{}_system.xyz".format(prefix))
            self.system.save("{}_system.dcd".format(prefix))

        df, bonds = self.system.topology.to_dataframe()
        df.to_csv("{}_topology.csv".format(prefix))
        np.savetxt("{}_topology_bonds.dat".format(prefix), bonds)
        # pickle unhappy with fictional elements
        # import pickle
        # pickle.dump( self.system.topology, open("{}_topology.p".format(prefix),"wb") )
        # note: mdtraj may be unhappy loading system with unknown elements

