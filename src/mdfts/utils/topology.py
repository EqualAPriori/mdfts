"""
Read files (pdb, yaml) and create topology of the system
Check chain architecture and store chain parameters for FTS
    for comb structure: assumes the first detected linear segment is backbone        
Currently, register chains in pdb as new chain types everytime add_chain_from_pdb is called

To do:
    Where to get information about chain statistics to use in FTS, default is DGC
    Create topology from yaml file
    Export to sim
    A way to create topology from a full pdb (including all chains in the system)
        and recognize unique chain types. Could assign unique name for each chain to residue column of pdb

Save/Load:
df = pandas.read_csv(args.top[0])
bonds = np.loadtxt(args.top[1])
top_mdtraj = mdtraj.Topology.from_dataframe(df,bonds)
top = top_mdtraj.to_openmm()

top_mdtraj = mdtraj.Topology.from_openmm(top)
df,bonds = top_mdtraj.to_dataframe()
df.to_csv('{}_topology.csv'.format(prefix))
np.savetxt('{}_topology_bonds.dat'.format(prefix),bonds)

# Want the object to contain:
chain_defs, in pdb/mdtraj format
chain_defs, in pfts format

really, need to keep two chain defs formats
should this be the master topolgy? or a *list* of topologies?

store as graft_types? OR segment types, units, and grafts?

everything else
"""
from __future__ import print_function

# standard imports
from collections import OrderedDict
from operator import itemgetter
import sys
import collections

try:
    collectionsABC = collections.abc
except:
    collectionsABC = collections

# 3rd party imports
import mdtraj
import numpy as np
import networkx as nx

# local imports
from mdfts.utils import serial


def dict_to_string(fts_param, dict_name, n_indent=1):
    s = "  " * n_indent
    s += "{}".format(dict_name) + "{\n"
    for key, val in fts_param.items():
        if not isinstance(val, dict):
            if isinstance(val, list) or isinstance(val, tuple):
                val = " ".join([str(x) for x in val])
            s += "  " * (n_indent + 1) + "{} = {}\n".format(key, val)
        else:
            sub_s = dict_to_string(val, key, n_indent=n_indent + 1)
            s += "\n" + sub_s
    s += "  " * n_indent + "}\n"
    return s


def validate_list(val, types):
    if not isinstance(val, collectionsABC.Iterable):
        raise TypeError("expecting an iterable of type {}".format(types))
    tmp = []
    for x in val:
        if not isinstance(x, types):
            raise TypeError(
                "Values must be of type {}, but found {} of type {}".format(
                    types, x, type(x)
                )
            )
        else:
            tmp.append(x)
    return tmp


class Segment(serial.Serializable):
    """Segments represent *linear* pieces of chain"""

    def __init__(self, segdef=[], stat="DGC"):
        self._data = OrderedDict()
        self.statistics = stat
        self.clear()
        self._serial_vars = ["statistics", "block_species", "n_per_block"]
        self.update(segdef)

    def to_fts_dict(self):
        pass

    def clear(self):
        self._data["BlockSpecies"] = []
        self._data["NPerBlock"] = []

    def update(self, segdef):
        if isinstance(segdef, Segment):
            segdef = segdef.sequence_compact()
        elif not isinstance(segdef, collectionsABC.Iterable):
            raise TypeError(
                "segment definition must be an iterable of beadnames or (beadname,num) pairs"
            )
        self.clear()
        for blockdef in segdef:
            if isinstance(blockdef, str):
                self.add_block(blockdef, 1)
            elif isinstance(blockdef, (list, tuple)):
                name, num = blockdef[0], blockdef[1]
                self.add_block(name, num)

    def add_block(self, species_name, num=1):
        if isinstance(species_name, (int, str)):
            self._data["BlockSpecies"].append(species_name)
        else:
            raise TypeError("species_name must be str or int")
        if isinstance(num, int):
            self._data["NPerBlock"].append(num)
        else:
            raise TypeError("species_name must be str or int")

    def remove_block(self, ind=None):
        if ind is None:
            ind = -1
        species = self._data["BlockSpecies"].pop(ind)
        num = self._data["NPerBlock"].pop(ind)
        return species, num

    def to_compact(self):
        """Turn chain_def in block_species and n_per_block into a compact form"""
        seq = self.sequence_compact()
        self.update(seq)

    # === Properties ===
    @property
    def statistics(self):
        return self._data["Statistics"]

    @statistics.setter
    def statistics(self, val):
        if not isinstance(val, str):
            raise ValueError("statistics must be a string [DGC,FJC,Continuous]")
        if val.lower() in ["dgc"]:
            self._data["Statistics"] = "DGC"
        elif val.lower() in ["fjc"]:
            self._data["Statistics"] = "FJC"
        elif val.lower() in ["continuous"]:
            self._data["Statistics"] = "Continuous"

    @property
    def block_species(self):
        # will be a list of beadnames
        return self._data["BlockSpecies"]

    """
    @block_species.setter
    def block_species(self, val):
        self._data["BlockSpecies"] = validate_list(val, (str, int))
    """

    @property
    def n_per_block(self):
        return self._data["NPerBlock"]

    """
    @n_per_block.setter
    def n_per_block(self, val):
        self._data["NPerBlock"] = validate_list(val, (str, int))
    """

    @property
    def n_beads(self):
        return sum(self._data["NPerBlock"])

    @property
    def n_blocks(self):
        return len(self._data["BlockSpecies"])

    def sequence(self):
        """return sequence, verbose, i.e. ['A','A','A','B','B','B','B']

        Returns:
            tuple: bead names
        """
        seq = []
        for ind in range(self.n_blocks):
            name = self.block_species[ind]
            num = self.n_per_block[ind]
            seq.extend([name] * num)
        return tuple(seq)

    def sequence_compact(self):
        """returns sequence in most compact version of the format [('name',num),...]

        Returns:
            tuple: of ('name',num) tuples
        """
        full_seq = self.sequence()
        seq = []
        current_block_species = full_seq[0]
        num = 1
        if len(full_seq) > 1:
            for bead_name in full_seq[1:]:
                if bead_name == current_block_species:
                    num += 1
                else:
                    seq.append((current_block_species, num))
                    current_block_species = bead_name
                    num = 1
        seq.append((current_block_species, num))

        return tuple(seq)

    def sequence_compactest(self):
        """returns sequence in most compact version
        Returns:
            tuple: of ('name',num) tuples
        """
        full_seq = self.sequence()
        seq = []
        current_block_species = full_seq[0]
        num = 1
        if len(full_seq) > 1:
            for bead_name in full_seq[1:]:
                if bead_name == current_block_species:
                    num += 1
                else:
                    if num > 1:
                        seq.append((current_block_species, num))
                    else:
                        seq.append(current_block_species)
                    current_block_species = bead_name
                    num = 1

        if num > 1:
            seq.append((current_block_species, num))
        else:
            seq.append(current_block_species)
        return tuple(seq)

    def __eq__(self, other):
        if self.sequence() != other.sequence():
            return False
        if self.statistics != other.statistics:
            return False

    # === Serializable interface ===


class ArmType(Segment):
    """Essentially a *hashable* segment (based on object id) that is meant to be *unique*.
    This is what can go into networkx."""

    def __eq__(self, other):
        return id(self) == id(other)

    def __hash__(self):
        return hash(id(self) >> 4)


class Chain(serial.Serializable):
    """Container for a chain
    auto-detect architecture
    auto-detect num SideArmTypes

    c = chain(label)
    c.backbone = segment
    c.add_sidearm(unit, graftingpositions)

    u = unit(label)
    u.root = segment
    u.add_sidearm(unit or segment, graftingpositions)
    #main difference: backbone evaluates *architecture* as well
    #need to be careful: doesn't get circular

    another approach:
    c.add_segment(segment, prev segment to graft to, position, multiplicity/# grafts)
    i.e. storing a connectivity graph

    i.e. as lots of objects that store their local connections (like linked list), or one large object storing everyone's connections

    but nested multi-gen comb is not well-defined right now anyway...

    note for turning into pdb, will need multi-generation anyway...

    one way:
    list of [segments]
    tuples of (segment, grafting positions, multiplicity/graft)
        - always assuming graft is grafted at end-0
    then this can be enumerated!
    My's convention is new segment per graft, and 1 grafting position per connection, potentially with multiplicity

    but... can I reuse segments in different units?
    or does each unit need new segment definitions?
    e.g. can declare new units based on a repertoire of segments...
    and then units themselves are unique/non-repeated?

    i.e. in the format I proposed, segments are just shorthand for the blockdef... that a unit can initialize itself with
    so not that the segments/sequences actually pre-compute anything.
    most likely they need to be specially calculated anyway to account for grafts, etc.!

    multiplicity is grafts/point

    while numarms is literally... the # of arms... weird why that has to be specified.

    advantage of putting the grafting position in the sidearmtype:
    very transparent that the same propagator can be used multiple times!

    i.e. is the "reverse propagator" perspective
    either way, need to allow for *multigraph* specification!

    can have that multigraph specification be *outside* the graft definitions...
    OK! So in the end, define as multigraph,
    which should be agnostic and easily translatable to wherever we choose to
    define grafts in the PolyFTS definition!

    NOTE: in the multi-arm definition,
    inherently does not allow for cyclics! which is counter to typical graph constructions
    this is because each propagation down a graft makes a new grapt,
    instead of connecting to existing graft
    (using an existing graft only makes sense if it's a unique graft that doesn't repeat)

    i.e.
    chains = particular in []
    segments = []
    units = things with backbones made from segments, really just for shorthand
    directed edges (u1,u2)['graft']=[list of ordered grafting pairs/indices]
    """

    def __init__(self, label=None):
        if label is None:
            self.label = "CHAIN"
        self._serial_vars = ["label"]

    @property
    def architecture(self):
        pass

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, val):
        self._label = val

    @property
    def num_sidearm_types(self):
        pass


class FTSTopology(serial.Serializable):
    """
    Note:
        example usage:

        FT = FTSTopology()

        c = FT.add_chaintype([('A',10)])
        a1 = FT.add_armtype([('B',3),('C',2)])
        FT.add_graft(c,a1,[1,2,4,8],0,multiplicity = 2)
        #^^ is equivalent to: FT.add_graft(c,a1,[1,1,2,2,4,4,8,8],0,multiplicity = 1)
        a2 = FT.add_armtype(['A',2])
        FT.add_graft(a1,a2,[0,2,4])
        FT.check() #True, i.e. no infinite recursion chains

        Note, c is also an ArmType.
        Calling add_chaintype() merely tells the code that this an armtype
        that should be taken as a root/backbone when enumerating chain architectures.
    """

    def __init__(self):
        self._serial_vars = []
        self.segments = []
        self.arm_types = []
        self.chain_types = OrderedDict()
        self.g = nx.MultiDiGraph()

    def add_segment(self, seg_def):
        """add a linear segment

        Args:
            segdef (list,tuple): tuple of tuples or names

        Note:
            segdef can be ['A','B',('C',2),'D'] etc.
        """
        seg = Segment(seg_def)
        self.segments.append(seg)
        return seg

    def add_armtype(self, arm_def):
        """Given arm definition, create a *new* armtype instance and store it."""
        if isinstance(arm_def, int):
            """if int, arm_def should be index of the segdef we want"""
            seg_def = self.segments[arm_def]
            arm = ArmType(seg_def)
        else:
            arm = ArmType(arm_def)
            seq = arm.sequence()
            no_match = True
            for s in self.segments:
                if s.sequence() == seq:
                    no_match = False
            if no_match:
                self.add_segment(seq)

        self.arm_types.append(arm)
        self.g.add_node(self.arm_types.index(arm))
        return arm

    def add_chaintype(self, chain_def, name=None):
        arm = self.add_armtype(chain_def)
        if name is None:
            name = "Chain{}".format(len(self.chain_types) + 1)
        self.chain_types[name] = arm
        return arm

    def add_graft(
        self,
        base_arm,
        graft_arm,
        base_graft_points,
        graft_graft_point=0,
        multiplicity=1,  # number of times each graft enumerated in base_graft_points is enumerated
    ):
        """the arms must be pre-existing!"""
        validate_list(base_graft_points, int)
        if graft_graft_point not in [0, -1]:
            raise ValueError(
                "The grafted chain should be attached from one of its ends. Default is 0."
            )
        graft_attr = OrderedDict(
            {
                "graft_at": base_graft_points,
                "graft_from": graft_graft_point,
                "multiplicity": multiplicity,
            }
        )
        if isinstance(base_arm, ArmType):
            u = self.arm_types.index(base_arm)
        elif isinstance(base_arm, int):
            u = base_arm
        else:
            raise TypeError("base_arm must be specified by index or by instance")
        if isinstance(graft_arm, ArmType):
            v = self.arm_types.index(graft_arm)
        elif isinstance(graft_arm, int):
            v = graft_arm
        else:
            raise TypeError("graft_arm must be specified by index or by instance")
        self.g.add_edge(u, v, **graft_attr)

    def get_armtype_index(self, arm_def, add_if_not_found=False):
        """return armtype index

        Args:
            arm_def (ArmType,int): anything that can be made into an arm
            add_if_not_found (bool, optional): Whether or not to add arm_type if definition is not found

        Raises:
            ValueError: if index exceeds current list
            ValueError: if arm not found and add_if_not_found = False
            TypeError: if not a va

        Returns:
            int: index of requested arm_type

        Note:
            search through existing arm_types and return index if found.
            if add_if_not_found, create arm_type if not found.
        """
        if isinstance(arm_def, int):
            if arm_def > len(self.arm_types) - 1:
                raise ValueError("index exceeds arm_types list")
            else:
                u = arm_def
        elif isinstance(arm_def, ArmType):
            if arm_def not in self.arm_types:
                if add_if_not_found:  # then add the arm definition
                    arm_def = self.add_armtype(arm_def)
                else:
                    raise ValueError("arm not found")
            u = self.arm_types.index(arm_def)
        else:  # i.e. a fresh definition
            if add_if_not_found:
                arm_def = self.add_armtype(arm_def)
                u = self.arm_types.index(arm_def)
            else:
                raise TypeError("arm must be specified by index or by instance")
        return u

    def add_path(self, root, *args, mode=0):
        """Shorthand for defining a chain, with only one new kind of chain attached to each branch

        nesting allowed!?

        Args:
            branch_defs (list): of lists or tuple of tuples

        Note:
            ASSUMES:
                always grafting from end 0, multiplicity = 1

            each element of branch_defs should be in the following format:
                (node,[graft_points],[potential more subgrafts])

            Nested syntax:
            [ (armtype1,graftingpoints,[(subarmtype,graftingpoints)]),
            (armtype2,graftingpoints) ]

            mode = 0:
                branch_defs should refer to EXISTING arm_types.
                this is potentially unsafe, as you might unwittingly change a branchdef that another existing chain depends on
            mode = 1:
                branch_defs should create *new* arm_types
                this is always safe, but may end up creating a lot of arms, as it doesn't detect for redundancy.
                exception is that the root always tries to build off of existing node
        """
        if mode not in [0, 1]:
            raise ValueError("Unknown mode for add_path()")
        branch_defs = args

        if mode == 0:
            try:
                u = self.get_armtype_index(root, add_if_not_found=False)
            except:
                raise ValueError("must use pre-defined armtype")
        elif mode == 1:
            u = self.get_armtype_index(root, add_if_not_found=True)

        for branch_def in branch_defs:
            print("branch_def {}".format(branch_def))
            proposed_graft = branch_def[0]
            base_graft_points = branch_def[1]

            if mode == 0:
                try:
                    v = self.get_armtype_index(proposed_graft, add_if_not_found=False)
                except:
                    raise ValueError("must use pre-defined armtype")
            elif mode == 1:  # always create new graft arm_type
                v = self.add_armtype(proposed_graft)

            self.add_graft(u, v, base_graft_points)

            if len(branch_def) > 2:  # need to use recursion to add future generations
                print("found subgraft")
                sub_grafts = branch_def[2:]
                for sg in sub_grafts:
                    print(sub_grafts)
                    self.add_path(v, sg, mode=mode)

    def to_verbose(self):
        """Outputs a topology/graph where every branch is explicitly enumerated"""
        pass

    def to_pdb(self):
        pass

    def visualize(self):
        pass

    def isvalid(self):
        try:
            cycle_paths = nx.find_cycle(self.g)
            return True
        except:
            return False


class Topology(mdtraj.core.topology.Topology):
    """Create a new topology based on mdtraj topology class"""

    def __init__(self):
        """Create a new Topology object based on mdtraj
        if fts_chains exist, make pdbs out of them"""
        mdtraj.core.topology.Topology.__init__(self)
        self._n_per_chaintype = {}
        self._atom_types = []
        self._chain_types = []
        self._fts_param = {}
        self._fts_segments = []
        self._fts_units = []
        self._fts_chains = []

    def save(self, prefix=None):
        df, bonds = self.to_dataframe()
        df.to_csv("{}_topology.csv".format(prefix))
        np.savetxt("{}_topology_bonds.dat".format(prefix), bonds)

    def load(self, beadfile, bondfile):
        df = pandas.read_csv(beadfile)
        bonds = np.loadtxt(bondfile)
        top_mdtraj = mdtraj.Topology.from_dataframe(df, bonds)
        self = top_mdtraj  # does this erase everything?

    def add_chain_from_pdb(self, pdb, n_copies, chain_name=None):
        """Create new chains and add to the topology from pdb file
            Register molecule in pdb as a new chain type

        Args:
            pdb (str): pdb (with CONECT fields)
                recommended for a single type of molecule at a time

            n_copies (int): number of copies to add to the topology
        """
        topology = mdtraj.load(pdb).topology

        if not chain_name:
            chain_name = "chain{}".format(len(self._n_per_chaintype))
        self._n_per_chaintype.update({chain_name: n_copies})
        self._chain_types.append(chain_name)

        for i in range(n_copies):
            atom_map = {}
            for chain in topology.chains:
                c = self.add_chain()
                for residue in chain.residues:
                    r = self.add_residue(
                        residue.name, c, self.n_residues, residue.segment_id
                    )
                    for atom in residue.atoms:
                        self.add_atom(atom.name, atom.element, r, serial=self.n_atoms)
                        atom_map.update({atom: self.atom(-1)})
                        if not atom.name in self.atom_types:
                            self.atom_types.append(atom.name)

            for bond in topology.bonds:
                a1, a2 = bond
                a1 = atom_map[a1]
                a2 = atom_map[a2]
                self.add_bond(a1, a2, type=bond.type, order=bond.order)

        self._fts_param.update(
            {
                "chain{}".format(len(self._chain_types)): self.get_fts_chain_param(
                    topology, chain_name
                )
            }
        )

    def add_chains_from_yaml(self, yaml):
        """Create new chains and add to the topology from yaml file

        Args:
            yaml (str): _description_
        """

    def get_fts_chain_param(self, chain_top, chain_name, chain_stat="DGC"):
        """Check chain architecture and return FTS parameters for this chain

        Args:
            chain_top (Topology): topology of a chain
            chain_name (str)
            chain_stat (str, optional): chain statistics. Defaults to 'DGC'.

        Returns:
            dict: dictionary of FTS parameters for this chain
        """
        g = chain_top.to_bondgraph()
        nx.draw_networkx(g)

        # check for closed-loop in chain
        has_cycle = False
        try:
            if nx.algorithms.cycles.find_cycle(g):
                has_cycle = True
        except:
            pass
        if has_cycle:
            print("Detect cycle in the topology, not compatible with PolyFTS")
            return None

        # check architecture of chain
        visited = {x: False for x in g.nodes}
        segments = OrderedDict()
        while not all(visited.values()):
            # algorithm: find segments of the chain that are contiguous in the node-list-order
            avail_nodes = [x for x in g.nodes if not visited[x]]
            seg = []
            seg_name = "seg{}".format(len(segments))
            for i, node in enumerate(avail_nodes):
                if i == 0:
                    seg.append(node)
                    visited[node] = True
                    continue
                adj = list(g.adj[node])
                a_id = [a.index for a in adj]
                adj_ind = np.vstack((a_id, adj)).transpose()
                adj_ind = np.array(sorted(adj_ind, key=itemgetter(0)))
                adj_sorted = adj_ind[:, 1]
                if adj_sorted[0] == seg[-1]:
                    seg.append(node)
                    visited[node] = True
            segments.update({seg_name: tuple(seg)})

        structure = {x: [] for x in list(segments.keys())}
        for i in range(len(segments.keys()) - 1):
            # for each segment, find out what other segmets are grafted to it.
            # does this assume that grafts are always attached at their ends?
            # is that always true? and which end is grafted?
            seg_name1 = list(segments.keys())[i]
            for j in range(i + 1, len(segments.keys())):
                seg_name2 = list(segments.keys())[j]
                edges = nx.algorithms.boundary.edge_boundary(
                    g, segments[seg_name1], segments[seg_name2]
                )
                edges = [
                    x for x in edges
                ]  # should only have at most one edge between any two linear segments for acyclic graph
                if edges:
                    edges = edges[0]
                    if edges[0] in segments[seg_name1]:
                        graft_node = edges[0]
                    else:
                        graft_node = edges[1]
                    structure[seg_name1].append((seg_name2, graft_node))

        if len(structure) == 1:
            seg_name = list(structure.keys())[0]
            if len(segments[seg_name]) == 1:
                arch = "point"
            else:
                arch = "linear"
        else:
            seg0_name = list(structure.keys())[0]
            graft_pos0 = [x[1] for x in structure[seg0_name]]
            n_other_branches = sum(
                [len(structure[x]) for x in list(structure.keys())[1:]]
            )
            if (
                len(structure[seg0_name]) == len(segments) - 1
                and all(x == graft_pos0[0] for x in graft_pos0)
                and n_other_branches == 0
            ):
                arch = "star"
                join_node = graft_pos0[0]
            else:
                arch = "comb"

        # write dictionary for chain in FTS language
        def linear_param(seg, chain_stat):
            """FTS parameters for linear

            Args:
                seg (list): list of atoms in this segment
                chain_stat (str): chain statistics

            Returns:
                dict: dictionary of FTS parameters
            """
            param = {}
            param.update({"Statistics": chain_stat})
            param.update(
                {"BlockSpecies": [self.atom_types.index(a.name) + 1 for a in seg]}
            )
            param.update({"NBeads": len(param["BlockSpecies"])})
            param.update({"NBlocks": len(param["BlockSpecies"])})
            param.update({"NperBlock": [1] * len(param["BlockSpecies"])})
            return param

        def comb_param(seg, chain_stat, graft_pos, seg_structure):
            """FTS parameters for comb

            Args:
                seg (list): list of atoms in this segment
                chain_stat (str): chain statistics
                graft_pos (int): relative grafting position (base 1) on the segment where this arm is attached
                seg_structure (list): pairs of (name of arm attached to this segment, graft node on this segment)

            Returns:
                dict: dictionary of FTS parameters
            """
            param = linear_param(seg, chain_stat)
            param.update({"BackboneGraftingPositions": graft_pos})
            param.update({"NumArms": 1})
            if len(seg_structure):
                param.update({"NumSideArmTypes": len(seg_structure)})
            for arm_i, (arm_name, arm_graft_pos) in enumerate(seg_structure):
                arm_graft_pos = (
                    seg.index(arm_graft_pos) + 1
                )  # relative grafting position on the superior arm
                arm_param = comb_param(
                    segments[arm_name], chain_stat, arm_graft_pos, structure[arm_name]
                )
                param.update({"sidearmtype{}".format(arm_i + 1): arm_param})
            return param

        fts_param = {"Architecture": arch, "Label": chain_name}
        if arch == "point":
            fts_param.update({"Species": list(segments.values())[0]})
        elif arch == "linear":
            param = linear_param(list(segments.values())[0], chain_stat)
            fts_param = {key: val for d in (fts_param, param) for key, val in d.items()}
        elif arch == "star":
            fts_param.update({"NumArmTypes": len(segments)})
            fts_param.update(
                {"JoinBeadSpecies": self.atom_types.index(join_node.name) + 1}
            )
            for arm_i, seg in enumerate(segments.values()):
                if arm_i == 0:  # exclude join node
                    seg1 = seg[: seg.index(join_node)]
                    arm_param = linear_param(seg1, chain_stat)
                    arm_param.update({"ArmMultiplicity": 1})
                    fts_param.update({"arm1": arm_param})

                    seg2 = seg[seg.index(join_node) + 1 :]
                    arm_param = linear_param(seg2, chain_stat)
                    arm_param.update({"ArmMultiplicity": 1})
                    fts_param.update({"arm2": arm_param})
                else:
                    arm_param = linear_param(seg, chain_stat)
                    arm_param.update({"ArmMultiplicity": 1})
                    fts_param.update({"arm{}".format(arm_i + 2): arm_param})
        elif arch == "comb":
            fts_param.update({"NumSideArmTypes": len(segments) - 1})
            backbone_seg = list(segments.values())[0]
            backbone_param = linear_param(backbone_seg, chain_stat)
            backbone_structure = list(structure.values())[0]
            fts_param.update({"backbone": backbone_param})
            for arm_i, (arm_name, graft_pos) in enumerate(backbone_structure):
                graft_pos = backbone_seg.index(graft_pos) + 1
                print("Arms of {}: {}".format(arm_name, structure[arm_name]))
                arm_param = comb_param(
                    segments[arm_name], chain_stat, graft_pos, structure[arm_name]
                )
                fts_param.update({"sidearmtype{}".format(arm_i + 1): arm_param})
        return fts_param

    @property
    def n_per_chaintype(self):
        """Get the number of chains for each chain type in the Topology"""
        return self._n_per_chaintype

    @property
    def chain_types(self):
        """Get the number of chain types in the Topology"""
        return self._chain_types

    @property
    def atom_types(self):
        """Get the number of atom types in the Topology"""
        return self._atom_types

    @property
    def fts_param(self):
        """Get the dictionary of FTS parameters for all chain types"""
        return self._fts_param


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    pdb_list = [
        "test_point.pdb",
        "test_linear.pdb",
        "test_comb.pdb",
        "test_comb2.pdb",
        "test_star.pdb",
        "test_cycle.pdb",
    ]
    """
    mode = int(sys.argv[1])
    top = Topology()
    top.add_chain_from_pdb(pdb_list[mode], 3, chain_name="POLYMER")
    fts_param = top.fts_param
    print("{} atoms".format(top.n_atoms))
    print("{} chains".format(top.n_chains))
    print("n_per_chaintype: ", top.n_per_chaintype)
    s = dict_to_string(fts_param, "chain", n_indent=0)
    print("\n" + s)
    plt.show()
    """
    s = Segment(["A", "B", ("A", 3)])

    FT = FTSTopology()

    c = FT.add_chaintype([("A", 10)])
    a1 = FT.add_armtype([("B", 3), ("C", 2)])
    FT.add_graft(c, a1, [1, 2, 4, 8], 0, multiplicity=2)
    # ^^ is equivalent to: FT.add_graft(c,a1,[1,1,2,2,4,4,8,8],0,multiplicity = 1)
    a2 = FT.add_armtype(["A", 2])
    FT.add_graft(a1, a2, [0, 2, 4])
    FT.isvalid()  # True, i.e. no infinite recursion chains
