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
import sys, warnings
import collections

try:
    collectionsABC = collections.abc
except:
    collectionsABC = collections

# 3rd party imports
import mdtraj
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


# local imports
from mdfts.utils import serial
from mdfts.utils import yamlhelper


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

    @block_species.setter
    def block_species(self, val):
        self._data["BlockSpecies"] = validate_list(val, (str, int))

    @property
    def n_per_block(self):
        return self._data["NPerBlock"]

    @n_per_block.setter
    def n_per_block(self, val):
        self._data["NPerBlock"] = validate_list(val, (str, int))

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
        self._serial_vars = [
            "chain_types",
            "arm_types",
            "segments",
            "g",
        ]
        self.segments = serial.SerializableTypedList(Segment)
        self.arm_types = serial.SerializableTypedList(ArmType)
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
        self.g.add_node(self.arm_types.index(arm), arm=arm)
        return arm

    def add_chaintype(self, chain_def, name=None):
        arm = self.add_armtype(chain_def)
        if name is None:
            name = "Chain{}".format(len(self.chain_types) + 1)
        self.chain_types[name] = self.arm_types.index(arm)
        return arm

    def add_graft(
        self,
        base_arm,
        graft_arm,
        base_graft_points,
        graft_graft_point=0,
        multiplicity=1,  # number of times each graft enumerated in base_graft_points is enumerated
    ):
        """Add a graft to the graph (via directed edge), with grafting information
        Note:
            The arms must be pre-existing!
            Consider raising an error if an existing directed edge is already defined.
                This encourages cleaner definitions.
        """
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
        if (u, v) in self.g.edges():
            warnings.warn(
                "Trying to add edge definition to already existing edge. Probably wrong!"
            )
        self.g.add_edge(u, v, **graft_attr)

    def get_armtype_index(self, arm_def, add_if_not_found=False):
        """return armtype index

        Args:
            arm_def (ArmType,int): anything that can be made into an arm
            add_if_not_found (bool, optional): Whether or not to add arm_type if definition is not found

        Raises:
            ValueError: if index exceeds current list
            ValueError: if arm not found and add_if_not_found = False
            TypeError: if not a valid arm_type

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

    def add_path(self, root, *args, mode=0, as_chain_name=False):
        """Shorthand for defining a chain, with only one new kind of chain attached to each branch

        nesting allowed!?

        Args:
            branch_defs (list): of lists or tuple of tuples
            as_chain_name (bool,str): whether or not to label `root` as a chain

        Returns:
            (int): index of root where things started
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
        if as_chain_name:
            self.chain_types[as_chain_name] = u

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

        return u

    def get_grafts(self, root, enum=False):
        """Returns a nested list of grafts that FTSToplogy().add_path() can use
        Fully enumerates all multiplicities and each graft_point

        Args:
            self (FTSTopology): with pre-existing chains
            root (ArmType or int): root arm to start enumerating from

        Returns:
            list: of lists
        """
        if isinstance(root, ArmType):
            root = self.arm_types.index(root)

        grafts = [(e1, d) for e0, e1, d in self.g.edges(data=True) if e0 == root]

        graft_defs = []
        for e, d in grafts:
            if d["graft_from"] != 0:
                raise ValueError(
                    "Chain enumeration currently only works if grafting from bead 0"
                )

            for ii in range(d["multiplicity"]):
                for graft_point in d["graft_at"]:
                    graft_def = [
                        self.arm_types[e].sequence_compactest(),
                        [graft_point],
                    ]
                    graft_def.extend(self.get_grafts(e))
                    if len(graft_def) > 0:
                        graft_defs.append(graft_def)
        return graft_defs

    def fully_enumerate(self):
        """Outputs a topology/graph where every branch is explicitly enumerated

        only works if:
        - acyclic
        - no directed edge is defined more than once
        """
        ft = FTSTopology()

        for chain_name, arm_index in self.chain_types.items():
            arm_type = self.arm_types[arm_index]

            graft_def = self.get_grafts(arm_index)
            u = ft.add_path(arm_type.sequence_compactest(), *graft_def, mode=1)
            ft.chain_types[chain_name] = u

        return ft

    def expand_graft_to_beads(
        self, arm_ind, graft_index=-1, node_list=None, edge_list=None
    ):
        """
        add on arm_type(arm_ind) to node and edge lists
        then for each graft off of the arm, expand

        Assumes arms are all attached at their beginning index!
        """
        if graft_index == -1:  # this arm not grafted to anything, i.e. is a backbone
            fresh_chain = True
            node_list = []
            edge_list = []
        else:
            fresh_chain = False

        bead_names = self.arm_types[arm_ind].sequence()
        start_ind = len(node_list)
        end_ind = start_ind + len(bead_names)

        # add this new arm
        bead_inds = list(range(start_ind, end_ind))
        node_list.extend(bead_names)

        if graft_index <= -1:
            new_edges = []
        else:
            new_edges = [(graft_index, bead_inds[0])]
        new_edges.extend((ii - 1, ii) for ii in bead_inds[1:])
        edge_list.extend(new_edges)

        # start adding child generations
        grafts = [(e1, d) for e0, e1, d in self.g.edges(data=True) if e0 == arm_ind]

        for e, d in grafts:
            if d["graft_from"] != 0:
                raise ValueError(
                    "Chain enumeration currently only works if grafting from bead 0"
                )

            for ii in range(d["multiplicity"]):
                for graft_point in d["graft_at"]:
                    self.expand_graft_to_beads(
                        e, bead_inds[graft_point], node_list, edge_list
                    )

        if fresh_chain:
            return node_list, edge_list

    def expand_to_beads(self):
        """Generate bead-by-bead definition, i.e. full connectivity graph

        Notes:
            Todo: also generate random walk initial coordinates!
            Algorithm:
                1. fully enumerate all grafts and multiplicities
                2. start with backbone, write out all beads and add edges
                3. for each branch
                    write out the branch fully, add edges
                        for each subbranch, write fully, add edges
                            ... recursion?

            Question: do we want each chain indexed to zero (e.g. for an independent pdb)
                or do we want the index to build on one another?
        """
        chains = OrderedDict()
        for chain_name, arm_index in self.chain_types.items():
            ch = self.expand_graft_to_beads(arm_index)
            # chains.append(ch)
            chains[chain_name] = ch
        return chains

    def visualize(self, detailed=False):
        """Visualize the armtype information flow
        Does NOT do well if there are multiple edges defined on a pair of nodes!
        But typical graphs should technically *never* have that happen.
        """
        # SIMPLE plot:
        # nx.draw_networkx(self.g)
        # plt.show()
        edge_labels = {}
        # this assumes all the edges have the same labels 'marks' and 'cable_name'
        for u, v, data in self.g.edges(data=True):
            label = "@{} x {}\nfrom {}".format(
                data["graft_at"], data["multiplicity"], data["graft_from"]
            )
            if (u, v) in edge_labels:
                edge_labels[u, v] += "\n" + label

            else:
                edge_labels[u, v] = label
            # data["label"] = label  # for pydot

        node_labels = {}
        for u in self.g.nodes():
            if detailed:
                node_labels[u] = self.arm_types[u].sequence_compactest()
            else:
                node_labels[u] = u

        # pos = nx.spring_layout(self.g)
        # nx.draw_networkx(self.g, pos=pos, labels=node_labels)
        # nx.draw_networkx_edge_labels(self.g, pos=pos, edge_labels=edge_labels)

        # from digraph_helper import my_draw_networkx_edge_labels
        arc_rad = 0.25
        G = self.g
        pos = nx.kamada_kawai_layout(G)
        for u, v, d in G.edges(data=True):
            d["weights"] = 10.0
        pos = nx.spring_layout(G, pos=pos, k=2)
        for u, v, d in G.edges(data=True):
            d.pop("weights")

        fig, ax = plt.subplots()
        if detailed:
            nx.draw_networkx_nodes(G, pos, ax=ax)
            nx.draw_networkx_labels(G, pos, ax=ax, labels=node_labels)
        else:
            nx.draw_networkx_nodes(G, pos, ax=ax)
            nx.draw_networkx_labels(G, pos, ax=ax, labels=node_labels)
        curved_edges = [
            edge for edge in G.edges() if tuple(reversed(edge)) in G.edges()
        ]
        straight_edges = list(set(G.edges()) - set(curved_edges))
        nx.draw_networkx_edges(G, pos, ax=ax, edgelist=straight_edges)
        nx.draw_networkx_edges(
            G,
            pos,
            ax=ax,
            edgelist=curved_edges,
            connectionstyle=f"arc3, rad = {arc_rad}",
        )
        if detailed:
            curved_edge_labels = {edge: edge_labels[edge] for edge in curved_edges}
            straight_edge_labels = {edge: edge_labels[edge] for edge in straight_edges}
            my_draw_networkx_edge_labels(
                G, pos, ax=ax, edge_labels=curved_edge_labels, rotate=False, rad=arc_rad
            )
            nx.draw_networkx_edge_labels(
                G, pos, ax=ax, edge_labels=straight_edge_labels, rotate=False
            )
        plt.show()

    def isvalid(self):
        """Return False if there are cycles (i.e. infinite recursion)

        Note:
            Can consider returning False if a single directed edge is multiply defined.
        """
        try:
            cycle_paths = nx.find_cycle(self.g)
            return False
        except:
            return True

    def is_equivalent(self, other):
        # Todo
        raise NotImplementedError

    def save(self, filename="top.yaml"):
        yamlhelper.save_dict(filename, self.to_dict())

    def load(self, filename="top.yaml"):
        topdef = yamlhelper.load(filename)
        self.from_dict(topdef)

    # === SERIALIZABLE INTERFACE ===
    def custom_get(self, k):
        if k == "g":
            nodes = list(self.g.nodes)
            edges = list(self.g.edges(data=True))
            return {"nodes": nodes, "edges": edges}
        else:
            return super(FTSTopology, self).custom_get(k)

    def custom_set(self, k, val):
        if k == "g":
            if isinstance(val, collectionsABC.Mapping):
                if "nodes" in val and "edges" in val:
                    self.g = nx.MultiDiGraph()
                    self.g.add_nodes_from(val["nodes"])
                    self.g.add_edges_from(val["edges"])
                else:
                    raise ValueError("did not receive proper node, edge data")
            else:
                raise ValueError("did not receive proper node, edge data")
        else:
            super(FTSTopology, self).custom_set(k, val)


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


def my_draw_networkx_edge_labels(
    G,
    pos,
    edge_labels=None,
    label_pos=0.5,
    font_size=10,
    font_color="k",
    font_family="sans-serif",
    font_weight="normal",
    alpha=None,
    bbox=None,
    horizontalalignment="center",
    verticalalignment="center",
    ax=None,
    rotate=True,
    clip_on=True,
    rad=0,
):
    """Draw edge labels.
    https://stackoverflow.com/questions/22785849/drawing-multiple-edges-between-two-nodes-with-networkx

    Parameters
    ----------
    G : graph
        A networkx graph

    pos : dictionary
        A dictionary with nodes as keys and positions as values.
        Positions should be sequences of length 2.

    edge_labels : dictionary (default={})
        Edge labels in a dictionary of labels keyed by edge two-tuple.
        Only labels for the keys in the dictionary are drawn.

    label_pos : float (default=0.5)
        Position of edge label along edge (0=head, 0.5=center, 1=tail)

    font_size : int (default=10)
        Font size for text labels

    font_color : string (default='k' black)
        Font color string

    font_weight : string (default='normal')
        Font weight

    font_family : string (default='sans-serif')
        Font family

    alpha : float or None (default=None)
        The text transparency

    bbox : Matplotlib bbox, optional
        Specify text box properties (e.g. shape, color etc.) for edge labels.
        Default is {boxstyle='round', ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0)}.

    horizontalalignment : string (default='center')
        Horizontal alignment {'center', 'right', 'left'}

    verticalalignment : string (default='center')
        Vertical alignment {'center', 'top', 'bottom', 'baseline', 'center_baseline'}

    ax : Matplotlib Axes object, optional
        Draw the graph in the specified Matplotlib axes.

    rotate : bool (deafult=True)
        Rotate edge labels to lie parallel to edges

    clip_on : bool (default=True)
        Turn on clipping of edge labels at axis boundaries

    Returns
    -------
    dict
        `dict` of labels keyed by edge

    Examples
    --------
    >>> G = nx.dodecahedral_graph()
    >>> edge_labels = nx.draw_networkx_edge_labels(G, pos=nx.spring_layout(G))

    Also see the NetworkX drawing examples at
    https://networkx.org/documentation/latest/auto_examples/index.html

    See Also
    --------
    draw
    draw_networkx
    draw_networkx_nodes
    draw_networkx_edges
    draw_networkx_labels
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if ax is None:
        ax = plt.gca()
    if edge_labels is None:
        labels = {(u, v): d for u, v, d in G.edges(data=True)}
    else:
        labels = edge_labels
    text_items = {}
    for (n1, n2), label in labels.items():
        (x1, y1) = pos[n1]
        (x2, y2) = pos[n2]
        (x, y) = (
            x1 * label_pos + x2 * (1.0 - label_pos),
            y1 * label_pos + y2 * (1.0 - label_pos),
        )
        pos_1 = ax.transData.transform(np.array(pos[n1]))
        pos_2 = ax.transData.transform(np.array(pos[n2]))
        linear_mid = 0.5 * pos_1 + 0.5 * pos_2
        d_pos = pos_2 - pos_1
        rotation_matrix = np.array([(0, 1), (-1, 0)])
        ctrl_1 = linear_mid + rad * rotation_matrix @ d_pos
        ctrl_mid_1 = 0.5 * pos_1 + 0.5 * ctrl_1
        ctrl_mid_2 = 0.5 * pos_2 + 0.5 * ctrl_1
        bezier_mid = 0.5 * ctrl_mid_1 + 0.5 * ctrl_mid_2
        (x, y) = ax.transData.inverted().transform(bezier_mid)

        if rotate:
            # in degrees
            angle = np.arctan2(y2 - y1, x2 - x1) / (2.0 * np.pi) * 360
            # make label orientation "right-side-up"
            if angle > 90:
                angle -= 180
            if angle < -90:
                angle += 180
            # transform data coordinate angle to screen coordinate angle
            xy = np.array((x, y))
            trans_angle = ax.transData.transform_angles(
                np.array((angle,)), xy.reshape((1, 2))
            )[0]
        else:
            trans_angle = 0.0
        # use default box of white with white border
        if bbox is None:
            bbox = dict(boxstyle="round", ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0))
        if not isinstance(label, str):
            label = str(label)  # this makes "1" and 1 labeled the same

        t = ax.text(
            x,
            y,
            label,
            size=font_size,
            color=font_color,
            family=font_family,
            weight=font_weight,
            alpha=alpha,
            horizontalalignment=horizontalalignment,
            verticalalignment=verticalalignment,
            rotation=trans_angle,
            transform=ax.transData,
            bbox=bbox,
            zorder=1,
            clip_on=clip_on,
        )
        text_items[(n1, n2)] = t

    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        left=False,
        labelbottom=False,
        labelleft=False,
    )

    return text_items


if __name__ == "__main__":
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
    FT.add_graft(c, a1, [1, 2, 4, 8], 0, multiplicity=1)
    FT.add_graft(c, a1, [1, 2, 4, 8], 0, multiplicity=1)
    # ^^ is equivalent to: FT.add_graft(c,a1,[1,1,2,2,4,4,8,8],0,multiplicity = 1)
    a2 = FT.add_armtype([("A", 2)])
    FT.add_graft(a1, a2, [0, 2, 4])
    print(FT.isvalid())  # True, i.e. no infinite recursion chains
    FT.add_graft(a2, a1, [0], -1)  # this leads to an invalid, cyclic graph
    print(FT.isvalid())
    FT.visualize(True)

    # shorthand, should be same as above!
    FT2 = FTSTopology()
    FT2.add_path(
        [("Aaaa", 10)],
        [
            (("Bbbb", 3), ("Cccc", 2)),
            [1, 1, 2, 2, 4, 4, 8, 8],
            [[("Aaaa", 2)], [0, 2, 4]],
        ],
        mode=1,
    )
    FT2.add_graft(2, 1, [0], -1)
    FT2.visualize(True)

    # Test
    FT3 = FTSTopology()
    u = FT3.add_path(
        [("Aaaa", 10)],
        [
            (("Bbbb", 3), ("Cccc", 2)),
            [1, 2, 4, 8],
            [[("Aaaa", 2)], [0, 2, 4]],
        ],
        mode=1,
        as_chain_name="Ch",
    )

    FT4 = FTSTopology()
    graft_def = FT3.get_grafts(0)
    u = FT4.add_path(FT3.arm_types[0], *graft_def, mode=1, as_chain_name="Ch")

    FT4.to_dict() == FT3.to_dict()

    FT5 = FT4.fully_enumerate()
    FT5.to_dict() == FT4.to_dict()

    FT6 = FTSTopology()
    u = FT6.add_path(
        [("Aaaa", 10)],
        [
            (("Bbbb", 3), ("Cccc", 2)),
            [1, 1, 8, 8],
            [[("Aaaa", 2)], [0, 2, 4]],
        ],
        mode=1,
        as_chain_name="Ch",
    )
    FT7 = FT6.fully_enumerate()
    FT7.visualize()
