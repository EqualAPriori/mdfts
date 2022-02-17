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
"""
from __future__ import print_function

# standard imports
from collections import OrderedDict
from operator import itemgetter

# 3rd party imports
import mdtraj
import numpy as np
import networkx as nx

# local imports

def dict_to_string(fts_param, dict_name, n_indent=1):
    s = '  '*n_indent
    s += '{}'.format(dict_name) +'{\n'
    for key, val in fts_param.items():
        if not isinstance(val,dict):
            if isinstance(val,list) or isinstance(val,tuple):
                val = ' '.join([str(x) for x in val])
            s += '  '*(n_indent+1) + '{} = {}\n'.format(key,val)
        else:
            sub_s = dict_to_string(val, key, n_indent=n_indent+1)
            s += '\n' + sub_s 
    s += '  '*n_indent + '}\n'
    return s                

class Topology(mdtraj.core.topology.Topology):
    """Create a new topology based on mdtraj topology class
    """
    
    def __init__(self):
        """Create a new Topology object based on mdtraj
            Add molecules if system_pdb (with CONECT fields) is provided

        Args:
            system_pdb (str, optional): pdb describes system. Defaults to None.
        """
        mdtraj.core.topology.Topology.__init__(self)
        self._n_per_chaintype = {}
        self._atom_types = []
        self._chain_types = []
        self._fts_param = {}

    def add_chain_from_pdb(self, pdb, n_copies, chain_name = None):
        """Create new chains and add to the topology from pdb file 
            Register molecule in pdb as a new chain type

        Args:
            pdb (str): pdb (with CONECT fields)
                recommended for a single type of molecule at a time 
                
            n_copies (int): number of copies to add to the topology
        """
        topology = mdtraj.load(pdb).topology

        if not chain_name:
            chain_name = 'chain{}'.format(len(self._n_per_chaintype))
        self._n_per_chaintype.update({chain_name: n_copies})
        self._chain_types.append(chain_name)

        for i in range(n_copies): 
            atom_map = {} 
            for chain in topology.chains:
                c = self.add_chain()
                for residue in chain.residues:
                    r = self.add_residue(residue.name, c, self.n_residues, residue.segment_id)
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
        
        self._fts_param.update({'chain{}'.format(len(self._chain_types)): self.get_fts_chain_param(topology, chain_name)})

    def add_chains_from_yaml(self, yaml):
        """Create new chains and add to the topology from yaml file

        Args:
            yaml (str): _description_
        """

    def get_fts_chain_param(self, chain_top, chain_name, chain_stat='DGC'):
        """_summary_

        Args:
            chain_top (Topology): topology of a chain
            chain_name (str)
            chain_stat (str, optional): chain statistics. Defaults to 'DGC'.

        Returns:
            dict: dictionary of FTS parameters for this chain
        """
        g = chain_top.to_bondgraph()
        nx.draw_networkx(g)
         
        # check for close-loop in chain
        has_cycle = False
        try:
            if nx.algorithms.cycles.find_cycle(g):
                has_cycle = True
        except:
            pass
        if has_cycle:
            print('Detect cycle in the topology, not compatible with PolyFTS')
            return None
        
        # check architecture of chain
        visited = {x: False for x in g.nodes}
        segments = OrderedDict()
        while not all(visited.values()):
            avail_nodes = [x for x in g.nodes if not visited[x]]
            seg = []
            seg_name = 'seg{}'.format(len(segments))
            for i, node in enumerate(avail_nodes):
                if i == 0:
                    seg.append(node)
                    visited[node] = True
                    continue
                adj = list(g.adj[node])
                a_id = [a.index for a in adj]
                adj_ind = np.vstack((a_id, adj)).transpose()
                adj_ind = np.array(sorted(adj_ind, key=itemgetter(0)))
                adj_sorted = adj_ind[:,1]
                if adj_sorted[0] == seg[-1]:
                    seg.append(node)
                    visited[node] = True
            segments.update({seg_name: tuple(seg)})

        structure = {x: [] for x in list(segments.keys())}
        for i in range(len(segments.keys()) -1 ):
            seg_name1 = list(segments.keys())[i]
            for j in range(i+1, len(segments.keys())):
                seg_name2 = list(segments.keys())[j]
                edges = nx.algorithms.boundary.edge_boundary(g, segments[seg_name1], segments[seg_name2])
                edges = [x for x in edges] # should only have at most one edge between any two linear segments for acyclic graph
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
                arch = 'point'
            else:
                arch = 'linear'
        else:
            seg0_name = list(structure.keys())[0]
            graft_pos0 = [x[1] for x in structure[seg0_name]]
            n_other_branches = sum([len(structure[x]) for x in list(structure.keys())[1:]])
            if len(structure[seg0_name]) == len(segments)-1 and all(x == graft_pos0[0] for x in graft_pos0) and n_other_branches == 0:
                arch = 'star'
                join_node = graft_pos0[0]
            else:
                arch = 'comb'

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
            param.update({'Statistics': chain_stat})
            param.update({'BlockSpecies': [self.atom_types.index(a.name)+1 for a in seg]})
            param.update({'NBeads': len(param['BlockSpecies'])})
            param.update({'NBlocks': len(param['BlockSpecies'])})
            param.update({'NperBlock': [1] * len(param['BlockSpecies'])})
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
            param.update({'BackboneGraftingPositions': graft_pos})
            param.update({'NumArms': 1})
            if len(seg_structure):
                param.update({'NumSideArmTypes': len(seg_structure)})
            for arm_i, (arm_name, arm_graft_pos) in enumerate(seg_structure):
                arm_graft_pos = seg.index(arm_graft_pos) + 1 # relative grafting position on the superior arm
                arm_param = comb_param(segments[arm_name], chain_stat, arm_graft_pos, structure[arm_name])
                param.update({'sidearmtype{}'.format(arm_i+1): arm_param})
            return param

        fts_param = {'Architecture': arch, 'Label': chain_name}
        if arch == 'point':
            fts_param.update({'Species': list(segments.values())[0]})
            # fts_parm.update('Species': segments.items()[0].name)
        elif arch == 'linear':
            param = linear_param(list(segments.values())[0], chain_stat)
            fts_param = {key: val for d in (fts_param, param) for key, val in d.items()} 
        elif arch == 'star':
            fts_param.update({'NumArmTypes': len(segments)})
            fts_param.update({'JoinBeadSpecies': self.atom_types.index(join_node.name)+1 })
            for arm_i, seg in enumerate(segments.values()):
                if arm_i == 0: #exclude join node 
                    seg1 = seg[:seg.index(join_node)]
                    arm_param = linear_param(seg1, chain_stat)
                    arm_param.update({'ArmMultiplicity': 1})
                    fts_param.update({'arm1': arm_param})

                    seg2 = seg[seg.index(join_node)+1 :]
                    arm_param = linear_param(seg2, chain_stat)
                    arm_param.update({'ArmMultiplicity': 1})
                    fts_param.update({'arm2': arm_param})
                else:
                    arm_param = linear_param(seg, chain_stat)
                    arm_param.update({'ArmMultiplicity': 1})
                    fts_param.update({'arm{}'.format(arm_i+2): arm_param})
        elif arch == 'comb':
            fts_param.update({'NumSideArmTypes': len(segments)-1})
            backbone_seg = list(segments.values())[0]
            backbone_param = linear_param(backbone_seg, chain_stat)
            backbone_structure = list(structure.values())[0]
            fts_param.update({'backbone': backbone_param})
            for arm_i, (arm_name, graft_pos) in enumerate(backbone_structure):
                graft_pos = backbone_seg.index(graft_pos) + 1
                print('Arms of {}: {}'.format(arm_name, structure[arm_name]))
                arm_param = comb_param(segments[arm_name], chain_stat, graft_pos, structure[arm_name])
                fts_param.update({'sidearmtype{}'.format(arm_i+1): arm_param})
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

    pdb_list = ['test_point.pdb', 'test_linear.pdb',
                'test_comb.pdb','test_comb2.pdb','test_star.pdb','test_cycle.pdb']
    mode = 5
    top = Topology()
    top.add_chain_from_pdb(pdb_list[mode], 3, chain_name='POLYMER')
    fts_param = top.fts_param
    print('{} atoms'.format(top.n_atoms))
    print('{} chains'.format(top.n_chains))
    print('n_per_chaintype: ',top.n_per_chaintype)
    s = dict_to_string(fts_param, 'chain', n_indent=0)
    print('\n' + s)
    plt.show() 