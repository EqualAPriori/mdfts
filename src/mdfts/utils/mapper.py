# (C) Kevin Shen, 2021, kevinshen@ucsb.edu
# Utility functions for doing mappings
# 
# Main function types:
# 1) I/O
# 2) utilities for generating the index mappings
# 3) utililities for actually doing the mapping of a configuration/trajectory
#
# Todo:
# - add parallelization
# - currently needs to keep full indexing list in memory... can be expensive, not the best for large systems
#
# Envisioned workflow:
# 1) Map single chain:
#       python mapper.py singlechain_pdb params -prefix desiredprefix -single
#       and repeat for each individual system. should save out mapped pdb, mapping pdb, and yaml file for the chain
# 2) Map system:
#       python mapper.py system_dcd params -top top_for_pdb -prefix desiredprefix -system
#       technically I've also allowed for a shorthand s.t. can skip the single chain intermediary above if desired
#
# instead of setting up a whole new system file... can also take multiple arguments from command line to build up the system! even easier?
#

import numpy as np
import mdtraj
import yamlhelper as yaml
import topologify
import os
from collections import OrderedDict
import parsify

# ===== File I/O =====

def save(prefix,aa_indices_in_cg,cg_site_of_aa,shorthand=None,traj_mapped=None, traj_mapping=None, unmapped=None):
    mappingfile = prefix + '_mapping.yaml'
    if traj_mapped is not None:
        pdbfile = prefix + '_mapped.pdb'
        traj_mapped.save(pdbfile)
    if traj_mapping is not None:
        pdbfile_mapping = prefix + '_mapping.pdb'
        traj_mapping.save(pdbfile_mapping)
    else: pdbfile_mapping = None
    d = {
            'pdbfile_mapped': pdbfile,
            'pdbfile_mapping': pdbfile_mapping,
            'pdbfile_unmapped': unmapped,
            'aa_indices_in_cg': aa_indices_in_cg,
            'cg_site_of_aa': cg_site_of_aa,
            'shorthand': shorthand
        }

    yaml.save_dict( mappingfile, d, header = '# mapping specifications for molecule {}\n'.format(prefix) )
    #with open(mappingfile,'w') as f:
    #    f.write('# mapping specifications for molecule {}\n'.format(prefix))
    #    yaml.dump(d,f)


# ===== Mapping Utilities =====
def generate_pdb_mapping(traj,mapping):
    '''
    take in a cg_site_of_aa mapping array, and generate a new mapping topology for the chain, annotated with cgbead site
    Parameters
    ----------
    traj : mdtraj trajectory
    mapping : array
        should be cg_site_of_aa array format, of a single chain
    '''
    mapping_top = mdtraj.Topology()
    ch_previous = None
    sites = []
    site_to_res_map = {}
    for ia,atom in enumerate(traj.top.atoms):
        cg_beadname, cg_identifier = mapping[ia]
        if ia == 0:
            #if atom.chain != ch_previous: #techincally I envision this only working for single chain first
            ch_previous = mapping_top.add_chain()
        if atom.residue.chain != traj.top.atom(0).residue.chain:
            raise ValueError('Multiple chains detected, not implemented yet')
        if cg_identifier not in sites: #need to add a new residue
            sites.append(cg_identifier)
            res = mapping_top.add_residue(cg_beadname,ch_previous)
            site_to_res_map[cg_identifier] = res
        
        res = site_to_res_map[cg_identifier]
        new_atom = mapping_top.add_atom(atom.name, atom.element, res)

    for bond in traj.top.bonds:
        a1,a2 = mapping_top.atom(bond[0].index), mapping_top.atom(bond[1].index)
        mapping_top.add_bond( a1,a2 )

    mapping_traj = mdtraj.Trajectory( traj.xyz, mapping_top )
    return mapping_traj
 
def generate_single_mapping_shorthand(mode,shorthand):
    '''
    A helper function to generate requisite data structure for generate_single_mapping() from a more compact notation.

    Parameters
    ----------
    shorthand : str or list
    mode = option

    Returns
    -------
    aa_indices_in_cg : list of lists
        aa_indices_in_cg[i] returns list of aa indices in cg bead i
    cg_site_of_aa : list
        cg_site_of_aa[i] returns the cg site that atom i maps to, for us in generate_single_mapping

    Notes
    -----
    Shorthand versions:
    1) assuming AA atoms are in contiguous blocks for mapping
       [ naa1, naa2, ... ] each entry is # AA atoms in the next CG bead
    2) ... smiles? some other graph format? 
    3) I guess can technically wrap this up into generate_single_mapping
    '''
    if mode.lower() in ['simple','contiguous','short']:
        # expect data to be of the form [(cg_bead_name, n_aa_in_cg)]
        #cg_site_of_aa = [ ii for ii,num_in_cgbead in enumerate(shorthand) for jj in range(num_in_cgbead) ]
        nbead_so_far = 0
        aa_indices_in_cg = []
        cg_site_of_aa = []
        for cgindex,entry in enumerate(shorthand):
            if (isinstance(entry,list) or isinstance(entry,tuple)):
                print(entry)
                assert len(entry) == 2, "mapping entry for simple shorthand should be an integer or a 2-element tuple"
                cg_bead_name,num_in_cgbead = entry
            elif isinstance(entry,int):
                cg_bead_name,num_in_cgbead = 'DMY',entry 
            else:
                raise ValueError('Unrecognized mapping entry for simple shorthand: {}'.format(entry))
            indices = np.arange(nbead_so_far,nbead_so_far + num_in_cgbead).tolist()
            aa_indices_in_cg.append( (cg_bead_name,indices) ) 
            nbead_so_far += num_in_cgbead
            for ii in range(num_in_cgbead):
                cg_site_of_aa.append( (cg_bead_name,cgindex) )
        #cg_site_of_aa = []
        #aa_indices_in_cg = []
        #for ii,cgbead in enumerate(shorthand):
        #    aa_indices_in_cg.append( [ii]*cgbead )
        #    for jj in range( cgbead ):
        #        cg_site_of_aa.append(ii)  
    if mode.lower() in ['simplest','compactest','shortest']:
        # expect data to be of the form [(cg_bead_name, n_aa_in_cg, #cg bead)]
        nbead_so_far = 0
        ncg_so_far = 0
        aa_indices_in_cg = []
        cg_site_of_aa = []
        for cgindex,entry in enumerate(shorthand):
            if (isinstance(entry,list) or isinstance(entry,tuple)):
                if len(entry) == 3:
                    cg_bead_name,num_in_cgbead,num_cg_bead = entry
                elif len(entry) == 2:
                    cg_bead_name,num_in_cgbead,num_cg_bead = 'DMY',entry[0],entry[1]
            else:
                raise ValueError('Unrecognized mapping entry for simple shorthand: {}'.format(entry))
            for ii in range(num_cg_bead):
                indices = np.arange(nbead_so_far,nbead_so_far + num_in_cgbead).tolist()
                aa_indices_in_cg.append( (cg_bead_name,indices) ) 
                nbead_so_far += num_in_cgbead
                for jj in range(num_in_cgbead):
                    cg_site_of_aa.append( (cg_bead_name,ncg_so_far) )
                ncg_so_far += 1

    return aa_indices_in_cg, cg_site_of_aa


def generate_single_mapping(mode,extra=None):
    '''
    Generate mapping for a single chain

    Parameters
    ----------
    mode : str or None
        can be a [#atoms x 1] array to map from atom to CG site (i.e. cg site of atom i). Each entry should be tuple (cgbeadname,identifier)
    extra : several
        a mapping array of one of the accepted formats, or a topology for use with 121 or 1res21

    Return
    ------
    aa_indices_in_cg : list of lists
        aa_indices_in_cg[i] returns list of aa indices in cg bead i
    chain_indices_of_cg : list
        chain_id[i] returns the chain index of CG bead i. Hopefully is contiguous...
    '''
    if type(mode) is str:
        if isinstance(extra,mdtraj.Topology):
            top = extra
            if mode.lower() in ['121','1to1']:
                print('doing one-to-one mapping of AA bead to CG bead')
                cg_bead_names_len = top.n_atoms
                aa_indices_in_cg = [ (a.name,[a.index]) for a in top.atoms ]
                cg_site_of_aa = [ (a.name,a.index) for a in top.atoms ]
            elif mode.lower() in ['1res21','1resto1']:
                print('doing one residue to one CG bead')
                aa_indices_in_cg = [ (r.name,[a.index for a in r.atoms]) for r in top.residues ]
                cg_site_of_aa = [ (a.residue.name,a.residue.index) for a in top.atoms ]
        elif mode.lower() in ['cg_site_of_aa']:
            #if len(list(extra)) != top.n_atoms:
            #    raise ValueError('Fed in mapping array of CG site for each atom, should have dimensions #atoms x 1')
            print('using full custom mapping, extra input is [#atoms x 1] array specifying the CG site ID of each atom')
            cg_site_of_aa = extra
            cgsites = [] #ordered in sequence of occurrence
            aa_indices_in_cg = []
            for cgsite in extra: #this is only really efficient if a single molecule, since it doesn't assume contiguity.
                cgsite_identifier = cgsite[1]
                cgsite_name = cgsite[0]
                if cgsite_identifier not in cgsites:
                    cgsites.append(cgsite_identifier)
                    indices_in_cg = [ ii for ii,siteID in enumerate(extra) if siteID[1] == cgsite_identifier ]
                    aa_indices_in_cg.append( (cgsite_name,indices_in_cg) )
        elif mode.lower() in ['aa_indices_in_cg']:
            aa_indices_in_cg = extra
            num_aa = np.array([ len( entry[1] ) for entry in aa_indices_in_cg ]).sum()

            tmp = [[]] * num_aa
            aa_indices = np.zeros( num_aa )
            nbead_so_far = 0
            for isite, site_entry in enumerate(aa_indices_in_cg):
                for index in site_entry[1]:
                    tmp[nbead_so_far] = (site_entry[0],isite)
                    aa_indices[nbead_so_far] = index
                    nbead_so_far += 1
            #next, need to sort
            sorted_indices = np.argsort( aa_indices )
            cg_site_of_aa = [ tmp[index] for index in sorted_indices ]
             

        elif mode.lower() in ['simple','contiguous']:
            print('using compact mapping mode assuming contiguous sequencing of AA beads. Extra input should be array of # AA beads in each CG bead')
            aa_indices_in_cg, cg_site_of_aa = generate_single_mapping_shorthand('simple',extra)
            '''
            aa_indices_in_cg = []
            chain_index_of_cg = []
            nbead_so_far = 0
            for num_in_cgbead in extra:
                indices_in_cg = list(np.arange(nbead_so_far,nbead_so_far + num_in_cgbead)) 
                aa_indices_in_cg.append( indices_in_cg )
                nbead_so_far += num_in_cgbead
            '''

        elif mode.lower() in ['simplest','compactest']:
            print('using most compact mapping mode assuming contiguous sequencing of AA beads. Extra input should be array of # AA beads in each CG bead')
            aa_indices_in_cg, cg_site_of_aa = generate_single_mapping_shorthand('simplest',extra)

        #elif isinstance(mode,(np.ndarray,list)):
    else:
        raise ValueError('non-string mode not implemented yet')

    #sanity check. This function meant to be used on a single chain (molecule), so all the chain indices better be the same one!
    if isinstance(extra,mdtraj.Topology):
        chain_indices_of_cg = []
        for cg_index,entry in enumerate(aa_indices_in_cg):
            indices_in_cg = entry[1]
            chain_indices_in_cgbead = [ top.atom(ii).residue.chain.index for ii in indices_in_cg ] 
            if all( chain_index == chain_indices_in_cgbead[0] for chain_index in chain_indices_in_cgbead ):
                chain_indices_of_cg.append( chain_indices_in_cgbead[0] )
            else:
                 raise ValueError('chain indices {} of atoms in cg bead {} span more than one chain!'.format(chain_indices_in_cgbead,cg_index))
        
    #return aa_indices_in_cg, residue_indices_of_cg, chain_indices_of_cg
    return aa_indices_in_cg, cg_site_of_aa


def generate_system_mapping(mappings):
    ''' 
    Generate the AA->CG mapping for an entire system, from single-chain mappings. New atom indices are global atom indices instead of intrachain indices.

    Parameters
    ----------
    mapping: list of tuples
        each tuple should be (aa_index_in_cg_mapping_of_molecule, # occurrences of molecule)

    Return
    ------
    aa_indices_in_cg : list of lists
        aa_indices_in_cg[i] returns list of aa indices in cg bead i

    Notes:
    ------
    The sequence of mapping should correspond to that of the AA system!
    Assume each *molecule* is its own new chain.
    In the future, would be nice to have it s.t. just need to define mapping at residue level, knows how to handle patches, and can get the chain-level mapping. This will necessitate knowing how to manipulate the xml topology definitions.
    '''
    aa_indices_in_cg = []
    n_atoms_total = 0
    n_chains_total = 0
    
    for mapping in mappings:

        molecule_mapping = mapping[0]
        n_replicates = mapping[1] #of replicates
        n_beads_in_mapping = np.sum([ len(cgbead_mapping[1]) for cgbead_mapping in molecule_mapping ]) #i.e. # beads in chain/molecule

        cg_bead_names = [cgbead_mapping[0] for cgbead_mapping in molecule_mapping]
        for ii in range(n_replicates):
            new_indices = [ ( np.array(cgbead_mapping[1]) + n_atoms_total ).tolist() for cgbead_mapping in molecule_mapping ]

            aa_indices_in_cg.append( list(zip(cg_bead_names,new_indices)) )  #this way, can later iterate through by molecule/chain!
            n_atoms_total += n_beads_in_mapping

    system_aa_indices_in_cg = aa_indices_in_cg
    return system_aa_indices_in_cg


#def map(top,mapping_AA2CG,beadtypes):
    '''
    Parameters
    ----------
    top : mdtraj topology
    mapping : array
        following the output of my generate() helper functions, i.e. returns [ [aa_indices_in_cgbead1], [aa_indices_in_cgbead2], ... ], grouped by chain

    Return
    ------
    newtop : mdtraj topology

    Notes
    -----
    Current implementation best for mapping single chain/small systems. For larger systems, will have to do a parallel implementation

    Nice thing is that this code will figure out the bonds automatically!
    Only check for bonding if intrachain.
    '''
#    for ic in chains:
#        for ir in residues: #I don't use residues in CG model?
#            for ia in atomlistofresidue:
#                add new atom
#

# ===== Processing, Mapping and Synthesizing =====
def process_pdbfile(filename):
    '''
    For if the mapping uses the residue name column of pdbfile as the cg bead
    '''
    t = mdtraj.load(filename)
    cg_site_of_aa = [(a.residue.name,a.residue.index) for a in t.top.atoms]
    aa_indices_in_cg,_ = generate_single_mapping(mode='cg_site_of_aa',extra=cg_site_of_aa)

    return aa_indices_in_cg, cg_site_of_aa

def process_mappingfile(filename,mode,customfield=None):
    '''
    For if the (single chain) mapping uses a yaml/text-based specification of the array
    '''
    mapping = yaml.load(filename)
   
    if mode.lower() == 'pdb':
        my_mapping = mapping['pdbfile_mapping']
        aa_indices_in_cg, cg_site_of_aa = process_pdbfile( my_mapping )
        shorthand=None
    elif mode.lower() in ['aa','aa_indices_in_cg','aa_index','aa_indices']:
        my_mapping = mapping['aa_indices_in_cg']
        aa_indices_in_cg, cg_site_of_aa = generate_single_mapping( mode='aa_indices_in_cg',extra=my_mapping )
        shorthand=None
    elif mode.lower() in ['cg','cg_site_of_aa','cg_site','cg_sites']:
        my_mapping = mapping['cg_site_of_aa']
        aa_indices_in_cg, cg_site_of_aa = generate_single_mapping( mode='cg_site_of_aa',extra=my_mapping )
        shorthand=None
    elif mode.lower() in ['short','shorthand','simple','contiguous']:
        if customfield is None:
            my_mapping = mapping['shorthand']
        else:
            my_mapping = mapping[customfield]
        aa_indices_in_cg, cg_site_of_aa = generate_single_mapping_shorthand( mode='short',shorthand=my_mapping )
        shorthand = my_mapping
    elif mode.lower() in ['shortest','shortesthand','simplest']:
        if customfield is None:
            my_mapping = mapping['shorthand']
        else:
            my_mapping = mapping[customfield]
        aa_indices_in_cg, cg_site_of_aa = generate_single_mapping_shorthand( mode='shortest',shorthand=my_mapping)
        shorthand = my_mapping
    else:
        print('caution, trying to use a custom mode {} that may not be implemented yet').format(mode)
        if customfield is None:
            my_mapping = mapping[mode]
        else:
            my_mapping = mapping[customfield]
        aa_indices_in_cg, cg_site_of_aa, shorthand = generate_single_mapping( mode=mode,extra=my_mapping)
        #raise ValueError('specified mode has not been implemented yet')
       
    return aa_indices_in_cg, cg_site_of_aa, shorthand

def process_mapping_system(topdef):
    """
    Create a system_aa_indices_in_cg mapping with *global* atom indices, as well as a new cg-system topology.
    Modified such that can take a filename, or a system_spec list of lists!
    """
    if type(topdef) is str:
        topdef = yaml.load(topdef)
        system_spec = topdef['system']
    else: #treat topdef as dictionary of specifications
        system_spec = topdef['system']

    mappings = []
    mappingfiles = []
    tops = []
    for entry in system_spec:
        mappingfile = entry[0]
        mappingfiles.append(mappingfile)

        print('\n---> generating system mapping for {}'.format(mappingfile))
        if 'paths' in topdef: #use smart pathing to find file
          mappingfile_full = parsify.findpath(mappingfile, topdef['paths'])
          #print('updating {} to full path {}'.format(mappingfile,mappingfile_full))
          mappingfile = mappingfile_full

        n_replicates = entry[1]
        if len(entry) == 3:
            custom = entry[2]
        else:
            custom = None

        if mappingfile.endswith('pdb') and custom is None:
            print('Recognizing .pdb mapping specification for {}'.format(mappingfile))
            aa_indices_in_cg, cg_site_of_aa = process_pdbfile( mappingfile )
            traj_mapped = map_single( mdtraj.load(mappingfile),aa_indices_in_cg )
        elif mappingfile.endswith('pdb') and isinstance(custom,list):
            print('Assuming the pdb is the unmapped pdb, and the customlist is a shortesthand notation')
            aa_indices_in_cg,_ = generate_single_mapping_shorthand( mode='shortest',shorthand=custom )
            traj_mapped = map_single( mdtraj.load(mappingfile),aa_indices_in_cg )
        elif mappingfile.endswith('yaml'):
            chain_mapping_spec = yaml.load(mappingfile)
            if custom is None:
                print('Recognizing .yaml mapping specification for {}, using the aa_indices_in_cg section'.format(mappingfile))
                aa_indices_in_cg = chain_mapping_spec['aa_indices_in_cg'] 
                traj_mapped = mdtraj.load(parsify.findpath(chain_mapping_spec['pdbfile_mapped'], topdef['paths']))
            else:
                #(e.g. if using shorthand)
                print('Recognizing .yaml mapping specification for {}, using the {} section'.format(mappingfile,custom))
                aa_indices_in_cg,_,_ = process_mappingfile(mappingfile,mode=custom,customfield=custom)
                traj_mapped = map_single( mdtraj.load(parsify.findpath(chain_mapping_spec['pdbfile_unmapped'], topdef['paths'])),aa_indices_in_cg )
        else:
            raise ValueError('chain/molecule mapping format not recognized')
            
        mappings.append( (aa_indices_in_cg, n_replicates) )
        tops.append( (traj_mapped.top, n_replicates) )

    system_aa_indices_in_cg = generate_system_mapping(mappings)
    
    #create new topology
    new_topology = topologify.create_system(tops)
    bead_types = list(set( [ a.name for a in new_topology.atoms ] )) #not most efficient, can optimize later
    bead_types = [ [aname,1.0,0.0] for aname in bead_types ]

    #create necessary files for loading the mapped trajectory (compatible if # atoms > 99999)
    #ultimately just a system section where each line is [chain_pdb_file, # chains], don't try to infer CG bead types yet
    _trajs = []
    _nreplicates = []
    _chain_names = []
    _chain_files = []
    for ie,entry in enumerate(tops):
        top = entry[0]
        chain_name = os.path.splitext( os.path.basename(mappingfiles[ie]) )[0]
        chain_file = chain_name + '_mapped_top.pdb'
        _chain_files.append(chain_file)
        _trajs.append( mdtraj.Trajectory( np.zeros([1,top.n_atoms,3]), top ) )
        _nreplicates.append( entry[1] )
        _chain_names.append( chain_name ) 

        _trajs[-1].save( chain_file )
    #moltype_defs = dict( list(zip( _chain_names, _chain_files )) )
    moltype_defs = [ {'name':n, 'def':d} for (n,d) in zip(_chain_names,_chain_files) ]
    system_def = list(zip( _chain_names, _nreplicates ))

    mapped_def = yaml.YAML.comments.CommentedMap( [ ('bead_types',bead_types), ('res_types',None), ('mol_types',moltype_defs), ('system',system_def) ] )
    #if 'paths' in topdef:
    #  mapped_def['paths'] = topdef['paths']
    mapped_def['paths'] = os.getcwd() #note that in the above, the chain_files are newly generated single-chain mapped files, in case the given definition was not in .pdb format

    print('\n---> final mapped_def:')
    print(mapped_def)

    return system_aa_indices_in_cg, new_topology, mapped_def


def map_multiple(traj,cgtop,system_mapping):
    '''
    main difference from map_single is that these mappings are now organized by chain
    Notes:
    ------
    Technically... didn't need to generate the system_cg_mapping to use global atom indices, since can reindex the atoms in a chain s.t. can use the intrachain indexing again! But this way keeps things fairly simple.

    In the future, can either put in multiprocessing into this function, or wrap it
    '''

    num_cg = np.array([ len( entry ) for entry in system_mapping ]).sum()
    xyz = np.zeros( [traj.n_frames,cgtop.n_atoms,3] )
    num_cg_so_far = 0
    for ichain,chain in enumerate(system_mapping):
        if (ichain%100) == 0:
          print('mapping molecule {}'.format(ichain))
        for ibead,cgbead_entry in enumerate(chain):
            aa_indices = cgbead_entry[1]
            masses = np.array([traj.top.atom(ia).element.mass for ia in aa_indices])

            com = np.sum(traj.xyz[:,aa_indices,:] * masses[None,:,None],1) / masses.sum()
            xyz[:,num_cg_so_far,:] = com
            num_cg_so_far += 1
    print('finished mapping molecule {}'.format(ichain))
    new_traj = mdtraj.Trajectory(xyz,cgtop,unitcell_lengths = traj.unitcell_lengths, unitcell_angles = traj.unitcell_angles)
    return new_traj

def map_single(traj,mapping):
    '''
    take in a mapping array, and generate a new mapped topology for the chain, retaining bonds!
    Parameters
    ----------
    traj : mdtraj trajectory
    mapping : array
        should be aa_indices_in_cgbead array format, of a single chain
    '''
    cg_top = mdtraj.Topology()
    #=== first add all the species ===
    print('creating topology: adding species')
    ch = cg_top.add_chain()
    n_beads = len(mapping)
    xyz = np.zeros([n_beads,3])
    cg_beads = []
    for ic,cgbead_entry in enumerate(mapping):
        cgbead_name,aa_indices = cgbead_entry
        resname = traj.top.atom(aa_indices[0]).residue.name #use first atom's residue as the resname. alternative: use most common occurrence
        masses = np.array([traj.top.atom(ia).element.mass for ia in aa_indices])
        com = np.sum(traj.xyz[0,aa_indices,:] * masses[:,None],0) / masses.sum()
        #print( mdtraj.compute_center_of_mass(traj.atom_slice(aa_indices))[0] )
        #print( com )
        xyz[ic,:] = com

        res = cg_top.add_residue( resname,ch )
        element = mdtraj.element.virtual #for now... later on can create elements for CG simulation from cg bead name
        bead = cg_top.add_atom(cgbead_name, element, res)
        cg_beads.append(bead)

    #=== figure out connectivity ===
    #loop over each cg bead, then over other cg beads, then check bonding
    print('creating topology: figuring out bonding')
    bondgraph = traj.top.to_bondgraph()
    for ii in range(n_beads):
        aa_indices1 = mapping[ii][1]
        for jj in range(ii+1,n_beads):
            aa_indices2 = mapping[jj][1]
            bonded = any([ bondgraph.has_edge(traj.top.atom(a1),traj.top.atom(a2)) for a1 in aa_indices1 for a2 in aa_indices2 ])
            if bonded:
                cg_top.add_bond(cg_beads[ii],cg_beads[jj])

    #=== create mapped trajectory and save ===
    cg_traj = mdtraj.Trajectory([xyz],cg_top)

    return cg_traj

# ===== Main =====
if __name__ == "__main__":
    ''' example calls
    python mapper.py -traj pS20.pdb -prefix test -style single -params pS20_mapping.yaml -mode cg
    python mapper.py -traj pS20.pdb -prefix test -style single -params pS20_mapping.yaml
    python mapper.py -traj pS20.pdb -prefix test -style single -params pS20_mapping.pdb

    python mapper.py -traj pS20-2-SDS-3.pdb -prefix test_system -style system -params system_map.yaml
    python mapper.py -traj pS20-2-SDS-3.pdb -prefix test_system -style system -params pS20_mapping.pdb 2 SDS_mapping.pdb 3
    python mapper.py -traj pS20-2-SDS-3.pdb -prefix test_system -style system -params pS20_mapping.pdb 2 SDS.yaml 3
    '''
    import argparse as ap
    parser = ap.ArgumentParser(description='general mapper utility')
    parser.add_argument('-traj', type=str, required=True, help='traj to map')
    parser.add_argument('-top', type=str, nargs='+', default=None, help='topology for reading non-pdb files')
    parser.add_argument('-prefix', type=str, default=None, help='optional prefix. default is to use traj filename root.')
    parser.add_argument('-style', type=str, required=True, choices=['single','system'], help='whether to use single-chain or system parsing')
    parser.add_argument('-mode', type=str, choices=['aa','cg','pdb','short','shortest'], default='aa', help='mapping specification format')
    parser.add_argument('-params', nargs='+', help = 'files specifying mappings. If more than one, can collate together into a system.')
    parser.add_argument('-stride', default=1, type=int, help='stride for processing trajectory')
    parser.add_argument('-skip', default=0, type=int, help='stride for processing trajectory')
    args = parser.parse_args()

    print(args.top)

    if args.prefix is None:
        prefix = os.path.splitext(args.traj)
    else:
        prefix = args.prefix

    if args.style == 'single':
        if len( args.params ) > 1:
            raise ValueError('single mode selected, but #params more than one: {}'.format(args.params))
        t = mdtraj.load(args.traj)
        filemap = args.params[0]
        mode = args.mode
        
        if filemap.endswith('pdb'):
            aa_indices_in_cg, cg_site_of_aa = process_pdbfile( filemap )
            shorthand = None
        else:
            aa_indices_in_cg, cg_site_of_aa, shorthand = process_mappingfile(filemap,mode)
        
        traj_mapped = map_single(t,aa_indices_in_cg)
        traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
        save(prefix, aa_indices_in_cg, cg_site_of_aa, 
                traj_mapped = traj_mapped, traj_mapping = traj_mapping, 
                shorthand = shorthand, unmapped = args.traj )
    elif args.style == 'system':
        #Note for simplicity the 'system' style does not take a mode specification. Only uses .pdb mappings or aa_indices_in_cg mappings.
        if args.top is None:
            t = mdtraj.load(args.traj)
            #t = mdtraj.load(args.traj, stride=args.stride)
        else:
            if isinstance(args.top,list) and len(args.top)==2:
              import pandas
              df = pandas.read_csv(args.top[0])
              bonds = np.loadtxt(args.top[1])
              top_mdtraj = mdtraj.Topology.from_dataframe(df,bonds)
              top = top_mdtraj
            else:
              top = args.top[0]
            t = mdtraj.load( args.traj, top=top )
            #t = mdtraj.load( args.traj, top=top, stride=args.stride )
        t = t[args.skip::args.stride]

        if len( args.params ) == 1:
            print('1 system parameter file received, assume contains system information')
            filemap = args.params[0]
            system_spec = yaml.load(filemap)
            #tmp = yaml.load(filemap)
            #system_spec = tmp['system']
            #if 'paths' in tmp:
              #for ie,entry in enumerate(system_spec):
              #  fname = entry[0]
              #  fname_full = parsify.findpath(fname, tmp['paths'])
              #  system_spec[ie][0] = fname_full
              #  print('updating {} to full path {}'.format(fname,fname_full))
        elif len( args.params ) >= 1:
            print('multiple parameters received, assume is pairs of: mapping_file, # of molecule')
            if len(args.params) % 2 != 0:
                raise ValueError('need even number of parameter arguments, but received odd amount')
            #now, essentially process the arguments into a system specification that I would've expected in a system mapping file        
            chain_mapping_files = args.params[::2]
            chain_numbers = [ int(n) for n in args.params[1::2] ]
            system_spec = list( zip( chain_mapping_files, chain_numbers ) )

        #do the mapping
        print('=== Generate mapping for the system ===')
        system_aa_indices_in_cg, new_topology, mapped_def = process_mapping_system( system_spec )
        print('\n=== Do the mapping ===')
        new_traj = map_multiple(t, new_topology, system_aa_indices_in_cg)
        print('\n=== Save ===')
        new_traj.save(prefix + '_mapped.dcd')
        new_traj[0].save(prefix + '_mapped.pdb')
        mapped_def['pdb'] = prefix + '_mapped.pdb'
        
        #save using prefix
        d = { 'traj_unmapped': args.traj,
              'top_unmapped': args.top,
              'traj_mapped': prefix + '_mapped.dcd',
              'top_mapped': prefix + '_mapped.pdb',
              'system': system_spec,
              'system_aa_indices_in_cg': system_aa_indices_in_cg
            }
        yaml.save_dict( prefix + '_mapping.yaml', d, header = '{} mapping summary'.format(prefix) ) 
        yaml.save_dict( prefix + '_mapped.yaml', mapped_def, header = '{} mapped system definition'.format(prefix) )
    else:
        raise ValueError('Unrecognized style {}'.format(args.style))



# ===== Tests =====
## Testing System
'''
systemfile = 'system_map.yaml'
filename = 'pS20-2-SDS-3.pdb'
prefix = os.path.splitext(filename)[0]
t = mdtraj.load(filename)
system_aa_indices_in_cg, new_topology = process_mapping_system(systemfile)
new_traj = map_multiple(t, new_topology, system_aa_indices_in_cg)
new_traj.save(prefix + '_mapped.pdb')
'''

## Testing Single Chain
test_singlechain = False
'''
print('--- maping pS ---')
filename = 'pS20.pdb'
filemap = 'pS20.yaml'
t = mdtraj.load(filename)
prefix,mode = 'pS20','shortest'
prefix,mode = 'pS20','pdb'
aa_indices_in_cg, cg_site_of_aa, shorthand = process_mappingfile(filemap,mode)
traj_mapped = map_single(t,aa_indices_in_cg)
traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping, shorthand=shorthand)
'''

if test_singlechain:
    #t = mdtraj.load('SDS.pdb')
    #try to figure out the mapping
    #rando = [(a.residue.name,a.residue.index) for a in t.top.atoms]
    #print('start test')
    #print(rando)
    #chain_mapping,cg_indices = generate_single_mapping(t.topology,mode='full',extra=rando)

    #t = mdtraj.load('pS20.pdb')
    #rando = np.random.randint(0,5,160)
    #generate_system_mapping( [(chain_mapping,5)] )
    #chain_mapping2,cg_indices2 = generate_single_mapping(t.topology,mode='simple',extra=[8]*t.n_residues)


    # Experiment 2: map from an annotated pdb file
    '''
    filename = 'SDS_ordered_2gro.pdb'
    t = mdtraj.load(filename)
    prefix = os.path.splitext(filename)[0]
    aa_indices_in_cg, cg_site_of_aa = process_pdbfile('SDS.pdb')
    traj_mapped = map_single(t,aa_indices_in_cg)
    traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
    save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping )
    '''

    # Experiment 3: map from a given mapping array
    filename = 'SDS_ordered_2gro.pdb'
    filemap = 'SDS.yaml'
    t = mdtraj.load(filename)

    ## Test pdbfile
    print('--- test0 ---')
    prefix,mode = 't0','pdb'
    aa_indices_in_cg, cg_site_of_aa, _ = process_mappingfile(filemap,mode)
    traj_mapped = map_single(t,aa_indices_in_cg)
    traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
    save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping )

    ## Test aa_indices_in_cg
    print('--- test1 ---')
    prefix,mode = 't1','aa'
    aa_indices_in_cg, cg_site_of_aa, _ = process_mappingfile(filemap,mode)
    traj_mapped = map_single(t,aa_indices_in_cg)
    traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
    save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping )

    ## Test cg_site_of_aa
    print('--- test2 ---')
    prefix,mode = 't2','cg'
    aa_indices_in_cg, cg_site_of_aa, _ = process_mappingfile(filemap,mode)
    traj_mapped = map_single(t,aa_indices_in_cg)
    traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
    save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping )

    ## Test shorthand
    print('--- test3 ---')
    prefix,mode = 't3','short'
    aa_indices_in_cg, cg_site_of_aa, shorthand = process_mappingfile(filemap,mode)
    traj_mapped = map_single(t,aa_indices_in_cg)
    traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
    save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping, shorthand=shorthand )

    ## Test shortesthand
    print('--- test4 ---')
    prefix,mode = 't4','shortest'
    #aa_indices_in_cg, cg_site_of_aa = process_mappingfile(filemap,mode)
    aa_indices_in_cg, cg_site_of_aa, shorthand = process_mappingfile(filemap,mode,customfield='shortesthand')
    traj_mapped = map_single(t,aa_indices_in_cg)
    traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
    save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping, shorthand=shorthand)

