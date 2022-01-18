# Helper functions for setting up a polyfts run from my custom "universal" data structures
# In this instance, would be writing out a polyfts input parameter file!
# (C) Kevin Shen, 2021
# will almost look like... setting up shorthand print notation for a lot of these things
# looks like field = values as way to enter things
#
# TODO:
# 1) conversions (units? or between file formats)
# 6) running commands
# 7) simple subsitutions of a script
# 8) parse PFTS script into a dictionary?
#     sequence of # --> list
#     = --> :
#     \n --> , (unless is last element of that dict)
#     ignore comments
#    can probably do line by line, "increment/decrement" dictionary level
# 
# Example:
# - with a Sys constructed from ffdef, topdef, settings
# - save/update the ffdef with Sys: 
#   export_sim.save_params_to_dict( Sys.ForceField, ffdef.processed_file )
#   ffdef.sim2md()
#   can also save directly without the Sys object, via
#     export_sim.update_ffdef_from_paramstring(ffdef,ffstringfile)

# - write out
#   ffdef.save(savename)
#
# - can then proceed to set up a pfts system!
#

#import sys,os
#import re,copy
#import context #makes sure we get the forcefield.py and optimize.py defined in utility instead of legacy versions in parent directory
import topologify,forcefield
from collections import OrderedDict
import yamlhelper as yaml
import networkx as nx
import numpy as np
from pfts import *

## Helpers
def blockify(sequence):
  """
  return "shorthand" that is (type, #) in a linear sequence
  i.e. AAABBAA --> (A,3), (B,2), (A,2)
  for parsing chain sequence

  Args:
    sequence: iterable, has len() method
  """
  n = len(sequence)
  blocks = []
  for ii,element in enumerate(sequence):
    if ii == 0:
      current = element
      count = 1
    elif element == current:
      count += 1
    elif element != current:
      blocks.append( (current,count) )
      current = element
      count = 1
  blocks.append((current,count))
  return blocks

def validate_moltype(bead_types,beads,bonds,name=None):
  """
  Make sure that a given molecule definition is PolyFTS-representable
  Args:
    bead_types: dictionary of bead_types
    beads: list of bead names in the chain
    bonds: list of bond tuples
  """
  nbonds = len(bonds)
  nbeads = len(beads)

  if nbeads == 1:
    return {'style':'smallmol','def':[(beads[0],1)]}

  mol_graph = nx.Graph()
  mol_graph.add_edges_from(bonds)
  
  polydef = None

  # check 1: are all bead types defined?
  for b in beads:
    if b not in bead_types:
      raise ValueError('bead type {} in molecule {} not defined'.format(b,name))
  # check 2: architecture. only check for linear right now, can extend later:w
  if not nx.is_connected(mol_graph):
    if nbonds < nbeads - 1:
      raise ValueError('too few bonds ({0}) defined for this {2} {1}-mer'.format(nbonds,nbeads,name))
    else: 
      raise ValueError('molecule {} is missing bonds, is not fully connected'.format(name))
  elif nx.is_tree(mol_graph):
    #check if is a linear
    n1 = [node for (node,val) in mol_graph.degree() if val == 1] #nodes with degree 1, i.e. ends
    if len(n1) == 2: #is linear, no branches
      path = nx.shortest_path(mol_graph, n1[0],n1[1]) 
      beads_ordered = [ beads[index] for index in path ]
      polydef = {'style':'linear', 'def':blockify(beads_ordered)}
    else:
      #TODO: check if is a feasible bottlebrush
      raise NotImplementedError('molecule {} is not unbranched, is not supported yet'.format(name))
  else:
    raise NotImplemented('molecule {} is not a tree, not supported yet'.format(name))

  return polydef



## Set up force field
def create_models(top,ff,settings):
  """
  Create the dictionary representation that will be easiest to write out into PFTS format.
  Initially, only use one model?

  PFTS is case-insensitive, but for consistency here I'll enforce CamelCase
  """
  # Load
  if isinstance(top,str):
    top = topologify.Topology(top)
  if isinstance(ff,str):
    ff = forcefield.ForceField(ff)
  if isinstance(settings,str):
    settings = yaml.load(settings)

  # === Set up skeleton ===
  MODELS = OrderedDict()
  MODELS['NumModels'] = 1
  MODELS['ModelType'] = 'PolymerAF'
  MODELS['Monomers'] = OrderedDict()
  MODELS['Chains'] = OrderedDict()
  MODELS['SmallMolecules'] = OrderedDict()
  MODELS['Models'] = []

  Simulation = OrderedDict()
  PLL = OrderedDict()

  # === Set up ===
  result = preprocess_topology(top)
  MODELS['Monomers'] = result['beads']
 
  # TODO
  # 0) units
  # 1) monomers: how to add polarizability, Kuhn length, etc? get from topology file?
  # 3) interactions


def preprocess_topology(top):
  """
  Pre-process topology into a easily polyfts-translatable representation.
  """
  # === Figure out topology part first ===
  # work from the fully verbose format
  if isinstance(top,str):
    top = topologify.Topology(top)
  
  # get unique types
  unique_beads = set()
  mol_types = OrderedDict() 
  bead_types = OrderedDict()
  bead_in_sys = OrderedDict()
  mol_in_sys = OrderedDict()

  for bead_def in top.processed_file['bead_types']:
    #print(bead_def)
    if 'fields' in bead_def or 'default' in bead_def: continue
    bead_name = bead_def['name']
    bead_charge = bead_def['charge']
    if bead_name not in bead_types:
      bead_types[bead_name] = OrderedDict([('charge',bead_charge)])
    else:
      raise ValueError('bead type {} has conflicting name/prior definition'.format(bead_name))
  
  for mol_def in top.processed_file['mol_types']:
    mol_name = mol_def['name']
    beads = mol_def['beads']
    bonds = mol_def['bonds']
    pftsdef = validate_moltype(bead_types, beads, bonds, mol_name)

    if mol_name not in mol_types:
      mol_types[mol_name] = OrderedDict([('beads',beads),('bonds',bonds)])
      mol_types[mol_name].update(pftsdef)
    else:
      raise ValueError('mol type {} has conflicting name/prior definition'.format(mol_name))
   
  for mol_name,mol_num in top.processed_file['system']:
    if mol_name not in mol_types:
      raise ValueError('mol type {} in system is undefined'.format(mol_name))
    if mol_name not in mol_in_sys:
      mol_in_sys[mol_name] = mol_types[mol_name]
      mol_in_sys[mol_name]['num'] = mol_num
    else:
      mol_in_sys[mol_name]['num'] += mol_num

    unique_beads |= set( mol_types[mol_name]['beads'] )

  for bead_name,bead_def in bead_types.items():
    if bead_name in unique_beads:
      bead_in_sys[bead_name] = bead_def
  

  result = { 'beads': bead_in_sys, 'molecules':mol_in_sys } #only holds the types that are in the system, i.e. prune out unused types
  return result


def generate_input(topdef, ffdef, settings):
  '''
  set up the dictionaries and formats for PFTS.
  TODO: consider breaking this function up into smaller pieces. but some data needs to be shared across multiple sections...
  '''
  ### ===== Parsing and organizing =====
  print('\n===== LOADING DEFINITIONS =====')
  if isinstance(topdef,str):
    topdef = topologify.Topology(topdef)
  if isinstance(ffdef,str):
    ffdef = forcefield.ForceField(ffdef)
  if isinstance(settings,str) and settings.endswith('yaml'):
    settings = yaml.load(settings)
  elif isinstance(settings,str) and settings.endswith('in'): #pfts format
    settings = load(settings)


  processed_top = preprocess_topology(topdef) 
  n_bead_types = len(processed_top['beads'])
  beadID = OrderedDict( [(name,(ic+1)) for ic,(name,beaddef) in enumerate( processed_top['beads'].items() )] )

  # convention for FTS units, as well as sim-units!
  Rgref = 1.0 #what ever length unit is used in sim
  Nref  = 1
  bref  = Rgref * 6.0**0.5/float(Nref)

  ### ===== Start saving =====
  spec = OrderedDict()
  spec['InputFileVersion'] = 3
  
  ### ===== Models =====
  spec['Models'] = OrderedDict()
  tmp = spec['Models']
  tmp['NumModels'] = 1
  tmp['ModelType'] = 'POLYMERAF'
  tmp['Monomers'] = OrderedDict()
  tmp['Monomers']['NSpecies'] = n_bead_types
  charges = [ v['charge'] for v in processed_top['beads'].values() ]
  tmp['Monomers']['Charge'] = charges
  #Kuhn Len TODO
  #GaussSmearWidth TODO

  ## === chains and molecules ===
  print('\n===== DEFINING CHAINS AND MOLECULES =====')
  chains = OrderedDict()
  smallmols = OrderedDict()
  chain_beadnum = [] #total # of beads contributed by a chain to the system
  smallmol_beadnum = [] #total # of beads contributed by a small molecule to the system

  for mol_name,mol_def in processed_top['molecules'].items():
    if len(mol_def['beads']) == 1: #small molecule
      smallmols[mol_name] = mol_def 
    else:
      chains[mol_name] = mol_def
  
  tmp['Chains'] = OrderedDict([('NChains',len(chains)), ('Contourds',1), ('DiffuserMethod','SOS'), ('PolymerReferenceN',Nref)])
  for ic,(mol_name,mol_def) in enumerate(chains.items()):
    chainid = 'Chain{}'.format(ic+1)
    nblocks = len(mol_def['def'])
    nperblock = []
    blockspecies = []
    for block in mol_def['def']:
      blockspecies.append(beadID[block[0]])
      nperblock.append(block[1]) 
    nbeads = sum(nperblock)
    if nbeads != len(mol_def['beads']): raise ValueError("inferred bead number {} for molecule {} doesn't match explicit bead sequence".format(nbeads,mol_name))
    chain_beadnum.append( mol_def['num'] * nbeads )

    chain_pftsdef = [ ('Label',mol_name),
                      ('Architecture',mol_def['style']),
                      ('Statistics','TODO'),
                      ('NBlocks', nblocks),
                      ('NBeads', nbeads),
                      ('BlockSpecies',blockspecies),
                      ('NPerBlock',nperblock)
                    ]
    tmp['Chains'][chainid] = OrderedDict( chain_pftsdef )

  tmp['SmallMolecules'] = OrderedDict({'NSmallMoleculeTypes':len(smallmols)})
  for ic,(mol_name,mol_def) in enumerate(smallmols.items()):
    smallmolid = 'SmallMolecule{}'.format(ic+1)
    beadtype = mol_def['beads'][0]
    tmp['SmallMolecules'][smallmolid] = {'Species':beadID[beadtype]}
    smallmol_beadnum.append( mol_def['num'] )
  # final calculations
  total_beadnum = sum(chain_beadnum) + sum(smallmol_beadnum)
  chain_beadfrac = [ float(_n)/float(total_beadnum) for _n in chain_beadnum ]
  smallmol_beadfrac = [ float(_n)/float(total_beadnum) for _n in smallmol_beadnum ]
  #precision = 10
  #chain_beadfrac = [ float(np.round(_n/total_beadnum,decimals=precision+1)) for _n in chain_beadnum ]
  #smallmol_beadfrac = [ float(np.round(_n/total_beadnum,decimals=precision+1)) for _n in smallmol_beadnum ]

  ## === model ===
  tmp['Model1'] = OrderedDict()
  #  = cell. mandatory: CellLengths, NPW = 
  print('\n===== DEFINING CELL =====')
  tmp['Model1']['Cell'] = OrderedDict()
  box = parse2list(settings['Cell']['CellLengths'])
  settings['Cell']['CellLenths'] = box
  box_dim = len(box)
  npw = parse2list(settings['Cell']['NPW'])
  if len(npw) != box_dim:
    raise ValueError('dimensions of NPW {} incommensurate with box {}'.format(npw,box))
  if 'CellAngles' in settings['Cell']:
    cell_angles = parse2list(settings['Cell']['CellAngles'])
    settings['Cell']['CellAngles'] = cell_angles
    #if len(cell_angles) != box_dim:
    #  raise ValueError('Cell Angles {} incommensurate with box {}'.format(cell_angles,box))
  else:
    cell_angles = [90.] * box_dim
    settings['Cell']['CellAngles'] = cell_angles

  update(tmp['Model1']['Cell'],settings['Cell'],'Dim',default=box_dim)
  update(tmp['Model1']['Cell'],settings['Cell'],'CellScaling',default=1.0)
  update(tmp['Model1']['Cell'],settings['Cell'],'CellLengths',mandatory=True)
  update(tmp['Model1']['Cell'],settings['Cell'],'CellAngles',default=cell_angles)
  update(tmp['Model1']['Cell'],settings['Cell'],'NPW',mandatory=True)
  update(tmp['Model1']['Cell'],settings['Cell'],'SpaceGroupName',fill = False)
  update(tmp['Model1']['Cell'],settings['Cell'],'CenterToPrimitiveCell',
          fill = 'SpaceGroupName' in tmp['Model1']['Cell'], default = True)
  update(tmp['Model1']['Cell'],settings['Cell'],'Symmetrize',
          fill = 'SpaceGroupName' in tmp['Model1']['Cell'], default = 'off')

  # = interactions =
  # currently only supports, chi, electrostatics, BexclVol. Doesn't do a whole lot of validation!
  #allow for inputting Kuhn Lengths in normal length, instead of in relative Kuhn unit b = Rgref*sqrt(6)/Nref
  print('\n===== DEFINING INTERACTIONS =====')
  tmp['Model1']['Interactions'] = OrderedDict()
  if ffdef is None: #use standard definition format from settings.yaml
    update(tmp['Model1']['Interactions'],settings['Interactions'],'ApplyCompressibilityConstraint',default=True)
    update(tmp['Model1']['Interactions'],settings['Interactions'],'Compressibility_InvZetaN',
            fill = not tmp['Model1']['Interactions']['ApplyCompressibilityConstraint'], default=100.0)
    update(tmp['Model1']['Interactions'],settings['Interactions'],'EElecStatic', fill = False)
    for k,v in settings['Interactions'].items():
      if k.lower().startswith('chi') or k.lower().startswith('bexclvolume'):
        tmp['Model1']['Interactions'][k] = v

    MonPatches = settings['Interactions'].get('MonPatches',{})
    MolPatches = settings['Interactions'].get('MolPatches',{})
      
    kuhns = [1.0]*n_bead_types
    asmear = [None]*n_bead_types
    for beadname,beadparams in MonPatches.items():
      if beadname in beadID:
        kuhns[beadID[beadname]-1] = beadparams['Kuhn']
        asmear[beadID[beadname]-1] = beadparams.get('Smearing',None)

  else: #currently only support BExcl Volume
    #ugh messy... in general can have FJC bonds, Kuhn bonds, and they might not be consistent with one another! 
    #need to make some reducing assumptions, e.g. all bead types are the same

    print('Given force field file, ignoring the Interactions section of settings file')
    print('Turning off (In)compressibility constraint')
    tmp['Model1']['Interactions']['ApplyCompressibilityConstraint'] = False
    MonPatches = OrderedDict([(k,{}) for k in beadID])
    MolPatches = OrderedDict([(k,{}) for k in chains])
    #aggregate BExclVol Matrix, asmear
    print(ffdef.processed_file['pair_ljg']['params_sim'])
    if 'pair_ljg' in ffdef.processed_file and 'params_sim' in ffdef.processed_file['pair_ljg']:
      print('\n--> Detected LJ-Gaussian definition for BExcludedVolume')
      Bsim = np.zeros([n_bead_types,n_bead_types])
      Kappa = np.ones([n_bead_types,n_bead_types])
      num_defs = np.zeros([n_bead_types,n_bead_types])

      for p in ffdef.processed_file['pair_ljg']['params_sim']:
        #print(p)
        num_new_defs = np.zeros([n_bead_types,n_bead_types])
        for sp1 in p['species'][0]:
          if sp1 not in beadID:
            continue
          ind1 = beadID[sp1]-1
          for sp2 in p['species'][1]:
            if sp2 not in beadID:
              continue
            ind2 = beadID[sp2]-1
            if p['Epsilon']['val'] != 0:
              raise NotImplementedError('... non-zero LJ epsilon detected in potential {}, not supported'.format(p))
            if p['Dist0']['val'] != 0:
              raise NotImplementedError('... non-zero Gaussian Dist0 detected in potential {}, not supported'.format(p))
            #store and symmetrize
            Bsim[ind1,ind2] = p['B']['val']
            Bsim[ind2,ind1] = p['B']['val']
            Kappa[ind1,ind2] = p['Kappa']['val']
            Kappa[ind2,ind1] = p['Kappa']['val']
            num_new_defs[ind1,ind2] = 1
            num_new_defs[ind2,ind1] = 1
        num_defs += num_new_defs

      print('... number of LJG definitions for each pair:\n{}'.format(num_defs))
      print('... bead IDs: {}'.format(beadID))
      if np.any( num_defs > 1 ):
        #print('some nonbonded definitions were defined more than once, ambiguous')
        raise ValueError('... some nonbonded definitions were defined more than once, ambiguous')
      print('... Detected B matrix (in sim units):\n{}'.format(Bsim))
      print('... Detected Kappa matrix (in sim units):\n{}'.format(Kappa))
 
      #do unit conversions
      abar2= 0.5/Kappa
      u0   = Bsim * (2*np.pi*abar2)**1.5 
      Bfts = u0 * Nref**2.0/Rgref**3.0
      print('... Processed BExclVol matrix (in pfts units):\n{}'.format(Bfts))

      #figure out asmear
      asmear = np.diag( np.sqrt(abar2/(2.0)) )
      abar2estimate = asmear[:,None]**2.0 + asmear[None,:]**2.0
      abar2estimate[Bsim==0.] = abar2[Bsim==0.] #if Bsim=0, irrelevant whether abar2 matches (i.e. may not have been defined)
      print('... abar^2 = ai^2 + aj^2 from Kappa:\n{}'.format(abar2))
      print('... abar^2 inferred from diagonal components of Kappa: (ignoring Bsim=0.0 entries)\n{}'.format(abar2estimate))
      
      if not np.isclose(abar2,abar2estimate).all():
        raise ValueError('... Kappas are not consistent with smearing model, ambiguous how to resolve')
      else:
        print('... Great! Kappas consistent with smearing length model.')
        asmear = asmear.tolist()
      for mon_name,bead_id in beadID.items():
        MonPatches[mon_name]['Smearing'] = asmear[bead_id-1]

      #save interaction matrix
      for ind1 in range(n_bead_types):
        for ind2 in range(ind1,n_bead_types):
          bexclvname = 'BExclVolume{}{}'.format(ind1+1,ind2+1)
          tmp['Model1']['Interactions'][bexclvname] = float(Bfts[ind1,ind2])

    else:
      raise ValueError('... No excluded volume interactions defined, are you sure?')

    #detect electrostatics
    if 'coulomb_smeared' in ffdef.processed_file and 'params_sim' in ffdef.processed_file['coulomb_smeared']:
      print('\n--> Detected Electrostatics')
      Coefs = []
      for p in ffdef.processed_file['coulomb_smeared']['params_sim']:
        Coefs.append(p['Coef']['val'])
      Coefs = np.array(Coefs)
      coefs_same = np.isclose(Coefs,Coefs[0])
      if all(coefs_same):
        tmp['Model1']['Interactions']['EElecStatic'] = float(Coefs[0]*4.0*np.pi * Nref**2.0 / Rgref)
        print('... Coef {} --> EElecStatic {}'.format(Coefs[0], tmp['Model1']['Interactions']['EElecStatic']))
      else:
        raise ValueError('... charge coefficients {} in ffdef are not all the same!'.format(Coefs))
    #TODO: technically, need to check that smearing is consistent with Excluded Volume above...

    #figure out bonding
    #do some inference... if Dist0=0, treat as DGC
    #if Dist0!=0 and FConst>=100.0, treat as FJC
    kuhns = [1.0, 1.0, 1.0]
    if 'bond_harmonic' in ffdef.processed_file and 'params_sim' in ffdef.processed_file['bond_harmonic']:
      print('\n--> Detected Bonding')
      FConst = np.zeros([n_bead_types,n_bead_types])
      Dist0 = np.zeros([n_bead_types,n_bead_types])
      num_defs = np.zeros([n_bead_types,n_bead_types])
      isDGC = np.full([n_bead_types,n_bead_types],True).astype(bool)
      bondl = np.zeros([n_bead_types,n_bead_types])

      for p in ffdef.processed_file['bond_harmonic']['params_sim']:
        print(p)
        num_new_defs = np.zeros([n_bead_types,n_bead_types])
        for sp1 in p['species'][0]:
          ind1 = beadID[sp1]-1
          for sp2 in p['species'][1]:
            ind2 = beadID[sp2]-1
            #store and symmetrize
            FConst[ind1,ind2] = p['FConst']['val']
            FConst[ind2,ind1] = p['FConst']['val']
            Dist0[ind1,ind2] = p['Dist0']['val']
            Dist0[ind2,ind1] = p['Dist0']['val']
            num_new_defs[ind1,ind2] = 1
            num_new_defs[ind2,ind1] = 1
        num_defs += num_new_defs
      if np.any( num_defs > 1 ):
        raise ValueError('... some bonded definitions were defined more than once, ambiguous how to resolve')
      print('... FConst:\n{}'.format(FConst))
      print('... Dist0:\n{}'.format(Dist0))
      for ind1 in range(n_bead_types):
        for ind2 in range(n_bead_types):
          if Dist0[ind1,ind2] == 0.0: #assume DGC, FConst = 3/2b^2 --> b = sqrt(1.5/FConst)
            isDGC[ind1,ind2] = True
            bondl[ind1,ind2] = np.sqrt(1.5/FConst[ind1,ind2]) if FConst[ind1,ind2] > 0.0 else 0.0
          elif FConst[ind1,ind2] >= 100.0: #assume FJC
            isDGC[ind1,ind2] = False
            bondl[ind1,ind2] = Dist0[ind1,ind2]
          else:
            print('CAUTION, Dist0!=0, but FConst is also not really large enough to maintain bondl=Dist0')
            isDGC[ind1,ind2] = False
            bondl[ind1,ind2] = Dist0[ind1,ind2]
      #try to impute missing (i,i) bonds (PolyFTS currently only supports (i,i) bonds    
      #TODO: should try to do this more intelligently... but may just have to do this manually afterwards
      print('... Read in inferred bond lengths:\n{}'.format(bondl))
      for ind1 in range(n_bead_types):
        if bondl[ind1,ind1] == 0.0:
          print('... Imputing bond length for missing ({},{}) type'.format(ind1,ind1))
          bondl[ind1,ind1] = np.max(bondl[ind1,:]) #currently just use the largest junction bond
      kuhns = np.diag(bondl)
      print('... Final estimated intrinsic bondlengths: {}'.format(kuhns))
    elif len(chains) > 0:
      raise ValueError('... chains defined but no bonding defined')
 
    for chain_name,chain_def in chains.items(): #need to infer chain statistics
      #count bond types
      print('... estimating statistics for chain {}'.format(chain_name))
      _beads = chain_def['beads']
      _bondcount = {}
      for bond in chain_def['bonds']:
        typetuple = ( _beads[bond[0]],_beads[bond[1]] )
        typetuple_reverse = ( _beads[bond[1]],_beads[bond[0]] )
        if typetuple not in _bondcount and typetuple_reverse not in _bondcount:
          _bondcount[typetuple] = 1
        elif typetuple in _bondcount:
          _bondcount[typetuple] += 1
        elif typetuple_reverse in _bondcount:
          _bondcount[typetuple] += 1
        else:
          raise ValueError('... bond {} of type {} not counted properly'.format(bond,typetuple))
      _bondcounts = [ v for v in _bondcount.values() ]
      if sum(_bondcounts) != len(chain_def['bonds']):
        raise ValueError('... bond count {} wrong, expecting {}'.format(sum(_bondcounts),len(chain_def['bonds'])))
      #now check whether is chain smore FJC or more DGC
      _numDGCbonds = 0
      _numFJCbonds = 0
      for bondtype,bondnum in _bondcount.items():
        ind1 = beadID[bondtype[0]]-1
        ind2 = beadID[bondtype[1]]-1
        if isDGC[ind1,ind2]:
          _numDGCbonds += bondnum
        else:
          _numFJCbonds += bondnum
      print('...... estimated {} DGC bonds, {} FJC bonds'.format(_numDGCbonds,_numFJCbonds))
      if _numDGCbonds >= _numFJCbonds:
        print('...... approximating chain as overall DGC')
        MolPatches[chain_name] = {'Stat':'DGC'}
      else:
        print('...... approximating chain as overall FJC')
        MolPatches[chain_name] = {'Stat':'FJC'}


  #patch some missing monomer and chain definitions when using MD topologies
  print('\n--> Final patching of interactions and molecule definitions')
  kuhns = [ float(bk)/float(bref) for bk in kuhns ]
  if 'KuhnLen' not in tmp['Monomers']: 
    tmp['Monomers']['KuhnLen'] = kuhns
  if 'GaussSmearWidth' not in tmp['Monomers'] and len(set([type(_sma) for _sma in asmear])) != 1:    
    raise ValueError('Detected incomplete asmear definition {}'.format(asmear))
  if asmear[0] != None:
    tmp['Monomers']['GaussSmearWidth'] = asmear
  for ic,(mol_name,mol_def) in enumerate(chains.items()):
    chainid = 'Chain{}'.format(ic+1)
    if mol_name in MolPatches:
      tmp['Chains'][chainid]['Statistics'] = MolPatches[mol_name]['Stat']
    else:
      raise ValueError('chain {} has undefined chain statistics'.format(mol_name))



  # = composition = 
  print('\n===== DEFINING COMPOSITIONS =====')
  tmp['Model1']['Composition'] = OrderedDict()
  update(tmp['Model1']['Composition'],settings['Composition'],'Ensemble',default='canonical')
  #first calculate what expected CChainDensity, volfracs are if using box and system topology
  precision = settings['Composition'].get('Precision',10)
  beadfracs = chain_beadfrac + smallmol_beadfrac #first group together
  beadfracs = [np.round(f,decimals=precision) for f in beadfracs]
  beadfracs[-1] = np.round(1.0-sum(beadfracs[:-1]),decimals=precision) #adjust to ensure sums to one, arbitrarily choose last species
  chain_beadfrac = [ float(f) for ii,f in enumerate(beadfracs) if ii < len(chain_beadfrac) ]
  smallmol_beadfrac = [ float(f) for ii,f in enumerate(beadfracs) if ii >= len(chain_beadfrac) ]
  #TODO: need to adjust for charge neutrality
  #Get box volume. Precedence is CellLengths, then settings.Composition.Box
  CC = settings['Composition'].get('CChainDensity',None)
  inputchainfrac = settings['Composition'].get('ChainVolFrac',None)
  autochainfrac = (isinstance(inputchainfrac,str) and inputchainfrac.lower() in ['auto,automatic']) or (inputchainfrac==None)
  inputsmallmolfrac = settings['Composition'].get('SmallMoleculeVolFrac',None)
  autosmallmolfrac = (isinstance(inputsmallmolfrac,str) and inputsmallmolfrac.lower() in ['auto,automatic']) or (inputsmallmolfrac==None)
  if CC is None or (isinstance(CC,str) and CC.lower() in ['auto','automatic']): #need to infer box volume
    print('automatically inferring CChainDensity')
    if box_dim == 3:
      V = float(box[0]) * float(box[1]) * float(box[2]) * tmp['Model1']['Cell']['CellScaling']**3.0
    elif 'Box' in settings['Composition']:
      _refbox = parse2list(settings['Composition']['Box'])
      print('CellLengths dimension !=3, relying on backup box definition {}'.format(_refbox))
      if len(_refbox) == 3:
        V = float(_refbox[0]) * float(_refbox[1]) * float(_refbox[2])
      elif len(_refbox) == 1:
        V = float(_refbox[0])
      else:
        raise ValueError('reference box for calculating CChain Density should have dimension 3, or 1 (if just Volume)')
    CC = float(np.round(total_beadnum/V, decimals=precision))
    print('... {}'.format(CC))
  if tmp['Model1']['Composition']['Ensemble'].lower() == 'canonical':
    tmp['Model1']['Composition']['CChainDensity'] = CC
    if autochainfrac and not isinstance(chain_beadfrac,float) and len(chain_beadfrac) > 0:
      tmp['Model1']['Composition']['ChainVolFrac'] = chain_beadfrac
    elif not autochainfrac and inputchainfrac is not None > 0:
      tmp['Model1']['Composition']['ChainVolFrac'] = inputchainfrac
    if autosmallmolfrac and not isinstance(smallmol_beadfrac,float) and len(smallmol_beadfrac) > 0:
      tmp['Model1']['Composition']['SmallMoleculeVolFrac'] = smallmol_beadfrac 
    elif not autosmallmolfrac and inputsmallmolfrac is not None > 0:
      tmp['Model1']['Composition']['SmallMoleculeVolFrac'] = inputsmallmolfrac
  if tmp['Model1']['Composition']['Ensemble'].lower() == 'grand':
    inputchainactivity = parse2list(settings['Composition']['ChainActivity'])
    inputsmallmolactivity = parse2list(settings['Composition']['SmallMoleculeActivity'])
    if len(inputchainactivity) != len(chain_beadfrac):
      raise ValueError('chain activity list length incommensurate with number of chains ({})'.format(len(chain_beadfrac)))
    elif len(inputchainactivity) > 0:
      tmp['Model1']['Composition']['ChainActivity'] = inputchainactivity
    if len(inputsmallmolactivity) != len(smallmol_beadfrac):
      raise ValueError('chain activity list length incommensurate with number of chains ({})'.format(len(chain_beadfrac)))
    elif len(inputsmallmolactivity) > 0:
      tmp['Model1']['Composition']['SmallMoleculeActivity'] = inputsmallmolactivity

  # = operators #fairly straightforward =
  print('\n===== DEFINING OPERATORS =====')
  tmp['Model1']['Operators'] = OrderedDict()
  operator_settings = settings.get('Operators',{})
  update(tmp['Model1']['Operators'],operator_settings,'CalcHamiltonian',default=True)
  update(tmp['Model1']['Operators'],operator_settings,'CalcPressure',default=False)
  update(tmp['Model1']['Operators'],operator_settings,'CalcStressTensor',default=False)
  update(tmp['Model1']['Operators'],operator_settings,'CalcChemicalPotential',default=False)
  update(tmp['Model1']['Operators'],operator_settings,'CalcDensityOperator',default=False)
  update(tmp['Model1']['Operators'],operator_settings,'IncludeIdealGasTerms',default=True)
  update(tmp['Model1']['Operators'],operator_settings,'CalcStructureFactor',fill=False)
  update(tmp['Model1']['Operators'],operator_settings,'CalcOrientationCorrelator',fill=False)
  update(tmp['Model1']['Operators'],operator_settings,'OrientationCorr_SpatialAverageRange',
         fill = 'CalcOrientationCorrelator' in tmp['Model1']['Operators'],default=0.25)
   
  # = init fields = 
  print('\n===== DEFINING INIT FIELDS =====')
  tmp['Model1']['InitFields'] = OrderedDict()
  update(tmp['Model1']['InitFields'],settings['InitFields'],'ReadInputFields',default=False)
  update(tmp['Model1']['InitFields'],settings['InitFields'],'InputFieldsFile',
          fill = 'ReadInputFields' in tmp['Model1']['InitFields'], default='fields0_k.bin')
  for ib in range(n_bead_types):
    #TODO, e.g. special options. Right now, no validation, just assume that things are entered as an appropriate dict
    update(tmp['Model1']['InitFields'],settings['InitFields'],'InitField{}'.format(ib+1), default={'InitType':'URNG'})

  ### ===== simulation #mandatory: JobType, FieldUpdater, DT, LambdaForceScale, NumTimeStepsPerBlock, NumBlocks =====
  print('\n===== DEFINING SIMULATION PARAMETERS =====')
  spec['Simulation'] = OrderedDict()
  lambda_force = parse2list(settings['Simulation']['LambdaForceScale'])
  settings['Simulation']['LambdaForceScale'] = lambda_force
  if len(lambda_force) < n_bead_types:
    print('CAUTION: lambda_force {} underspecified for {} bead types'.format(lambda_force,n_bead_types))
  if len(lambda_force) > n_bead_types:
    print('CAUTION: lambda_force {} overspecified for {} bead types'.format(lambda_force,n_bead_types))
  update(spec['Simulation'],settings['Simulation'],'JobType',mandatory=True)
  update(spec['Simulation'],settings['Simulation'],'FieldUpdater',mandatory=True)
  update(spec['Simulation'],settings['Simulation'],'CellUpdater',default='Broyden')
  update(spec['Simulation'],settings['Simulation'],'TimeStepDT',mandatory=True)
  update(spec['Simulation'],settings['Simulation'],'LambdaForceScale',mandatory=True)
  update(spec['Simulation'],settings['Simulation'],'SCFTForceStoppingTol',default=1.0e-4)
  update(spec['Simulation'],settings['Simulation'],'VariableCell',default=False)
  update(spec['Simulation'],settings['Simulation'],'LambdaStressScale',
                            fill='VariableCell' in spec['Simulation'], default=1.0)
  update(spec['Simulation'],settings['Simulation'],'SCFTStressStoppingTol',
                            fill='VariableCell' in spec['Simulation'], default=1.0e-4)
  update(spec['Simulation'],settings['Simulation'],'NumTimeStepsPerBlock',mandatory=True)
  update(spec['Simulation'],settings['Simulation'],'NumBlocks',mandatory=True)
  #spec['Simulation']['JobType'] = settings['Simulation']['JobType']
  #spec['Simulation']['FieldUpdater'] = settings['Simulation']['FieldUpdater']
  #spec['Simulation']['TimeStepDT'] = settings['Simulation']['TimeStepDT']
  #spec['Simulation']['LambdaForceScale'] = lambda_force
  #spec['Simulation']['SCFTForceStoppingTol'] = settings['Simulation'].get('SCFTForceStoppingTol',1.0e-4)
  #spec['Simulation']['VariableCell'] = settings['Simulation'].get('VariableCell',False)
  #spec['Simulation']['LambdaStressScale'] = settings['Simulation'].get('LambdaStressScale',1.0)
  #spec['Simulation']['SCFTStressStoppingTol'] = settings['Simulation'].get('SCFTStressStoppingTol',1.0e-4)
  #spec['Simulation']['NumTimeStepsPerBlock'] = settings['Simulation'].get('NumTimeStepsPerBlock')
  #spec['Simulation']['NumBlocks'] = settings['Simulation']['NumBlocks']
  #spec['Simulation']['RandomSeed'] = settings['Simulation'].get('RandomSeed',0)

  ## === IO === 
  print('\n===== DEFINING I/O =====')
  spec['Simulation']['IO'] = OrderedDict()
  io_settings = settings.get('IO',{})
  update(spec['Simulation']['IO'], io_settings, 'KeepDensityHistory', fill=False)
  update(spec['Simulation']['IO'], io_settings, 'KeepFieldHistory', fill=False)
  update(spec['Simulation']['IO'], io_settings, 'DensityOutputByChain', fill=False)
  update(spec['Simulation']['IO'], io_settings, 'OutputFormattedFields', fill=False)
  update(spec['Simulation']['IO'], io_settings, 'OutputFields', default='HFields')
  update(spec['Simulation']['IO'], io_settings, 'FieldOutputSpace', default='both')
  #spec['Simulation']['IO']['KeepDensityHistory'] = io_settings.get('KeepDensityHistory',False)
  #spec['Simulation']['IO']['KeepFieldHistory'] = io_settings.get('KeepFieldHistory',False)
  #spec['Simulation']['IO']['DensityOutputByChain'] = io_settings.get('DensityOutputByChain',False)
  #spec['Simulation']['IO']['OutputFormattedFields'] = io_settings.get('OutputFormattedFields',False)
  #spec['Simulation']['IO']['OutputFields'] = io_settings.get('OutputFields','HFields')
  #spec['Simulation']['IO']['FieldOutputSpace'] = io_settings.get('FieldOutputSpace','both')

  ### ===== pll specifications =====
  print('\n===== DEFINING PARALLEL PARAMETERS =====')
  spec['Parallel'] = OrderedDict()
  pll_settings = settings.get('Parallel',{})
  update(spec['Parallel'],pll_settings,'CUDA_SelectDevice',default=0)
  update(spec['Parallel'],pll_settings,'CUDA_ThreadBlockSize',default=64)
  update(spec['Parallel'],pll_settings,'OpenMP_NThreads',default=8)
  #spec['Parallel']['CUDA_selectdevice'] = pll_settings.get('CUDA_SelectDevice',0)
  #spec['Parallel']['CUDA_threadblocksize'] = pll_settings.get('CUDA_ThreadBlockSize',64)
  #spec['Parallel']['OpenMP_nthreads'] = pll_settings.get('OpenMP_NThreads',8)

  ### final
  return paramdict(spec)


def dict_to_str(spec):
  '''
  Expect to be dictionary of dictionaries (possibly of dictionaries)
  Assume each key is either a dict, value, or list (of floats/numbers)
  Convert to list of lines first

  how to track level of indentation? is the wrapper tracking indentation, or is the recursion call taking care of indentation?
  Args:
    spec : dict-like
  '''
  lines = []
  # parse
  for key,val in spec.items():
    if isinstance(val,(dict,OrderedDict)):
      lines.append('\n{} {{'.format(key))
      indented = indent(dict_to_str(val),1)
      if isinstance(indented,str): lines.append(indented)
      else: lines.extend(indented)
      lines.append('}')
    elif isinstance(val,(float,int,str)):
      lines.append(key.ljust(21) + ' = {}'.format(val))
    elif isinstance(val,(list,tuple,np.ndarray)):
      lines.append(key.ljust(21) + ' = {}'.format(strsimple(val)))
    else:
      raise TypeError('unsupported type for value: {}'.format(val))

  # convert lines to str
  spec_string = '\n'.join(lines)

  # final
  #return lines
  return spec_string


def write(filename,spec,header=None):
  import os
  from datetime import datetime
  prefix,extension = os.path.splitext(filename)
  heading = 'Initially generated {}'.format(datetime.now())
  if header is None:
    header = heading
  else:
    header = '{}; {}'.format(heading,header)
  if not isinstance(spec,str):
    #should be a dict
    if isinstance(spec,paramdict):
      spec = OrderedDict(spec)
    yaml.save_dict(prefix+'.yaml', spec, header=header)
    spec = dict_to_str(spec)
  with open(filename,'w') as f:
    f.write('#{}\n'.format(header))
    f.write(spec)
    f.write('\n') #for some reason, pfts needs extra linebreak in order to parse correctly...

#Test
def test():
  pt = preprocess_topology('tests_pfts/system_cg.yaml')
  spec = generate_input('tests_pfts/system_cg.yaml','tests_pfts/example_ff.yaml','tests_pfts/input.yaml')
  specstring = dict_to_str(spec)
  write('tests_pfts/input_generated.in',spec)
  #print(specstring)
def test_minimal():
  pt = preprocess_topology('tests_pfts/system_cg.yaml')
  spec = generate_input('tests_pfts/system_cg.yaml','tests_pfts/example_ff.yaml','tests_pfts/input_minimal.yaml')
  specstring = dict_to_str(spec)
  write('tests_pfts/inputminimal_generated.in',spec)
  #print(specstring)

#test()
#test_minimal()



# OTHER TODOs:
# Run system

### doing simple substitutions
# fill in later; temporary work-around is to use the generated yaml file, which can be easily written out in pfts format! Maybe more verbose, but also more intuitive and direct, without messing with regexp, can easily add new sections, etc.
# i.e. functions to manage an input script, and all the different sections


# and some setters that don't do any validation, relies on getting right key name:
# no key removals, need to set everything correctly!!!
def set_composition(spec,**kwargs):
  subdict = spec['Models']['Model1']['Composition'] 
  if 'ChainVolFrac' in kwargs:
    chvol = kwargs['ChainVolFrac']
  else:
    chvol = subdict.get('ChainVolFrac',[])
  if 'SmallMoleculeVolFrac' in kwargs:
    smvol = kwargs['SmallMoleculeVolFrac']
  else:
    smvol = subdict.get('SmallMoleculeVolFrac',[])
  precision = 10

  chaindefs = [ val for key,val in spec['Models']['Chains'].items() if key.startswith('Chain') ]
  if len(chvol) != len(chaindefs): 
    raise ValueError('input chain vol frac {} has diff. # species than defined')
  smdefs = [ val for key,val in spec['Models']['SmallMolecules'].items() if key.startswith('SmallMolecule') ]
  if len(chvol) != len(chaindefs): 
    raise ValueError('input smallmol vol frac {} has diff. # species than defined')
  charges = spec['Models']['Monomers'].get('Charge',None)
 
  bead_nums = []
  for chdef in chaindefs:
    bead_num = chdef['NPerBlock']
    bead_nums.append( sum(bead_num) )
  for smdef in smdefs:
    bead_nums.append(1)
  nmoltypes = len(chaindefs) + len(smdefs)
  bead_fracs = chvol + smvol
  for ii in range(nmoltypes): #now every beadfrac is a nice multiple of 1/10**precision/Nbead
    bead_fracs[ii] = 10.**-precision * bead_nums[ii] * round( bead_fracs[ii]/10.**-precision/bead_nums[ii] )

  mol_charges = None
  if charges is not None:
    print('bead fracs before adjusting charges: {}'.format(bead_fracs))
    #adjust for charge neutrality. not fully general, can't handle all cases, where numbers are not easily expressed in decimal (e.g. chain charge frac of 1/3, 1/7)
    #get chain species charges
    ch_charges = []
    sm_charges = []
    for chdef in chaindefs:
      bead_charges = [charges[ind-1] for ind in chdef['BlockSpecies']]
      bead_num = chdef['NPerBlock']
      charge_tuple = zip(bead_charges,bead_num)
      netcharge = sum( [_chg*_n for _chg,_n in charge_tuple] )
      ch_charges.append(netcharge)
    for smdef in smdefs:
      sm_charges.append( charges[smdef['Species']-1] )

    mol_charges = ch_charges + sm_charges
    charge_fracs = zip( bead_fracs, bead_nums, mol_charges )
    charge_fracs = [ _vfrac*_chg/float(_bn) for _vfrac,_bn,_chg in charge_fracs ]

    #now adjust charges for neutrality: adjust amount of the last charged species w/ nonzero beadfrac
    adjustindex = None
    print(charge_fracs)
    for revindex in range(1,nmoltypes+1):
      if charge_fracs[-revindex] != 0.0:
        adjustindex = nmoltypes - revindex
        break
    if adjustindex != None: 
      othercharges = sum( charge_fracs[:adjustindex] )
      charge_fracs[adjustindex] = - othercharges
      bead_fracs[adjustindex] = round(charge_fracs[adjustindex]/mol_charges[adjustindex],precision)
      print('adjusting charge fraction with index {}'.format(adjustindex))
    print('bead fracs after adjusting charges: {}'.format(bead_fracs))

  #Adjust for incompressibility as needed, rely on neutral species
  adjustindex = None
  for revindex in range(1,nmoltypes+1):
    if bead_fracs[-revindex] != 0.0:
      if mol_charges is None or mol_charges[-revindex]==0.0:
        adjustindex = nmoltypes - revindex
        break
  if adjustindex is not None:
    print('adjusting total bead fraction with index {}'.format(adjustindex))
    bead_fracs[adjustindex] = round(1.0 - sum([_bf for ii,_bf in enumerate(bead_fracs) if ii!=adjustindex]),precision)
  print('bead fracs after adjusting total bead fraction: {}'.format(bead_fracs))

  chvol = [ _bf for ii,_bf in enumerate(bead_fracs) if ii < len(chaindefs) ]
  smvol = [ _bf for ii,_bf in enumerate(bead_fracs) if ii >= len(chaindefs) ]
  if len(chvol) > 0:
    kwargs['ChainVolFrac'] = chvol
  if len(smvol) > 0:
    kwargs['SmallMoleculeVolFrac'] = smvol

  for key in kwargs:
    update(subdict,kwargs,key)

