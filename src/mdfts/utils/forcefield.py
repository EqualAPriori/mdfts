#! General helper functions for processing and manging force field definitions/files
#  In future, can include general helpers and function expressions, e.g. so that it's easy to construct an openmm custom force.
# TODO
# 1) save forcefield out in shorthand as well, b/c the processed file is quite verbose and not quite readable, plus still has some funky whitespace issues. Essentially, mainly want to dump specific potentials in one-line format
# 2) checks/tests/validation?
# 3) convert defaults sections, not just the params sections
#

import numpy as np
import yamlhelper as yaml
import os,copy
import parsify
from collections import OrderedDict

verbose = True
def vprint(string):
  if verbose: print(string)

class ForceField():
  '''
  Parameters
  ----------
  string or dict

  Notes
  -----
  Need to be careful not to store numpy floats and ints, because ruamel doesn't know how to process them.

  Examples
  --------
  >>> ff = forcefield.ForceField('example_ff.yaml')
  no section found for standard parameters
  >>> ff.sim2md()
  >>> ff.save('example_ff')
  >>> ff.processed_file['bond_harmonic']['params_md'][0]['k']['val']=10.0
  >>> ff.md2sim()
  >>> ff.save('example_ff')
  '''
  def __init__(self,ffdef=None):
    self.loaded_file = None
    self.processed_file = None

    if ffdef is not None:
      self.load(ffdef)
 
  def load(self,ffdef):
    if isinstance(ffdef,str):
      self.loaded_file = yaml.load(ffdef)
      self.processed_file = yaml.load(ffdef)
      self.process_ff_dict()
 
  def save(self,prefix):
    filename = '{}_processed_ff.yaml'.format(prefix)
    yaml.save_dict(filename, self.processed_file, header='#Processed force field definition.\n')


  def process_ff_dict(self):
    #create aliases for easy access
    loadedff = self.loaded_file
    outputff = self.processed_file

    #process sections one by one
    for ff_type_name in ff_types:
      if ff_type_name in loadedff:
        ff0 = loadedff[ff_type_name]
        ff1 = outputff[ff_type_name]
        ff_class_dict[ff_type_name].expand(ff0,ff1)

  def sim2md(self):
    #create aliases for easy access
    loadedff = self.loaded_file
    outputff = self.processed_file

    #process sections one by one
    for ff_type_name in ff_types:
      if ff_type_name in loadedff:
        ff1 = outputff[ff_type_name]
        ff_class_dict[ff_type_name].sim2md(ff1)

  def md2sim(self):
    #create aliases for easy access
    loadedff = self.loaded_file
    outputff = self.processed_file

    #process sections one by one
    for ff_type_name in ff_types:
      if ff_type_name in loadedff:
        ff = loadedff[ff_type_name]
        ff1 = outputff[ff_type_name]
        ff_class_dict[ff_type_name].md2sim(ff1)


# ===== potential-specific sections, for organization =====
def parse_default(ffdict,outdict,section,field,fixable=True,defaultdefault=0.0):
  '''
  Notes
  -----
  Need to distinguish b/t parameters that can be optimized/toggled, and those that can't.
  '''
  if section in ffdict:
    if field in ffdict[section]:
      parsed = parsify.parse_entry( ffdict[section][field] )
      if 'val' not in parsed:
        parsed['val'] = defaultdefault
      if 'fixed' not in parsed and fixable:
        parsed['fixed'] = True
    else:
      vprint('Tried to access parameter/field {} in {}, but does not exist. Filling in with ({},{})'.format(field,section,defaultdefault,True))
      parsed = {'val':defaultdefault}
      if fixable:
        parsed['fixed'] = True

  else:
    vprint('Tried to access parameter/field {} in {}, but does not exist.'.format(field,section))
    parsed = {'val':defaultdefault}
    if fixable:
      parsed['fixed'] = True
    if section not in outdict:
      outdict[section] = yaml.YAML.comments.CommentedMap()
      ffdict[section] = yaml.YAML.comments.CommentedMap()

  outdict[section][field] = parsed

  return parsed
  
def fill_defaults(ffdict,defaults,paramfields,options=None):
  for param_name in paramfields:
    if param_name not in ffdict:
      ffdict[param_name] = copy.deepcopy(defaults[param_name])
    else:
      for k,v in defaults[param_name].items(): #i.e. check value, fixed
        if k not in ffdict[param_name]:
          ffdict[param_name][k] = copy.deepcopy(v)

  for option_name in options:
    if option_name not in ffdict:
      ffdict[option_name] = copy.deepcopy(defaults[option_name])

  return ffdict

class base_potential():
  '''To use as a reference. In addition, the static methods are able to do most of the common heavy lifting.
  '''
  fields_sim = {'param1':{'fixable':True, 'defaultdefault':0.0}}
  nbody = 1
  prefix = 'base'

  @staticmethod
  def expand(cls,ffdict,outdict):
    '''Expands shorthand of both sim and md sections of parameters, filling in defaults
    Notes
    -----
    Required variables in defaults: FConst, Dist
    '''
    vprint('--->Reading Defaults')
    #read in defaults first
    section = 'defaults_sim'
    sim_default = {}
    for field,v in cls.fields_sim.items():
      sim_default[field] = parse_default( ffdict,outdict, section,field, fixable = v['fixable'], defaultdefault = v['defaultdefault'] )
    
    section = 'defaults_md'
    md_default = {}
    for field,v in cls.fields_md.items():
      md_default[field] = parse_default( ffdict,outdict, section,field, fixable = v['fixable'], defaultdefault = v['defaultdefault'] )

    #read in "special options" that don't fit the standard parameter format above
    #for special options, don't do special formatting/parsing. require user to get it right
    #right now this is for setting up required options that may not have been defined
    if len(cls.options_sim) > 0: 
      sections = ['defaults_sim']
      for section in sections:
        for field,v in cls.options_sim.items():
          if section not in outdict:
            outdict[section] = yaml.YAML.comments.CommentedMap()
            ffdict[section] = yaml.YAML.comments.CommentedMap()
          if field not in ffdict[section]:
            outdict[section][field] = copy.deepcopy(v)
          if section == 'defaults_sim':
            sim_default[field] = copy.deepcopy(outdict[section][field])
    if len(cls.options_md) > 0: 
      sections = ['defaults_md']
      for section in sections:
        for field,v in cls.options_md.items():
          if section not in outdict:
            outdict[section] = yaml.YAML.comments.CommentedMap()
            ffdict[section] = yaml.YAML.comments.CommentedMap()
          if field not in ffdict[section]:
            outdict[section][field] = copy.deepcopy(v)
          if section == 'defaults_md':
            md_default[field] = copy.deepcopy(outdict[section][field])


              

    #now expand and fill in defaults
    vprint('--->Expanding and filling in defaults')
    section = 'params_sim'
    if section in ffdict:
      if ffdict[section] is not None:
        for ientry,entry in enumerate(ffdict[section]):
          if type(outdict[section][ientry]) is str or isinstance(outdict[section][ientry],list):
            outdict[section][ientry] = yaml.YAML.comments.CommentedMap()
          parsed = parsify.parse_potential_entry(entry,cls.nbody,store_dict = outdict[section][ientry], prefix = cls.prefix)
          fill_defaults( outdict[section][ientry], sim_default, cls.fields_sim, cls.options_sim )
      else:
        vprint('no section found for sim parameters')

    section = 'params_md'
    if section in ffdict:
      if ffdict[section] is not None:
        for ientry,entry in enumerate(ffdict[section]):
          if type(outdict[section][ientry]) is str or type(outdict[section][ientry]) is list:
            outdict[section][ientry] = yaml.YAML.comments.CommentedMap()
          parsed = parsify.parse_potential_entry(entry,cls.nbody,store_dict = outdict[section][ientry], prefix = cls.prefix)
          fill_defaults( outdict[section][ientry], md_default, cls.fields_md, cls.options_md )
      else:
        vprint('no section found for standard parameters')


  @staticmethod
  def convert(cls,ffdict,source,target):
    '''Convert parameters between conventions.
    Notes
    -----
    In this base class, just do the essentials
    '''
    source = 'params_{}'.format(source)
    target = 'params_{}'.format(target)
    if source in ffdict and ffdict[source] is not None:
      #if target not in ffdict or ffdict[target] is None:
      #  ffdict[target] = yaml.YAML.comments.CommentedSeq()
      #  for ir in range(len(ffdict[source])):
      #    ffdict[target].append( yaml.YAML.comments.CommentedMap() )
      ffdict[target] = yaml.YAML.comments.CommentedSeq()
      for ir in range(len(ffdict[source])):
        ffdict[target].append( yaml.YAML.comments.CommentedMap() )
      for ip,potential in enumerate(ffdict[source]):
        ffdict[target][ip]['species'] = copy.deepcopy(potential['species'])
        ffdict[target][ip]['name'] = copy.deepcopy(potential['name'])
    else:
      raise ValueError('{} parameters not defined for potential of type {}'.format(source,cls))
 
  @staticmethod
  def sim2md(cls,ffdict):
    '''Convert sim-convention parameters to standard-convention parameters
    Notes
    -----
    In this base class, just do the essentials
    '''
    source = 'params_sim'
    target = 'params_md'
    if source in ffdict:
      if target not in ffdict or ffdict[target] is None:
        ffdict[target] = yaml.YAML.comments.CommentedSeq()
      for ir in range(len(ffdict[source])):
        ffdict[target].append( yaml.YAML.comments.CommentedMap() )
      for ip,potential in enumerate(ffdict[source]):
        ffdict[target][ip]['species'] = copy.deepcopy(potential['species'])
        ffdict[target][ip]['name'] = copy.deepcopy(potential['name'])
    else:
      raise ValueError('sim parameters not defined for potential of type {}'.format(cls))
 
  @staticmethod
  def md2sim(cls,ffdict):
    '''Convert sim-convention parameters to standard-convention parameters
    Notes
    -----
    In this base class, just do the essentials (i.e. copy species, name).
    '''
    source = 'params_md'
    target = 'params_sim'
    if source in ffdict:
      if target not in ffdict or ffdict[target] is None:
        ffdict[target] = yaml.YAML.comments.CommentedSeq()
      for ir in range(len(ffdict[source])):
        ffdict[target].append( yaml.YAML.comments.CommentedMap() )
      for ip,potential in enumerate(ffdict[source]):
        ffdict[target][ip]['species'] = copy.deepcopy(potential['species'])
        ffdict[target][ip]['name'] = copy.deepcopy(potential['name'])
    else:
     raise ValueError('md parameters not defined for potential of type {}'.format(cls))
 


class bond_harmonic():
  '''
  Notes
  -----
  Sim: see sim/potential/bond.py
    U = FConst * ( r - Dist0 )^2

  OpenMM/HooMD: http://docs.openmm.org/6.1.0/userguide/theory.html
    U = 0.5 * k * ( r - r0 )^2
  '''
  options_sim = {}
  options_md = {}

  fields_sim = OrderedDict({'FConst':{'fixable':True, 'defaultdefault':1.0},
            'Dist0':{'fixable':True, 'defaultdefault':0.0}})

  fields_md = OrderedDict({'FConst':{'fixable':True, 'defaultdefault':1.0},
            'Dist0':{'fixable':True, 'defaultdefault':0.0}})
  nbody = 2
  prefix = 'bond'

  @classmethod
  def expand(cls,ffdict,outdict):
    '''Expands shorthand of both sim and md sections of parameters, filling in defaults
    Parameters
    ----------
    cls : class
    ffdict : dict
      for bond_harmonic
    outdict: dict
      for bond_harmonic

    Notes
    -----
    Required variables in defaults: FConst, Dist
    '''
    vprint('\n=== Parsing Harmonic Bond ===')
    base_potential.expand(cls,ffdict,outdict)

  @classmethod
  def sim2md(cls,ffdict):
    '''Converts sim interaction parameters to "standard" interaction parameters. Assume forcefields are fully defined already (e.g. filled out via expand()).
    Notes
    -----
    r0 = Dist0
    k = 2*Fconst
    '''
    vprint('\n=== Converting Harmonic Bond from Sim to MD ===')
    base_potential.convert(cls,ffdict,'sim','md')
    source = 'params_sim'
    target = 'params_md'
    if source in ffdict:
      for ip,potential in enumerate(ffdict[source]):
        ffdict[target][ip]['k'] = { 'val':2.0*potential['FConst']['val'], 'fixed':potential['FConst']['fixed'] }
        ffdict[target][ip]['r0'] = copy.deepcopy(potential['Dist0'])
    else:
      raise ValueError('sim parameters not defined for potential of type {}'.format(cls))

  @classmethod
  def md2sim(cls,ffdict):
    '''Converts sim interaction parameters to "standard" interaction parameters. Assume forcefields are fully defined already.
    Notes
    -----
    r0 = Dist0
    k = 2*Fconst
    '''
    vprint('\n=== Converting Harmonic Bond from MD to Sim ===')
    base_potential.convert(cls,ffdict,'md','sim')
    source = 'params_md'
    target = 'params_sim'
    if source in ffdict:
      for ip,potential in enumerate(ffdict[source]):
        ffdict[target][ip]['FConst'] = { 'val':0.5*potential['k']['val'], 'fixed':potential['k']['fixed'] }
        ffdict[target][ip]['Dist0'] = copy.deepcopy(potential['r0'])
    else:
      raise ValueError('md parameters not defined for potential of type {}'.format(cls))


class pair_ljg():
  '''
  Notes
  -----
  Sim: see sim/potential/ljgaussian.py
    U = 4 * Epsilon * ( (Sigma/r)^12 - (Sigma/r)^6 ) + B * exp( - Kappa* ( r - Dist0 )^2 )

  OpenMM/HooMD: http://docs.openmm.org/6.1.0/userguide/theory.html
    U = 4 * Epsilon * ( (Sigma/r)^12 - (Sigma/r)^6 ) + (u0/8pi^1.5a^3)* exp( - 0.5 ( r - Dist0 )^2 / 4a^2 )
    U = 4 * Epsilon * ( (Sigma/r)^12 - (Sigma/r)^6 ) + (u0/(2pi)^1.5s^3)* exp( - 0.5 ( r - Dist0 )^2 / 2s^2 )
    where 2a^2= ai^2 + aj^2
    and the "traditional" sigma^2 = 2a^2

    B = u0/(2pi)^1.5 s^3
  '''
  options_sim = {'Shift': False}
  options_md = {'shift': False}
  
  fields_sim = OrderedDict([('B',{'fixable':True, 'defaultdefault':0.0}),
                ('Kappa',{'fixable':True, 'defaultdefault':1.0}),
                ('Dist0',{'fixable':True, 'defaultdefault':0.0}),
                ('Epsilon',{'fixable':True, 'defaultdefault':0.0}),
                ('Sigma',{'fixable':True, 'defaultdefault':1.0}),
                ('Cut',{'fixable':False, 'defaultdefault':1.0})
               ])
  fields_md = OrderedDict([('u0',{'fixable':True, 'defaultdefault':0.0}),
                ('sigma_g',{'fixable':True, 'defaultdefault':1.0}),
                ('r0',{'fixable':True, 'defaultdefault':0.0}),
                ('epsilon',{'fixable':True, 'defaultdefault':0.0}),
                ('sigma',{'fixable':True, 'defaultdefault':0.0}),
                ('cut',{'fixable':True, 'defaultdefault':1.0})
              ])

  nbody = 2
  prefix = 'ljg'


  @classmethod
  def expand(cls,ffdict,outdict):
    '''Expands shorthand of both sim and md sections of parameters, filling in defaults'''
    vprint('\n=== Parsing LJ-Gaussian ===')
    base_potential.expand(cls,ffdict,outdict)
    # special processing of the smearing parameters? or assume a particular given form.
    # At this point Kappa and sigma should be defined already. If they weren't provided in the loaded ff (i.e. they were default-filled), then need to override with some smearing  
    # Kappa = 1/2(ai^2 + aj^2)
    vprint('--->Accounting for mixing rules')
    style = 'sim' 
    section = 'params_'+style
    if 'a_smear' in ffdict['defaults_'+style]:
      asmears = ffdict['defaults_'+style]['a_smear']
      for ientry,entry in enumerate(ffdict[section]): #override if Kappa wasn't defined in original ff definition
        loaded_entry = parsify.parse_potential_entry(entry,cls.nbody,prefix = cls.prefix)
        vprint('entry {}, {}'.format(ientry,entry))
        if 'Kappa' in loaded_entry  and 'val' in loaded_entry['Kappa']:
          #then Kappa was defined, don't override 
          vprint('  ljg entry defined, no asmear overide')
        else:
          #no type-checking, i.e. if smearing bead type not defined, raise error
          asmears1 = [ asmears[beadtype] for beadtype in loaded_entry['species'][0] ]
          asmears2 = [ asmears[beadtype] for beadtype in loaded_entry['species'][1] ]

          if ( all( [ a == asmears1[0] for a in asmears1 ] ) and
               all( [ a == asmears2[0] for a in asmears2 ] )
             ):
            Kappa = 0.5/(asmears1[0]**2.0 + asmears2[0]**2.0)
            outdict[section][ientry]['Kappa']['val'] = Kappa
            vprint('  Using asmear mixing rule, Kappa = {}'.format(Kappa))
            outdict[section][ientry]['Cut']['val'] = 5.0*(0.5*(asmears1[0]**2.0+asmears2[0]**2.0))**0.5
            vprint('  with cutoff = 5abar = {}'.format( outdict[section][ientry]['Cut']['val'] ))
          else:
            vprint('  CAUTION: Multiple bead pairs defined, but asmears are not consistent, cannot define unique Kappa for this potential. Using default.'.format())

    style = 'md' 
    section = 'params_'+style
    if 'a_smear' in ffdict['defaults_'+style]:
      asmears = ffdict['defaults_'+style]['a_smear']
      for ientry,entry in enumerate(ffdict[section]): #override if Kappa wasn't defined in original ff definition
        loaded_entry = parsify.parse_potential_entry(entry,cls.nbody,prefix = cls.prefix)
        vprint('entry {}, {}'.format(ientry,entry))
        if 'sigma_g' in loaded_entry  and 'val' in loaded_entry['sigma_g']:
          #then Kappa was defined, don't override 
          vprint('  ljg entry defined, no asmear overide')
        else:
          #no type-checking, i.e. if smearing bead type not defined, raise error
          asmears1 = [ asmears[beadtype] for beadtype in loaded_entry['species'][0] ]
          asmears2 = [ asmears[beadtype] for beadtype in loaded_entry['species'][1] ]

          if ( all( [ a == asmears1[0] for a in asmears1 ] ) and
               all( [ a == asmears2[0] for a in asmears2 ] )
             ):
            sigma_g = (asmears1[0]**2.0 + asmears2[0]**2.0)**0.5
            outdict[section][ientry]['sigma_g']['val'] = sigma_g
            vprint('  Using asmear mixing rule, sigma_g = {}'.format(sigma_g))
            outdict[section][ientry]['cut']['val'] = 5.0*(0.5*(asmears1[0]**2.0+asmears2[0]**2.0))**0.5
            vprint('  with cutoff = 5abar = {}'.format( outdict[section][ientry]['Cut']['val'] ))
          else:
            vprint('  CAUTION: Multiple bead pairs defined, but asmears are not consistent, cannot define unique sigma_g for this potential. Using default.'.format())


  @classmethod
  def sim2md(cls,ffdict):
    '''Have to convert B,Kappa to u0,sigma_g. 
    Notes
    -----
    Also copy over default smearing lengths, check that Kappa is fixed.
    '''
    vprint('\n=== Converting LJGaussian from Sim to MD ===')
    base_potential.convert(cls,ffdict,'sim','md')
    source = 'params_sim'
    target = 'params_md'
    if source in ffdict:
      for ip,potential in enumerate(ffdict[source]):
        Kappa = potential['Kappa']['val']
        B = potential['B']['val']

        sigma_g = (0.5/Kappa)**0.5
        u0 = B * sigma_g**3.0 * (2.0*np.pi)**1.5

        ffdict[target][ip]['u0'] = { 'val':u0, 'fixed':potential['Kappa']['fixed'] }
        ffdict[target][ip]['sigma_g'] = { 'val':sigma_g, 'fixed':potential['B']['fixed'] }
        ffdict[target][ip]['r0'] = copy.deepcopy(potential['Dist0'])
        ffdict[target][ip]['sigma'] =  copy.deepcopy(potential['Sigma'])
        ffdict[target][ip]['epsilon'] = copy.deepcopy(potential['Epsilon'])
        ffdict[target][ip]['rcut'] = copy.deepcopy(potential['Cut'])
    else:
      raise ValueError('sim parameters not defined for potential of type {}'.format(cls))
    ffdict['defaults_md']['a_smear'] = copy.deepcopy(ffdict['defaults_sim']['a_smear'])


  @classmethod
  def md2sim(cls,ffdict):
    '''Have to convert u0,sigma_g to B,Kappa'''
    vprint('\n=== Converting LJGaussian from MD to Sim ===')
    base_potential.convert(cls,ffdict,'md','sim')
    source = 'params_md'
    target = 'params_sim'
    if source in ffdict:
      for ip,potential in enumerate(ffdict[source]):
        sigma_g = potential['sigma_g']['val']
        u0 = potential['u0']['val']
        Kappa = 0.5/sigma_g**2.0
        B = u0/(2.0*np.pi)**1.5/sigma_g**3.0
      
        ffdict[target][ip]['B'] = {'val':B}
        ffdict[target][ip]['Kappa'] = {'val':Kappa}
        ffdict[target][ip]['Dist0'] = copy.deepcopy(potential['r0'])
        ffdict[target][ip]['Sigma'] = copy.deepcopy(potential['sigma'])
        ffdict[target][ip]['Epsilon'] = copy.deepcopy(potential['epsilon'])
        ffdict[target][ip]['rcut'] = copy.deepcopy(potential['rcut'])
    else:
      raise ValueError('md parameters not defined for potential of type {}'.format(cls))
    ffdict['defaults_sim']['a_smear'] = copy.deepcopy(ffdict['defaults_md']['a_smear'])


class external_sin():
  '''
  Notes
  -----
  U = UConst*sin( (r - PlaneLoc)/(2 pi NPeriods) ), where `r` is along PlaneAxis
  '''
  options_sim = {}
  options_md = {}
  
  fields_sim = OrderedDict([('UConst',{'fixable':True, 'defaultdefault':0.0}),
                ('NPeriods',{'fixable':True, 'defaultdefault':1.0}),
                ('PlaneLoc',{'fixable':False, 'defaultdefault':0.0}),
                ('PlaneAxis',{'fixable':False, 'defaultdefault':0}),
               ])
  fields_md = OrderedDict([('UConst',{'fixable':True, 'defaultdefault':0.0}),
                ('NPeriods',{'fixable':True, 'defaultdefault':1.0}),
                ('PlaneLoc',{'fixable':False, 'defaultdefault':0.0}),
                ('PlaneAxis',{'fixable':False, 'defaultdefault':0}),
              ])

  nbody = 1
  prefix = 'extsin'

  @classmethod
  def expand(cls,ffdict,outdict):
    '''Expands shorthand of both sim and md sections of parameters, filling in defaults
    Parameters
    ----------
    cls : class
    ffdict : dict
    outdict: dict

    Notes
    -----
    Required variables in defaults: FConst, Dist
    '''
    vprint('\n=== Parsing External Sin ===')
    base_potential.expand(cls,ffdict,outdict)

  @classmethod
  def sim2md(cls,ffdict):
    '''Converts sim interaction parameters to "standard" interaction parameters. Assume forcefields are fully defined already (e.g. filled out via expand()).
    '''
    vprint('\n=== Converting External Sin from Sim to MD ===')
    base_potential.convert(cls,ffdict,'sim','md')
    source = 'params_sim'
    target = 'params_md'
    if source in ffdict:
      for ip,potential in enumerate(ffdict[source]):
        ffdict[target][ip]['UConst'] = copy.deepcopy(potential['UConst'])
        ffdict[target][ip]['NPeriods'] = copy.deepcopy(potential['NPeriods'])
        ffdict[target][ip]['PlaneLoc'] = copy.deepcopy(potential['PlaneLoc'])
        ffdict[target][ip]['PlaneAxis'] = copy.deepcopy(potential['PlaneAxis'])
    else:
      raise ValueError('sim parameters not defined for potential of type {}'.format(cls))

  @classmethod
  def md2sim(cls,ffdict):
    '''Converts sim interaction parameters to "standard" interaction parameters. Assume forcefields are fully defined already.
    '''
    vprint('\n=== Converting External Sin from MD to Sim ===')
    base_potential.convert(cls,ffdict,'md','sim')
    source = 'params_md'
    target = 'params_sim'
    if source in ffdict:
      for ip,potential in enumerate(ffdict[source]):
        ffdict[target][ip]['UConst'] = copy.deepcopy(potential['UConst'])
        ffdict[target][ip]['NPeriods'] = copy.deepcopy(potential['NPeriods'])
        ffdict[target][ip]['PlaneLoc'] = copy.deepcopy(potential['PlaneLoc'])
        ffdict[target][ip]['PlaneAxis'] = copy.deepcopy(potential['PlaneAxis'])
    else:
      raise ValueError('md parameters not defined for potential of type {}'.format(cls))

class coulomb_smeared():
  '''
  Notes
  -----
  automatically turns on ewald, then adds the correction potential
  if BornA = 0, then same as unsmeared Coulomb. Even if no interactions defined, still add ewald.
  
  first implementation: only overrides if BornA not defined below. Use same smearing convention as Gaussian above, not Born-A convention.

  in Fourier space: 4pi * Coef * exp(-k^2(ai^2+aj^2)/2)/k^2 --> 4pi * Coef * exp(-k^2a^2)/k^2
     Real space: Coef * erf(r/2a)/r

  i.e. a^2 = (ai^2 + aj^2)/2
       aborn = a*sqrt(pi)

  Asmears are given in the standard convention. Note that Sim implemented with the BornA convention.

  '''
  options_sim = {'ExcludeBondOrd':0, 'Shift':True}
  options_md = {'ExcludeBondOrd':0, 'Shift':True}
  
  fields_sim = OrderedDict([('Coef',{'fixable':True, 'defaultdefault':0.0}),
                ('BornA',{'fixable':True, 'defaultdefault':1.0}),
                ('Cut',{'fixable':False, 'defaultdefault':1.0}),
               ])
  fields_md = OrderedDict([('lb',{'fixable':True, 'defaultdefault':0.0}),
                ('a_smear',{'fixable':True, 'defaultdefault':1.0}),
                ('rcut',{'fixable':False, 'defaultdefault':1.0}),
              ])

  nbody = 2
  prefix = 'coulomb_smeared'

  @classmethod
  def expand(cls,ffdict,outdict):
    '''Expands shorthand of both sim and md sections of parameters, filling in defaults
    Parameters
    ----------
    cls : class
    ffdict : dict
    outdict: dict

    Notes
    -----
    '''
    vprint('\n=== Parsing Smeared Coulomb ===')
    base_potential.expand(cls,ffdict,outdict)

    vprint('--->Accounting for mixing rules')
    style = 'sim' 
    section = 'params_'+style
    if 'a_smear' in ffdict['defaults_'+style]:
      asmears = ffdict['defaults_'+style]['a_smear']
      for ientry,entry in enumerate(ffdict[section]): #override if BornA wasn't defined in original ff definition
        loaded_entry = parsify.parse_potential_entry(entry,cls.nbody,prefix = cls.prefix)
        vprint('entry {}, {}'.format(ientry,entry))
        if 'BornA' in loaded_entry  and 'val' in loaded_entry['BornA']:
          #then BornA was defined, don't override 
          vprint('  smeared coulomb BornA defined, no asmear overide')
        else:
          #no type-checking, i.e. if smearing bead type not defined, will get a key error
          asmears1 = [ asmears[beadtype] for beadtype in loaded_entry['species'][0] ]
          asmears2 = [ asmears[beadtype] for beadtype in loaded_entry['species'][1] ]

          if ( all( [ a == asmears1[0] for a in asmears1 ] ) and
               all( [ a == asmears2[0] for a in asmears2 ] )
             ):
            BornA = np.pi**0.5 * (0.5*(asmears1[0]**2.0 + asmears2[0]**2.0))**0.5
            outdict[section][ientry]['BornA']['val'] = BornA
            vprint('  Using asmear mixing rule, BornA = {}'.format(BornA))
          else:
            vprint('  CAUTION: Multiple bead pairs defined, but asmears are not consistent, cannot define unique smearing for this potential. Using default.'.format())

    style = 'md' 
    section = 'params_'+style
    if 'a_smear' in ffdict['defaults_'+style]:
      asmears = ffdict['defaults_'+style]['a_smear']
      for ientry,entry in enumerate(ffdict[section]): #override if SmearA wasn't defined in original ff definition
        loaded_entry = parsify.parse_potential_entry(entry,cls.nbody,prefix = cls.prefix)
        vprint('entry {}, {}'.format(ientry,entry))
        if 'SmearA' in loaded_entry  and 'val' in loaded_entry['SmearA']:
          #then SmearingA was defined, don't override 
          vprint('  smeared coulomb SmearA entry defined, no asmear overide')
        else:
          #no type-checking, i.e. if smearing bead type not defined, will get a key error
          asmears1 = [ asmears[beadtype] for beadtype in loaded_entry['species'][0] ]
          asmears2 = [ asmears[beadtype] for beadtype in loaded_entry['species'][1] ]

          if ( all( [ a == asmears1[0] for a in asmears1 ] ) and
               all( [ a == asmears2[0] for a in asmears2 ] )
             ):
            a_smear = (0.5*(asmears1[0]**2.0 + asmears2[0]**2.0))**0.5
            outdict[section][ientry]['a_smear']['val'] = a_smear
            vprint('  Using asmear mixing rule, a_smear = {}'.format(a_smear))
          else:
            vprint('  CAUTION: Multiple bead pairs defined, but asmears are not consistent, cannot define unique smearing for this potential. Using default.'.format())


  @classmethod
  def sim2md(cls,ffdict):
    '''Converts sim interaction parameters to "standard" interaction parameters. Assume forcefields are fully defined already (e.g. filled out via expand()).
    '''
    vprint('\n=== Converting Smeared Coulomb from Sim to MD ===')
    base_potential.convert(cls,ffdict,'sim','md')
    source = 'params_sim'
    target = 'params_md'
    if source in ffdict:
      for ip,potential in enumerate(ffdict[source]):
        ffdict[target][ip]['lb'] = copy.deepcopy(potential['Coef'])
        ffdict[target][ip]['a_smear'] = copy.deepcopy(potential['BornA'])
        ffdict[target][ip]['a_smear']['val'] = potential['BornA']['val']/np.pi**0.5
        ffdict[target][ip]['rcut'] = copy.deepcopy(potential['Cut'])
    else:
      raise ValueError('sim parameters not defined for potential of type {}'.format(cls))

  @classmethod
  def md2sim(cls,ffdict):
    '''Converts sim interaction parameters to "standard" interaction parameters. Assume forcefields are fully defined already.
    '''
    vprint('\n=== Converting Smeared Coulomb from MD to Sim ===')
    base_potential.convert(cls,ffdict,'md','sim')
    source = 'params_md'
    target = 'params_sim'
    if source in ffdict:
      for ip,potential in enumerate(ffdict[source]):
        ffdict[target][ip]['Coef'] = copy.deepcopy(potential['lb'])
        ffdict[target][ip]['BornA'] = copy.deepcopy(potential['a_smear'])
        ffdict[target][ip]['BornA']['val'] = potential['a_smear']['val']*np.pi**0.5
        ffdict[target][ip]['Cut'] = copy.deepcopy(potential['rcut'])
    else:
      raise ValueError('md parameters not defined for potential of type {}'.format(cls))


# ===== Helpful Bits =====
ff_types = ['bond_harmonic','pair_ljg','external_sin','coulomb_smeared']
ff_class_dict={'bond_harmonic':bond_harmonic, 'bond':bond_harmonic, 'pair_ljg':pair_ljg, 'ljg':pair_ljg,'external_sin':external_sin,'extsin':external_sin,'coulomb_smeared':coulomb_smeared,'coul_smeared':coulomb_smeared,'coulsmear':coulomb_smeared}

