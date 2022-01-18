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

#import topologify,forcefield
from collections import OrderedDict
import yamlhelper as yaml
#import networkx as nx
import numpy as np
import pandas as pd

## Helpers
def indent(text,n):
  '''
  sets indentation level.
  don't left-strip, to not ruin prior indentation
  note, even handles if have a list of multi-line entries, i.e. if an entry in a list has multiple lines, i.e. recursively go in to indent them!
  '''
  #print('trying to indent {}, {}'.format(type(text),text))
  if isinstance(text,str):
    t = text.split('\n')
    t = [('  ' * n) + line for line in t]
    t = '\n'.join(t)
  elif isinstance(text,(list,tuple)): #assume is list/tuple of line strings
    t = [('  ' * n) + line for line in text]

  return t

def strsimple(mylist,nspaces=1):
  '''
  turn [a,b,c] list into string: a b c
  only works if list is of numbers
  '''
  if isinstance(mylist,str):
    return mylist
  else:
    return (nspaces*' ').join(precisionprint(x) for x in mylist)
def precisionprint(x):
  if isinstance(x,float):
    return '{:.15f}'.format(x)
  else:
    return str(x)

def parse2list(s):
  '''assumes flat list. does not handle some of the nested syntax that parsify module defines'''
  if isinstance(s,(int,float,bool)):
    return [s]
  else:
    if isinstance(s,(list,tuple)):
      tmp = s
      #return s
    elif isinstance(s,str):
      tmp = s.split()
    else:
      raise ValueError('Unprocessable input {} of type {}'.format(s,type(s)))
    result = []
    for entry in tmp:
      if entry.lower() in ['true']:
        result.append(True)
      elif entry.lower() in ['false']:
        result.append(False)
      else:
        try: 
          x = int(entry)
          result.append(x)
          continue
        except:
          pass
        try:
          x = float(entry)
          result.append(x)
          continue
        except:
          pass
        result.append(entry)
    return result

def update(dict1,dict2,key,default=None,fill=True,mandatory=False):
  '''
  dict[key] = dict2[key]
  use default if key not in dict2
  don't insert key if fill = False

  mandatory=True raises error if key not in dict2. refers to whether mandatory in dict 2. fill refers to whether needs to be in dict 1
  '''
  if mandatory: #raise error if key not in dict2
    dict1[key] = dict2[key]
  elif fill: #not mandatory to have key in dict2, substitute default value
    dict1[key] = dict2.get(key,default)
  elif not fill and key in dict2 and not mandatory: #not mandatory to have key in dict2, not necessary to keep in dict1
    dict1[key] = dict2[key]

## loading a PFTS-type input
# 8) parse PFTS script into a dictionary?
#     sequence of # --> list
#     = --> :
#     \n --> , (unless is last element of that dict)
#     ignore comments
#    can probably do line by line, "increment/decrement" dictionary level
#
def parseline(line):
  tmp = line.strip().split('#')
  content=tmp[0]
  if len(tmp) > 1:
    comments='#'.join(tmp[1:]) #later can think about putting comments back
  return content

def getname(name):
  namemaps = {'inputfileversion':'InputFileVersion', 
      'nummodels':'NumModels','modeltype':'ModelType',
      'monomers':'Monomers','nspecies':'NSpecies','charge':'Charge',
      'kuhnlen':'KuhnLen','gausssmearwidth':'GaussSmearWidth',
      'chain':'Chain','nchains':'NChains','diffusermethod':'DiffuserMethod','contourds':'Contourds',
      'polymerreferencen':'PolymerReferenceN','label':'Label','architecture':'Architecture',
      'statistics':'Statistics','nblocks':'NBlocks','nbeads':'NBeads',
      'blockspecies':'BlockSpecies','nperblock':'NPerBlock',
      'smallmolecules':'SmallMolecules','nsmallmoleculetypes':'NSmallMoleculeTypes',
      'species':'Species',
      'dim':'Dim',
      'cell':'Cell','cellscaling':'CellScaling','celllengths':'CellLengths',
      'cellangles':'CellAngles','npw':'NPW','spacegroupname':'SpaceGroupName',
      'centertoprimitivecell':'CenterToPrimitiveCell','symmetrize':'Symmetrize',
      'interactions':'Interactions','applycompressibilityconstraint':'ApplyCompressibilityConstraint',
      'eelecstatic':'EElecStatic','inversezetan':'InverseZetaN',
      'composition':'Composition','ensemble':'Ensemble','cchaindensity':'CChainDensity',
      'chainvolfrac':'ChainVolFrac','smallmoleculevolfrac':'SmallMoleculeVolFrac',
      'chainactivity':'ChainActivity','smallmoleculeactivity':'SmallMoleculeActivity',
      'operators':'Operators','calchamiltonian':'CalcHamiltonian','calcpressure':'CalcPressure',
      'calcstresstensor':'CalcStressTensor','calcchemicalpotential':'CalcChemicalPotential',
      'calcdensityoperator':'CalcDensityOperator','includeidealgasterms':'IncludeIdealGasTerms',
      'calcstructurefactor':'CalcStructureFactor','calcorientationcorrelator':'CalcOrientationCorrelator',
      'orientationcorr_spatialaveragerange':'OrientationCorr_SpatialAverageRange',
      'initfields':'InitFields','readinputfields':'ReadInputFields','inputfieldsfile':'InputFieldsFile',
      'inittype':'InitType','initparameters':'InitParameters',
      'simulation':'Simulation','jobtype':'JobType','fieldupdater':'FieldUpdater','cellupdater':'CellUpdater','timestepdt':'TimeStepDT', 'dt':'TimeStepDT',
      'lambdaforcescale':'LambdaForceScale','scftforcestoppingtol':'SCFTForceStoppingTol', 
      'lambdaforce':'LambdaForceScale', 'scfttol':'SCFTForceStoppingTol','tolscft':'SCFTForceStoppingTol',
      'lambdastress':'LambdaStressScale', 'stresstol':'SCFTStressStoppingTol','tolstress':'SCFTStressStoppingTol',
      'variablecell':'VariableCell','lambdastressscale':'LambdaStressScale',
      'scftstressstoppingtol':'SCFTStressStoppingTol','numtimestepsperblock':'NumTimeStepsPerBlock',
      'numblocks':'NumBlocks','io':'IO','outputfields':'OutputFields','fieldoutputspace':'FieldOutputSpace',
      'parallel':'Parallel','cuda_selectdevice':'CUDA_SelectDevice',
      'cuda_threadblocksize':'CUDA_ThreadBlockSize','openmp_nthreads':'OpenMP_NThreads',
      'randomseed':'RandomSeed','keepdensityhistory':'KeepDensityHistory','keepfieldhistory':'KeepFieldHistory',
      'densityoutputbychain':'DensityOutputByChain','outputformattedfields':'OutputFormattedFields'}
  keys = list(namemaps.keys())
  extranames = {'model':'Model','chain':'Chain','smallmolecule':'SmallMolecule',
      'chi':'Chi','bexclvolume':'BExclVolume','initfield':'InitField'}
  foundkey = False
  name = name.strip()
  #if not foundkey:
  #  for k in keys:
  #    if k.lower() == name.lower():
  #      foundkey = True
  #      return namemaps[name.lower()]
  if name.lower() in namemaps:
    foundkey = True
    return namemaps[name.lower()]
  else:
    for k in extranames:
      if name.lower().startswith(k):
        suffix = name.lower().split(k)[1]
        foundkey = True
        return extranames[k] + suffix
  if not foundkey: #unrecognized option
    print('unrecognized field name: {}'.format(name))
    return name

def load(fname):
  '''Load PFTS-type input'''
  spec = OrderedDict()
  levels = [spec]
  with open(fname,'r') as f:
    l = f.readline()
    while l:
      l = parseline(l)
      if '{' in l: #start a new dictionary entry
        name = l.strip().split('{')[0]
        name = getname(name) #standardized format
        levels[-1][name] = OrderedDict()
        levels.append(levels[-1][name])
      if '}' in l: #end current dictionary entry
        levels.pop()
      elif '=' in l:
        name,entry = l.strip().split('=') #if line has more than one equal sign, is a problem!
        name = getname(name)
        entry = parse2list(entry) #TODO: if list of numbers, cast from string to int/float as appropriate
        if len(entry) == 1:
          entry = entry[0]
        levels[-1][name] = entry
      l = f.readline()
  #print('ending levels: {}'.format(levels))
  #print(len(levels))
  return paramdict(spec)

def load_operator(dirname='.'):
  df = pd.read_csv(dirname+'/operators.dat',delimiter='\s+',escapechar='#',error_bad_lines=False) 
  return df

class paramdict(OrderedDict):
  '''To enable lazy, case-insensitive access'''
  def iget(self,key):
    res = list(self.nested_lookup(self,key))
    if len(res) == 0:
      raise KeyError(key)
    if len(res) == 1:
      return res[0][0][res[0][1]]
    if len(res) > 1:
      print('Multiple ({}) entries found!'.format(len(res),res))
      extractedresults = [ r[0][r[1]] for r in res ] #unsure if this is what I want for output. more succinct, but less usable 
      return res
  def iset(self,key,val,setall=False):
    res = list(self.nested_lookup(self,key))
    if len(res) == 0:
      self[key] = val
    if len(res) == 1:
      res[0][0][res[0][1]] = val
    if len(res) > 1:
      if setall:
        print('chose to overwrite multiple ({}) entries found'.format(len(res)))
        for r in res:
          r[0][r[1]] = val
      else:
        raise KeyError('Multiple ({}) entries found, unclear which to set: {}'.format(len(res),res))

  @staticmethod
  #returns recursively subnested values
  #see gen_dict_extract from
  def nested_lookup(d,key): 
    '''
    may be slow for large dicts: nested, and in py2.7 creates full lists instead of using iterators when iterating.
    '''
    if hasattr(d,'iteritems'):
      for k,v in d.items():
        if k.lower() == key.lower():
          #print('found {}'.format(d))
          yield (d,k)
        if isinstance(v, dict):
          for result in paramdict.nested_lookup(v,key):
            yield result
        elif isinstance(v, (list,tuple)):
          for e in v:
            for result in paramdict.nested_lookup(e,key):
              yield result

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
    elif isinstance(val,float):
      lines.append(key.ljust(21) + ' = {:.15f}'.format(val))
    elif isinstance(val,(int,str)):
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

### doing simple substitutions
# fill in later; temporary work-around is to use the generated yaml file, which can be easily written out in pfts format! Maybe more verbose, but also more intuitive and direct, without messing with regexp, can easily add new sections, etc.
# i.e. functions to manage an input script, and all the different sections

def setkey(d,key,val,remove=False):
  if val is not None:
    d[key] = val
  if remove and val is None:
    d.pop(key,None)

def set_cell(spec,Ls=None,CellScaling=None,CellAngles=None,NPW=None,SpaceGroup=None,CenterToPrimitive=None,Symmetrize=None,remove=False):
  '''
  remove : bool
    whether or not to remove fields that have "None"
  '''
  subdict = spec['Models']['Model1']['Cell']
  if Ls is not None:
    if isinstance(Ls,(float,int)):
      Dim = 1
    elif isinstance(Ls,(tuple,list,np.array)):
      Dim = len(Ls)
    else:
      Dim = subdict['Dim']
  else:
    Dim = subdict['Dim']

  setkey(subdict,'Dim',Dim,remove=False)
  setkey(subdict,'CellScaling',CellScaling,remove=remove)
  setkey(subdict,'CellLengths',Ls,remove=False)
  setkey(subdict,'CellAngles',CellAngles,remove=remove)
  setkey(subdict,'NPW',NPW,remove=False)
  setkey(subdict,'SpaceGroupName',SpaceGroup,remove=remove)
  setkey(subdict,'CenterToPrimitiveCell',CenterToPrimitive,remove=remove)
  if isinstance(Symmetrize,str) and Symmetrize.lower() == 'on':
    Symmetrize = 'on'
  elif isinstance(Symmetrize,bool) and Symmetrize:
    Symmetrize = 'on'
  elif Symmetrize is not None:
    Symmetrize = 'off'
  setkey(subdict,'Symmetrize',Symmetrize,remove)

def set1d(spec,L=None,NPW=None):
  if L is not None and not isinstance(L,(float,int)) and len(L)!=1:
    raise ValueError('1d system must have box L be one dimension instead of {}'.format(L))

  set_cell(spec,Ls=[1.0],CellScaling=L,CellAngles=[90.],NPW=NPW,
      SpaceGroup=None,CenterToPrimitive=None,Symmetrize=None,remove=True)

def set2d(spec,Ls=None,NPWs=None,Angle=90.,HEX=True,Symmetrize=False):
  if Ls is not None:
    if isinstance(Ls,(float,int)) or len(Ls)!=2:
      raise ValueError('2d system must have box L be 2d instead of {}'.format(Ls))
  if isinstance(NPWs,(int,float)):
    NPWs = [int(NPWs)]*2
  if HEX:
    Angle = 120.
    SpaceGroup = 'p6mm'
    CenterToPrimitive = True
  else:
    Angle = 90.
    SpaceGroup = None
    CenterToPrimitive = None
    Symmetrize = None
  set_cell(spec,Ls,1.0,Angle,NPWs,SpaceGroup,CenterToPrimitive,Symmetrize,remove=True)

def set3d(spec,Ls=None,NPWs=None,Angles=None,SpaceGroup=None,Symmetrize=False):
  """still incomplete logic, i.e. sometimes want to override, sometimes dont...
  safest is to always specify SpaceGroupName, which will then update angles accordingly"""
  if Ls is not None:
    if isinstance(Ls,(float,int)) or len(Ls)!=3:
      raise ValueError('3d system must have box L be 3d instead of {}'.format(Ls))
  if isinstance(NPWs,(int,float)):
    NPWs = [int(NPWs)]*3
  if SpaceGroup is not None:
    if SpaceGroup.lower().startswith('i'):
      Angles = [109.4712206345, 109.4712206345, 109.4712206345] #arccos(-1/3), for the primitive cell
    else:
      Angles = [90.,90.,90.]# will be overriden by CenterToPrimitive anyway
    CenterToPrimitive = True
  else:
    CenterToPrimitive = False
  remove = False
  if SpaceGroup is False:
    SpaceGroup = None
    CenterToPrimitive = None 
    Symmetrize = None
    remove = True

  set_cell(spec,Ls,1.0,Angles,NPWs,SpaceGroup,CenterToPrimitive,Symmetrize,remove=remove)


# and some setters that don't do any validation, relies on getting right key name:
# no key removals, need to set everything correctly!!!
def set_cell_simple(spec,**kwargs):
  subdict = spec['Models']['Model1']['Cell'] 
  for key in kwargs:
    update(subdict,kwargs,getname(key))

def set_operators(spec,**kwargs):
  subdict = spec['Models']['Model1']['Operators'] 
  for key in kwargs:
    update(subdict,kwargs,getname(key))

def set_init(spec,**kwargs):
  subdict = spec['Models']['Model1']['InitFields'] 
  for key in kwargs:
    update(subdict,kwargs,getname(key))

def set_simulation(spec,**kwargs):
  subdict = spec['Simulation'] 
  for key in kwargs:
    update(subdict,kwargs,getname(key))

def set_pll(spec,**kwargs):
  subdict = spec['Parallel'] 
  for key in kwargs:
    update(subdict,kwargs,getname(key))


# More custom methods
def set_init_file(spec,filename):
  if filename in [None,False]:
    set_init(spec,ReadInputFields='False')
  elif filename == True:
    set_init(spec,ReadInputFields='HFields')
  elif isinstance(filename,str):
    set_init(spec,ReadInputFields='HFields',InputFieldsFile=filename)
  else:
    raise ValueError('Unrecognized file specification for initial fields: {}'.format(filename))


