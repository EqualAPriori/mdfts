#! My parsing helper functions
# (C) Kevin Shen, 2021
import re
import copy
import os

# ===== Pathing =====
def findpath(fname, paths):
    """Return path to first instance of path/fname, from given list of paths

  Args:
      fname (str): file to search for
      paths (str or list): paths to search through

  Raises:
      ValueError: if file not found in paths

  Returns:
      str: first instance of path/fname
  """
    if paths is None:
        paths = ""
    if isinstance(paths, str):  # assume none of [\s,:;] in the path name
        paths = re.split(r"[\s,:;]+", paths)
    for path in paths:
        fname_full = os.path.abspath(path + "/" + fname)
        if os.path.exists(fname_full):
            print("located file {} at path {}: {}".format(fname, path, fname_full))
            return fname_full
    raise ValueError("file {} not found in paths {}".format(fname, paths))


# ===== Atomic Parsing Functions =====
def isfloat(value):
    """Check if value (string) can be cast to float. Right now clunky cuz float turns alphabet strings to float!"""
    if type(value) is bool:
        return False
    try:
        float(value)
        return True
    except:
        return False


def isbool(value):
    """Check if value (string) is true/false"""
    if isinstance(value, str) and value.lower() in ["false", "true", "fixed", "free"]:
        return True
    else:
        return False


def tobool(value):
    """Check if value (string) is true/false. Else return original value"""
    if isinstance(value, str) and value.lower() in ["false", "free"]:
        return False
    elif isinstance(value, str) and value.lower() in ["true", "fixed"]:
        return True
    else:
        return value


def parse_beadtypes(entry, delim=";"):
    """Parses an entry for a bead type filter. 

    Args:
        entry (str or list): of bead info
        delim (str, optional): Defaults to ";".

    Raises:
        ValueError: if beadentry is not supported

    Returns:
        list: bead entry rectified to guarantee it's a list
    """
    if isinstance(entry, str):
        # result = entry.split(delim)
        result = re.split(r"[,:{}]+".format(delim), entry)
    elif isinstance(entry, list):
        result = entry
    else:
        raise ValueError("Unrecognized beadtype entry {}".format(entry))

    return result


def parse_entry(entry, delim=";"):
    """Parses an entry for a parameter

    Args:
        entry (str, list, dict): data
        delim (str, optional): [description]. Defaults to ";".

    Returns:
        processed_entry : dict or value
          namestring (anything that can't be split or cast to float or bool)--> string
          float --> {'val:float}, bool --> {'fixed:bool}, namedparam set --> {param: name or {'val': float, 'fixed: bool}}
    """

    """
    Notes
    -----
    "name" is the only special field that is left untouched
    as constructed, don't worry about iterative parsing. trying to keep structures flat, parse_line handles multiple entries. Each entry will either be some key:value pair or only contain values

    Different use cases (non-exhaustive):
    entry = 
    - "name:somename"
    - "name;somename"
    - "paramname;float;bool"
    - "float"
    - "bool"
    - [paramname,float,bool]
    - [paramname,"float;bool"]
    - [paramname,"float bool"]
    - {paramname: "float;bool"}
    - {paramname: [float,bool]}
    - {paramname: {val:, fixed:}} #don't need to do anything

    i.e. entry should be atomic, single parameter, not multiple parameters in a container

    Examples
    --------
    >>> parsify.parse_beadtypes('A;B;CC')
    ['A', 'B', 'CC']

    >>> parsify.parse_entry(1.0)
    {'val': 1.0}

    >>> parsify.parse_entry(False)
    {'fixed': False}
    
    >>> parsify.parse_entry([1.0,False])
    {'fixed': False, 'val': 1.0}

    >>> parsify.parse_entry(['B',1.0,False])
    {'B': {'fixed': False, 'val': 1.0}}

    >>> parsify.parse_entry(['B',[1.0,False]])
    {'B': {'fixed': False, 'val': 1.0}}

    >>> parsify.parse_entry(['B','1.0;False'])
    {'B': {'fixed': False, 'val': 1.0}}

    >>> parsify.parse_entry(['B','1.0 False']) #also supports white space delimiting, if appropriate
    {'B': {'fixed': False, 'val': 1.0}}

    parsify.parse_entry({'name':'lala'})
    {'name': 'lala'}

    >>> parsify.parse_entry('name;lala')
    {'name': 'lala'}

    >>> parsify.parse_entry({'B':"1.0 False"})
    {'B': {'fixed': False, 'val': 1.0}}

    >>> parsify.parse_entry({'B':"1.0;False"})
    {'B': {'fixed': False, 'val': 1.0}}

    >>> parsify.parse_entry({'B':[1.0,'False']})
    {'B': {'fixed': False, 'val': 1.0}}

    >>> parsify.parse_entry({'B':{'val':1.0, 'fixed':False}})
    {'B': {'fixed': False, 'val': 1.0}}
    """
    # print('received {}'.format(entry))
    if isinstance(entry, list) and len(entry) == 1:
        entry = entry[0]
    if isinstance(entry, str):
        # entry = [ e.strip() for e in entry.split(delim) ]
        entry = [e.strip() for e in re.split(r"[\s,:{}]+".format(delim), entry)]
    if isinstance(entry, float) or isinstance(entry, bool) or isinstance(entry, int):
        entry = [entry]

    processed_entry = {}
    if isinstance(entry, list):
        for ie, el in enumerate(entry):
            if isbool(el):
                entry[ie] = tobool(el)
            elif isinstance(el, int):
                entry[ie] = el
            elif isfloat(el):
                entry[ie] = float(el)
            else:
                entry[ie] = el

        if isinstance(entry[0], str) and entry[0].lower() in ["name"]:
            processed_entry = {"name": entry[1]}
        elif isinstance(entry[0], str):
            if (
                len(entry) == 1
            ):  # envision this being the case where just getting the name value
                return entry[0]
            else:  # this is the case when we have a key:value pair
                paramname = entry[0]
                field = entry[1:]
                processed_entry = {paramname: parse_entry(field)}
        else:  # this is the case when we just get a value
            for el in entry:
                if isinstance(el, bool):
                    processed_entry["fixed"] = el
                elif isinstance(el, int):
                    processed_entry["val"] = el
                elif isinstance(el, float):
                    processed_entry["val"] = el
        entry = processed_entry

    if isinstance(entry, dict):
        if "val" in entry or "fixed" in entry or "name" in entry:
            # assume is inner-most type definition, no more validation
            processed_entry = entry
        else:
            for k, v in entry.items():
                processed_entry[k] = parse_entry(v)
    # print(processed_entry)
    return processed_entry


def parse_potential_entry(entry, nbody, store_dict=None, prefix=""):
    """Parse a line entry for a potential to create
    Parameters
    ----------
    entry : str,list,dict
    nbody : # atoms involved to specify filters for this potential

    Notes
    -----
    This only does the initial shorthand conversions.
    String format: one line entry of everything.
    List: MUST have format [ #body1, ..., #body-n, strings, lists, or dicts, per parameter ], i.e. does not make sense to have a sublist with multiple interactions
        or [bead1, bead2, {dict of all parameters}]
        i.e. if more than one dict, then must have all entries collected together!
        expect each individual entry to be:
        paramname;val
        [paramname, valuestring]
        [paramname, valuelist]
        [paramname, valuedict]
        {paramname: valuestring, valuelist, valuedict}

    try to use storage container, i.e. if storage is Ordered (e.g. Yaml object), then hopefully the resulting output will be as well.

    Examples
    --------
    Vanilla example with string definition of potential 

    >>> parsify.parse_potential_entry('A B name;123 B;1.0;Fixed', 2, prefix='ljg')
    {'B': {'fixed': True, 'val': 1.0}, 'species': ['A', 'B'], 'name': 123.0}

    #auto-named!
    
    >>> parsify.parse_potential_entry('A B B;1.0;Fixed', 2, prefix='ljg')
    {'B': {'fixed': True, 'val': 1.0}, 'species': ['A', 'B'], 'name': 'ljg_A_B'}


    #multiple parameters!
    
    >>> parsify.parse_potential_entry('A B B;1.0;Fixed Kappa;1.0;False', 2, prefix='ljg')
    {'B': {'fixed': True, 'val': 1.0}, 'Kappa': {'fixed': False, 'val': 1.0}, 'species': ['A', 'B'], 'name': 'ljg_A_B'}

    #list as input, with mix of parameter definition styles

    >>> parsify.parse_potential_entry(['A','B;C','Kappa;1.0;False',['B',1.0,True],{'Dist0': '1.0;fixed'}], 2, prefix='ljg')
    process line
    {'B': {'fixed': True, 'val': 1.0}, 'Kappa': {'fixed': False, 'val': 1.0}, 'species': [['A'], ['B', 'C']], 'name': 'ljg_A_B;C', 'Dist0': {'fixed': True, 'val': 1.0}}

    #dict as input, with mix of paramter definition styles

    >>> parsify.parse_potential_entry({'species':'A B;C', 'Kappa':[1.0,False], 'B':'1.0 fixed'}, 2, prefix='ljg')
    {'B': {'fixed': True, 'val': 1.0}, 'Kappa': {'fixed': False, 'val': 1.0}, 'species': [['A'], ['B', 'C']], 'name': 'ljg_A_B;C'}

    #an ordered dict can keep things in order:

    >>> from collections import OrderedDict
    >>> od = OrderedDict()
    >>> parsify.parse_potential_entry({'species':'A B;C', 'Kappa':[1.0,False], 'B':'1.0 fixed'}, 2, prefix='ljg', store_dict=od)
    {'B': {'fixed': True, 'val': 1.0}, 'Kappa': {'fixed': False, 'val': 1.0}, 'species': [['A'], ['B', 'C']], 'name': 'ljg_A_B;C'}
    >>> od
    OrderedDict([('species', [['A'], ['B', 'C']]), ('B', {'fixed': True, 'val': 1.0}), ('Kappa', {'fixed': False, 'val': 1.0})])

    One last example
    >>> parsify.parse_potential_entry(['A','B,D,E', {'Kappa':[1.0,False], 'B':'1.0 fixed'}], 2, prefix='ljg', store_dict=od)
    {'B': {'fixed': True, 'val': 1.0}, 'Kappa': {'fixed': False, 'val': 1.0}, 'species': [['A'], ['B', 'D', 'E']], 'name': 'ljg_A_B;D;E'}
    >>> od
    OrderedDict([('species', [['A'], ['B', 'D', 'E']]), ('B', {'fixed': True, 'val': 1.0}), ('Kappa', {'fixed': False, 'val': 1.0})])
    """
    if isinstance(entry, str):
        # a string format
        entry = entry.split()

    if isinstance(entry, list):
        processed_entry = {"species": entry[:nbody]}

        # process rest:
        if len(entry[nbody:]) == 1 and isinstance(entry[nbody], dict):
            for k, v in entry[nbody].items():
                processed_entry[k] = v
        else:
            for ii, v in enumerate(entry[nbody:]):
                tmp_dict = parse_entry(v)
                for k, v in tmp_dict.items():
                    processed_entry[k] = v
        entry = processed_entry

    if isinstance(entry, dict):
        # process species
        if isinstance(entry["species"], str):
            entry["species"] = entry["species"].split()
        if isinstance(entry["species"], list):
            species = [parse_beadtypes(s) for s in entry["species"]]
            species_strings = [";".join(s) for s in species]
            proposed_name = prefix + "_" + "_".join(species_strings)
            entry["species"] = species
            if store_dict is not None:
                store_dict["species"] = copy.copy(species)
        if "name" not in entry:
            entry["name"] = proposed_name
            if store_dict is not None:
                store_dict["name"] = proposed_name
        for k, v in entry.items():
            if k not in ["species", "name"]:
                # at this point, each sub entry v might be string, list, or dict
                entry[k] = parse_entry(v)
                if store_dict is not None:
                    store_dict[k] = entry[k].copy()

    return entry


# ===== Parse full force field definition, fill in missing values =====

