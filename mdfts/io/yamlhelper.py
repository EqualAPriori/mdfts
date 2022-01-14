import ruamel.yaml as YAML
import json
from collections import OrderedDict
# Consider adding representers for numpy floats, ints, arrays. may have to be careful about bit version. dig into Representer code a bit more.

def create_yaml():
    '''
    Notes
    -----
    https://stackoverflow.com/questions/49669236/ruamel-yaml-bad-dump
    '''
    yaml = YAML.YAML()
    yaml.explicit_start = True
    yaml.default_flow_style = None 
    yaml.encoding = "utf-8"     # default when using YAML() or YAML(typ="rt")
    yaml.allow_unicode = True   # always default in the new API
    yaml.errors = "strict"
    yaml.Representer.add_representer(OrderedDict, yaml.Representer.represent_dict)
    return yaml

yaml = create_yaml()

def save_dict( filename, mydict, header=None ):
    import os
    with open( filename, 'w' ) as f:
        f.write('# {}\n'.format(header))
        yaml.dump( mydict, f )
    prefix,ext = os.path.splitext(filename)
    with open( prefix + '.json', 'w' ) as f:
      json.dump( mydict, f, indent=4 )

def load( filename ):
    with open(filename,'r') as stream:
        contents = yaml.load(stream)
    return contents
