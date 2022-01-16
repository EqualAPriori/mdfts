"""wrapper and helper functions to facilitate use of ruamel.yaml

By default, also outputs an equivalent .json file

Typical usage example:
    import yamlhelper as yml
    yml.load()
    yml.save_dict('filename',mydict)

Examples:
    >>> save_dict('test.yaml',{'asdf':1,'qwer':2},header='testheader')
    >>> load('test.yaml')
    ordereddict([('asdf', 1), ('qwer', 2)])
    >>> type(load('test.yaml'))
    <class 'ruamel.yaml.comments.CommentedMap'>

Todo:
    Consider adding representers for numpy floats, ints, arrays. may have to be careful about bit version. dig into Representer code a bit more.
"""
# Standard Imports
import os, json
from collections import OrderedDict

# 3rd party Imports
import ruamel.yaml as YAML


def create_yaml():
    """Create YAML object with default settings

    Returns:
        ruamel.yaml.YAML(): 
            a ruamel.yaml YAML() object w/ utf-8 encoding, unicode encoding, OrderedDict representation

    Note:
        https://stackoverflow.com/questions/49669236/ruamel-yaml-bad-dump
    """
    yaml = YAML.YAML()
    yaml.explicit_start = True
    yaml.default_flow_style = None
    yaml.encoding = "utf-8"  # default when using YAML() or YAML(typ="rt")
    yaml.allow_unicode = True  # always default in the new API
    yaml.errors = "strict"
    yaml.Representer.add_representer(OrderedDict, yaml.Representer.represent_dict)
    return yaml


_yaml = create_yaml()  # private default yaml object for loading and writing files


def save_dict(filename, mydict, header=None):
    """ Saves any yaml-writable object (e.g. list, dict-like, or yaml representation) to .json and .yaml

    Args:
        filename (str): filename for yaml file, preferably with '.yml' or '.yaml' extension
        mydict (yaml-representable): list, dict, str, or yaml representation
        header (str, optional): header comment for yaml file. Defaults to None.
    """
    with open(filename, "w") as f:
        f.write("# {}\n".format(header))
        _yaml.dump(mydict, f)
    prefix, ext = os.path.splitext(filename)
    with open(prefix + ".json", "w") as f:
        json.dump(mydict, f, indent=4)


def load(filename):
    """loads a filename

    Args:
        filename (str): file name, NOT a `file` object

    Returns:
        ruamel commented map: behaves like an ordered dictionary

    Todo:
        type-checking to allow feeding in a `file` object
        in fact, ruamel.yaml.YAML objects can also load a string whose contents is in a yaml spec!
    """
    with open(filename, "r") as stream:
        contents = _yaml.load(stream)
    return contents


### === TESTS ===
# test_string = "# testheader\n--- {asdf: 1, qwer: 2}"
# test_dict = OrderedDict([('asdf', 1), ('qwer', 2)])

