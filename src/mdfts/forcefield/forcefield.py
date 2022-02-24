"""
forcefield.py: Constructs a force field for a system. Contains bead types and
potentials between those bead types.
"""
from __future__ import absolute_import, division, print_function

from ast import literal_eval
from collections import defaultdict, OrderedDict
import inspect

from .potential import _Potential
from .potential_gaussian import Gaussian
from .potential_harmonic_bond import HarmonicBond
from .beadtype import BeadType
from mdfts.utils import serial
from mdfts.utils import yamlhelper

__all__ = ["load_from_sim_ff", "ForceField"]


# collect all classes of potentials and put them into an dictionary
_POTENTIAL_TYPES = {}
_key, _val = None, None
for _key, _val in locals().items():
    if inspect.isclass(_val) and issubclass(_val, _Potential):
        _POTENTIAL_TYPES[_key] = _val


class _SimPotentialSpecification(object):
    """Container that parses and stores information for a potential from a sim
    forcefield file"""

    def __init__(self):
        self.type = None
        self.bead_names = None
        self.parameters = OrderedDict()

    def from_string(self, s):
        # clean up string
        s = s.replace(">>> POTENTIAL", "")
        s = s.replace("\n", "")

        # split into name and data
        _potential_name, _potential_data = s.split("{", 1)
        _potential_name = _potential_name.strip()
        _potential_data = "{" + _potential_data
        _potential_data = _potential_data.strip()

        # determine potential type and bead names
        _potential_name_split = _potential_name.split("_")
        _potential_type = _potential_name_split[0]
        try:
            self.type = _POTENTIAL_TYPES[_potential_type]
        except KeyError:
            raise ValueError("{} is not a valid potential type".format(_potential_type))
        self.bead_names = _potential_name_split[1:]

        # set values of parameters
        self.parameters = literal_eval(_potential_data)


def load_from_sim_ff(filepath, kT=1.0):
    """Initializes a ForceField instance from a sim force field file"""
    # initialize ForceField instance
    ff = ForceField(kT=kT)

    # open sim force field file
    s = open(filepath, "r").read()

    # split string by potential
    potential_strings = [">>> POTENTIAL" + p for p in s.split(">>> POTENTIAL")[1:]]

    # separate out potentials by type
    potential_specification_dict = defaultdict(list)
    for ps in potential_strings:
        potential_specification = _SimPotentialSpecification()
        potential_specification.from_string(ps)
        potential_specification_dict[potential_specification.type].append(
            potential_specification
        )

    # check if there are Gaussian interactions; if there are, separate the like
    # Gaussian interactions from all the other ones
    if Gaussian in potential_specification_dict.keys():

        # separate Gaussian interactions between like and unlike
        like_gaussians = []
        other_gaussians = []
        for g in potential_specification_dict[Gaussian]:
            if g.bead_names[0] == g.bead_names[1]:
                like_gaussians.append(g)
            else:
                other_gaussians.append(g)

        # use like Gaussian interactions to compute smearing lengths
        for g in like_gaussians:
            bead_name = g.bead_names[0]
            if not ff.has_bead_type(bead_name):
                ff.add_bead_type(BeadType(bead_name))
            bead_type = ff.get_bead_type(bead_name)
            gaussian = Gaussian(bead_type, bead_type)
            gaussian.from_sim_specification(g, kT=kT)
            ff.add_potential(gaussian)

        # add other Gaussian interactions
        for g in other_gaussians:
            bead_type_1 = ff.get_bead_type(g.bead_names[0])
            bead_type_2 = ff.get_bead_type(g.bead_names[1])
            gaussian = Gaussian(bead_type_1, bead_type_2)
            gaussian.from_sim_specification(g, kT=kT)
            ff.add_potential(gaussian)

    # add all other interactions
    for p_type, p_spec_list in potential_specification_dict.items():
        if p_type is not Gaussian:
            for p_spec in p_spec_list:
                potential = p_type(*[ff.get_bead_type(bn) for bn in p_spec.bead_names])
                potential.from_sim_specification(p_spec, kT=kT)
                ff.add_potential(potential)

    return ff


@serial.serialize(["_kT", "_bead_types", "_potentials"])
class ForceField(object):
    """Container object for a force field

    A ForceField represents the interactions in a system. The ForceField stores
    BeadTypes and Potentials between those bead types. BeadTypes and
    Potentials can be added to the force field through various methods.
    """

    def __init__(self, ffname=None, kT=1.0):
        """Constructor for the ForceField class"""
        self.kT = kT
        self._bead_types = serial.SerializableTypedList(BeadType)
        self._potentials = serial.SerializableTypedDict()
        if isinstance(ffname, str):
            ffdef = yamlhelper.load(ffname)
            self.from_dict(ffdef)

    @property
    def kT(self):
        """Thermal energy of the system. Used when outputting parameters of
        potentials in terms of real energy units."""
        return self._kT

    @kT.setter
    def kT(self, value):
        """Set the thermal energy of the system"""
        try:
            self._kT = float(value)
        except ValueError:
            raise ValueError("must be able to convert new value of kT to float")

    @property
    def bead_types(self):
        """List of BeadTypes in the ForceField"""
        return list(self._bead_types)

    @property
    def bead_names(self):
        """List of bead names of the ForceField"""
        return [bt.name for bt in self.bead_types]

    @property
    def bead_type_dict(self):
        return {bt.name: bt for bt in self.bead_types}

    def add_bead_type(self, bead_type):
        """Add a BeadType to the ForceField"""
        # check that ForceField doesn't already contain a BeadType of the same name
        if bead_type.name in self.bead_names:
            raise ValueError(
                "ForceField instance already contains BeadType '{}'".format(
                    bead_type.name
                )
            )

        self._bead_types.append(bead_type)

    def has_bead_type(self, bead_name):
        """Check if BeadType of the specified name exists in the ForceField"""
        try:
            self.get_bead_type(bead_name)
            return True
        except ValueError:
            return False

    def get_bead_type(self, bead_name):
        """Returns a BeadType of the specified name"""
        # return bead type of specified name if it exists in ForceField
        for bt in self.bead_types:
            if bead_name == bt.name:
                return bt

        # raise error if BeadType with specified name doesn't exist in ForceField
        raise ValueError(
            "ForceField instance does not contain BeadType '{}'".format(bead_name)
        )

    def reorder_bead_types(self, reordered_bead_names):
        """Rearranges the order of BeadTypes in the ForceField. Useful when
        converting to PolyFTS format."""
        if set(reordered_bead_names) != set(self.bead_names):
            raise ValueError(
                "provided bead names don't match bead names in ForceField instance"
            )
        self._bead_types = serial.SerializableTypedList(
            BeadType, *[self.get_bead_type(bn) for bn in reordered_bead_names]
        )

    @property
    def potentials(self):
        """List of _Potentials in the ForceField"""
        _potential_list = []
        for v in self._potentials.values():
            _potential_list.extend(v)
        return _potential_list

    def add_potential(self, potential):
        """Add a _Potential to the ForceField"""
        # check that the potential inherits from base _Potential class
        if not isinstance(potential, _Potential):
            raise TypeError("potential must be of type _Potential or inherit from it")

        # create entry for potential type if it doesn't exist already
        potential_class = potential.__class__
        potential_class_name = potential_class.__name__
        if potential_class_name not in self._potentials.keys():
            self._potentials.add_entry_type(
                potential_class_name, potential_class, has_many=True
            )

        # add potential
        # TODO: add warning if _Potential using same BeadTypes already exists in ForceField
        self._potentials[potential_class_name].append(potential)

    def get_potential(self, potential_type, *args):
        if inspect.isclass(potential_type) and issubclass(potential_type, _Potential):
            potential_type = potential_type.__name__
        for p in self._potentials[potential_type]:
            if set(args) == set(p.bead_names):
                return p
        raise ValueError(
            ("no potential of type '{}' with following bead types:" + " {},").format(
                potential_type.__name__, *args
            )
        )

    def __str__(self):
        """Default method is overridden to allow for more transparent printing
        for debugging purposes"""
        s = "\n"
        s += "BeadTypes: {}\n".format([bt.__str__() for bt in self.bead_types])
        s += "Potentials: {}".format([p.__str__() for p in self.potentials])
        return s

    # ===== Helper methods for serialization
    def from_dict(self, d):
        self.infer_potential_schema(d)  # infer the potential schema

        serial.Serializable.from_dict(self, d)

        for ps in self._potentials.values():
            for p in ps:
                p._bead_types.align_to_dict(self.bead_type_dict)

    def infer_potential_schema(self, d):
        potentials = d["_potentials"]
        for k in potentials:
            self._potentials.add_entry_type(k, _POTENTIAL_TYPES[k], has_many=True)

    def save(self, filename="ff.yaml"):
        yamlhelper.save_dict(filename, self.to_dict())

    def load(self, filename="ff.yaml"):
        ffdef = yamlhelper.load(filename)
        self.from_dict(ffdef)
