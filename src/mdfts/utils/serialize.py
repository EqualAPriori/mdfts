""" Helper tools to create objects that automatically serialize and deserialize
By serializing to dictionary, can trivially read/write with (json).
(Provided that the dictionary only has basic data types, we don't do type checking here)

Note:
    I envision simply needing to have classes subclass Serializable and decorate...
    and then use decorators to track which variables need serialization/deseriazliation
        and also track complexity of the objects being serialized 

Todo:
    - decorators to track variables that can be saved to and read from dictionaries.
    - different use-methods:
        - subclass Serializable, and define _serial_vars and _extra_vars accordingly
        - use metaclass as aid to ensure _serial_vars is indeed defined
        - decorator that takes _serial_vars as argument
"""
from collections import OrderedDict
import collections

try:
    collectionsABC = collections.abc
except:
    collectionsABC = collections

import inspect


# def serialize(tracked_args):


class Meta(type):
    """to inject default class variables into a class

    Can also use to change the __str__ and __repr__ of a class.
    AH, e.g. I can use the metaclass in a decorator!
    """

    def __new__(cls, clsname, bases, clsdict):
        c = super().__new__(cls, clsname, bases, clsdict)
        c._serial_vars = []
        c._extra_vars = None
        # too clunky... will require classes to implement/define __clsname__
        # c.__clsname__ = ""
        # .__clsmodule__ = ""
        return c


class Serializable(object):
    """ Serializable object

    Note:
        serializes from and serializes *to* an OrderedDict
        Should be able to handle general objects that need to get serialized.
        Can override the to_dict and from_dict functions, e.g. to handle special formatting.
        (esp. for ForceField... think about if abbreviations should be handled by the object or by the parser.)

        In order to initialize from class, the dict either needs to contain enough 
        input variables to initialize (and no extra variables),
        or last reort the class constuctor must be able to take no arguments (i.e. all optional).

    Todo:
        support creating from lists, although this is not ideal since one will rely on 
        implicit knowledge of ordering of variables, an implementation detail.

        also: allow for lazy key matching

        __init_from_dict not tested/working yet... 
        or at least requires certain naming conventions of constructor for it to work out.
    """

    _serial_vars = []  # names of variables. should prob. default to self.__dict__...
    _extra_vars = None  # variables that were fed in but not expected for serialization

    def __init__(self):
        # If one makes a bare Serializable object, then need to make sure that it
        # has instance-level tracking of variables
        self._serial_vars = []  # MUST IMPLEMENT IN SUBCLASS, e.g. not as instance var
        self._extra_vars = None  # MUST IMPLEMENT IN SUBCLASS
        pass

    def to_dict(self):
        """always send out a dict
        
        Essentially, this dict is the *state* of the system!"""
        od = OrderedDict()
        for k in self._serial_vars:
            od[k] = self.custom_get(k)  # custom handle (int,float,str,bool,list,dict)?
        return od

    def from_dict(self, d):
        """set state/variables from dict
        
        Note:
            can also support list, but that is not preferable since it relies on knowing
            the sequence of tracked variables in the class.
        """
        print(self)
        if isinstance(d, collectionsABC.Mapping):  # dictlike
            for (k, v) in d.items():
                if k in self._serial_vars:
                    self.custom_set(k, v)
                else:
                    if self._extra_vars is None:
                        self._extra_vars = OrderedDict()
                        self._extra_vars[k] = v
        elif isinstance(d, collectionsABC.Iterable):  # list-like
            # assume is the tracked variables, in order. not need to be complete list.
            # in this case d stores the values
            # really not recommended to use this... since it requires detailed knowledge of implementation
            for iv, v in enumerate(d):
                if iv < len(self._extra_vars):
                    k = self._extra_vars[iv]
                    self.custom_set(k, v)
                else:
                    self._extra_vars = d[iv:]
        else:
            raise ValueError(
                "Unsupported input data type {}, can't serialize".format(type(d))
            )

    def custom_get(self, k):
        if isinstance(getattr(self, k), Serializable):
            return getattr(self, k).to_dict()
        else:
            return getattr(self, k)  # retrieve the getter/variable of choice

    def custom_set(self, k, v):
        # if serial_vars[k] points to a Serializable object,
        # then need to recursively update that object's values by calling from_dict!!
        if isinstance(getattr(self, k), Serializable):
            getattr(self, k).from_dict(v)
        else:  # type of self.k should be int,float,str,bool,list,tuple,something standard, and so should `v`
            setattr(self, k, v)
            # if self.k is a non-serializable, non-standard object...
            # then this could be problematic. we don't check for such eventualities here.

    @classmethod
    def init_from_dict(cls, d, *args, **kwargs):
        """ 1) Try basic initialization. This is the case where the dict doesn't carry non-constructor parameters.
        Todo:
            Not working yet
        """
        try:
            obj = cls(*args, **kwargs)
        except:
            print("Could not initialize {}".format(cls))

        # final update
        obj.from_dict(d)


def serialize(track_args):
    """Decorator that makes track_args/_serial_vars EXPLICITLY required.

    Args:
        track_args (list, tuple, Iterable): [description]

    Returns:
        cls: decorated class

    Note:
        Can make the str() and repr() modifications by decorator as well!
    """

    def deco(cls):
        old_init = cls.__init__

        class tmpMeta(type):
            """to change the __str__ and __repr__ a wrapped class uses to report itself.
            """

            def __str__(self):
                return "<serialize-wrapped class '{}.{}'>".format(
                    cls.__module__, cls.__name__,
                )

            def __repr__(self):
                return "<class '{}.{}-wrapped {}.{}'>".format(
                    serialize.__module__,
                    serialize.__name__,
                    cls.__module__,
                    cls.__name__,
                )

        class myclass(cls, Serializable, metaclass=tmpMeta):
            """
            demonstrates proper subclassing of Serializable
            defines new _serial/extra_vars class variables

            Todo:
                have to modify metaclass syntax for python 2.x
            """

            _serial_vars = track_args  # class variables
            _extra_vars = None

            def __init__(self, *args, **kwargs):
                old_init(self, *args, **kwargs)
                self.__class__.__name__ = cls.__name__  # retain info after decorating
                self.__class__.__module__ = cls.__module__

            def __str__(self):  # pretier reporting to get around manual decorator
                return "<{}.{} object at {}>".format(
                    self.__class__.__module__, self.__class__.__name__, hex(id(self))
                )

            def __repr__(self):  # pretier reporting to get around manual decorator
                return "<{}.{}-wrapped {}.{} object at {}>".format(
                    serialize.__module__,
                    serialize.__name__,
                    self.__class__.__module__,
                    self.__class__.__name__,
                    hex(id(self)),
                )

        return myclass

    return deco
