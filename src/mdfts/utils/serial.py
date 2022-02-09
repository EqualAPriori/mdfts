""" Helper tools to create objects that automatically serialize and deserialize

By serializing to dictionary, can trivially read/write with (e.g. json).
(Provided that the dictionary only has basic data types, we don't do type checking here)

Note:
    I envision simply needing to have classes subclass Serializable...
    or use decorators to track which variables need serialization/deseriazliation

    Note that nesting stops at whenever an object is no longer of type `Serializable`.
    I.e. if a variable has a *plain list* of serializable objects, the serialization will simply save a list of the `__repr__`'s of the objects. 
    Instead, one needs to use a SerializableTypedList to properly continue the nesting and have the Serializable objects handle their own serialization (and deserialization).

    For more heavy duty work, use the SerializableTypedDict as a kind of schema for tracking object types, etc.

Todo:
    - decorators to track variables that can be saved to and read from dictionaries.
    - different use-methods:
        - subclass Serializable, and define _serial_vars and _extra_vars accordingly
        - use metaclass as aid to ensure _serial_vars is indeed defined
        - decorator that takes _serial_vars as argument
    - once we migrate to python 3.4+, implement abstract base classes to ensure interface
    - change terminology, because custom_get and custom_set are more like encoders and decoders in SerializableTypedList and SerializableTypedDict,
        whereas they're more dict-like (non-encoding) in the simple Serializable class. Should probably try to make behavior consistent.
"""
from __future__ import absolute_import, division, print_function

from collections import OrderedDict
import collections
import inspect

try:
    collectionsABC = collections.abc
except:
    collectionsABC = collections

import inspect


# ===== Helpers =====
verbose = False


def vprint(msg):
    if verbose:
        print(msg)


def with_metaclass(mcls):
    """python 2 & 3 compatible metaclass decorator syntax
    
    see: https://stackoverflow.com/questions/22409430/portable-meta-class-between-python2-and-python3
    
    does not handle __slots__ variable that is in the `six` library
    """

    def decorator(cls):
        body = vars(cls).copy()
        # clean out class body
        body.pop("__dict__", None)
        body.pop("__weakref__", None)
        return mcls(cls.__name__, cls.__bases__, body)

    return decorator


# ===== Encoders =====
# special "decoders" and "encoders"
# note, the custom_get(), custom_set(), and init_from() methods in `Serializable`
# allow for modified encoding and decoding


_standard_types = (str, int, float, list, bool, dict)
import numpy as np


def encode(thing):
    """Helpers to encode to .json and .yaml types"""
    if isinstance(thing, _standard_types):
        return thing
    elif isinstance(thing, (np.int64, np.int32)):
        return int(thing)
    elif isinstance(thing, (np.float64, np.float32)):
        return float(thing)
    elif isinstance(thing, (np.array)):
        return list(thing)
    else:
        raise TypeError("unknown type {}, can not encode".format(type(thing)))


def decode(thing, target_type=None):
    """Helpers for use *within* a defined class that knows what types it wants to deserialize to/from.
    In practice, each object may want to handle its own custom encoding/decoding.
    """
    if target_type is None:
        if isinstance(thing, _standard_types):
            return thing
        else:
            raise TypeError("unknown type {}, can not decode".format(type(thing)))
    else:
        try:
            val = target_type(thing)
        except:
            raise TypeError(
                "error trying to create type {} with value {}".format(
                    target_type, thing
                )
            )
        return val


# ===== Main Code =====
class Serializable(object):
    """ Serializable object

    Note:
        deserializes from and serializes *to* an OrderedDict

        Should be able to handle general objects that need to get serialized.

        Can override the to_dict and from_dict functions, e.g. to handle special formatting.
        (esp. for ForceField... think about if abbreviations should be handled by the object or by the parser.)

        In order to initialize from class, the dict either needs to contain enough 
        input variables to initialize (and no extra variables),
        or last resort the class constuctor must be able to take no arguments (i.e. all optional).

        IMPORTANT, to_dict() is ~ supposed to give a flat, json/yaml-ready format. 
        Another potential method is to output actually a dictionary of the tracked variables, without flattening... I.e. maybe the encoding should be left to json. While this produces python-usable dictionaries with objects, it is not technically serializable if, for instance, one of the variables stores a non-basic type and non-serializable type.

    Todo:
        support creating from lists, although this is not ideal since one will rely on 
        implicit knowledge of ordering of variables, an implementation detail.

        also: allow for lazy key matching?

        __init_from_dict not tested/working yet... 
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
        # print(self)
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
                    k = self._serial_vars[iv]
                    self.custom_set(k, v)
                else:
                    self._extra_vars = d[iv:]
        else:
            raise ValueError(
                "Unsupported input data type {}, can't serialize".format(type(d))
            )

    def custom_get(self, k):
        """Functions as an *encoder*"""
        if isinstance(getattr(self, k), Serializable):
            return getattr(self, k).to_dict()
        else:
            return encode(getattr(self, k))  # retrieve the getter/variable of choice

    def custom_set(self, k, v):
        """Functions as a *decoder*"""
        # if serial_vars[k] points to a Serializable object,
        # then need to recursively update that object's values by calling from_dict!!
        if isinstance(getattr(self, k), Serializable):
            getattr(self, k).from_dict(v)
        else:  # type of self.k should be int,float,str,bool,list,tuple,something standard, and so should `v`
            setattr(self, k, v)
            # if self.k is a non-serializable, non-standard object...
            # then this could be problematic. we don't check for such eventualities here.
            # children classes can handle this case later

    def __eq__(self, obj):
        ret = isinstance(obj, self.__class__)
        if not ret:
            return False
        else:
            for k in self._serial_vars:
                if getattr(self, k) != getattr(obj, k):
                    return False
            return True

    @classmethod
    def init_from_dict(cls, d, *args, **kwargs):
        """
        Todo:
            Not tested yet
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

    Todo:
        metaclass conflicts with further mix-ins... temporarily ignore. 
        revisit later to get better naming/reporting for the *class*

    Example:
        >>> @serialize(['tracked1','tracked2','anothervar'])
        >>> class custom_class(object):
        >>>   __init__(self,x):
        >>>     self.tracked1 = x
        >>>     self.tracked2 = x**2
        >>>     self.anothervar = x**4
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

        # class myclass(cls, Serializable, metaclass=tmpMeta):
        # @with_metaclass(tmpMeta)
        class myclass(cls, Serializable):
            """
            demonstrates proper subclassing of Serializable
            defines new _serial/extra_vars class variables
            """

            _serial_vars = track_args  # class variables
            _extra_vars = None

            def __init__(self, *args, **kwargs):
                old_init(self, *args, **kwargs)

            def __str__(self):  # prettier reporting to get around manual decorator
                # ret = "<{}.{} object at {}>\n".format(
                #    self.__class__.__module__, self.__class__.__name__, hex(id(self))
                # )
                # ret += "\n{}".format(cls.__str__(self))
                ret = "<{}> {}".format(self.__class__.__name__, cls.__str__(self))
                return ret

            def __repr__(self):  # pretier reporting to get around manual decorator
                ret = "<{}.{}-wrapped {}.{} object at {}>".format(
                    serialize.__module__,
                    serialize.__name__,
                    self.__class__.__module__,
                    self.__class__.__name__,
                    hex(id(self)),
                )
                return ret

        myclass.__name__ = cls.__name__  # retain info after decorating
        myclass.__module__ = cls.__module__
        return myclass

    return deco


def serialize_list(cls):
    """ Take a `Serializable` classdef and make a Serializable list/container class for it.

    Primary job of this class is to provide standard to_dict, from_dict definitions 
    that know how to instantiate further Serializable objects

    Todo:
        Right now from_dict/list completely overwrites the list. Not sure if there's a more targeted way to handle this.

    Example:
        Let MyClass be a Serializable class. Then, can write:
        `MyClassList = serialize_list(MyClass)`
        And then can further add functionality via `MyClassList.new_method() = ...`

    """

    class myclass(collections.MutableSequence, cls, Serializable):
        """
        demonstrates proper subclassing of Serializable
        defines new _serial/extra_vars class variables
        """

        _serial_vars = []  # not used since we extend list and override to_dict()
        _extra_vars = None
        _content_type = cls  # for reference, what kind of object this class stores

        def __init__(self, *args):
            self.list = list()
            self.extend(list(args))

        def to_dict(self):
            return self.to_list()

        def to_list(self):
            return [el.to_dict() for el in self]

        def from_dict(self, d):
            """default behavior is to override"""
            if isinstance(d, collectionsABC.Mapping):  # dictlike
                d_iterator = d.values()
            elif isinstance(d, collectionsABC.Iterable):
                d_iterator = d
            else:
                raise ValueError("Unsupported input data type {}".format(type(d)))
            self.from_list(d_iterator)

        def from_list(self, l):
            del self[:]  # clear out the list for over-riding
            for v in l:
                if isinstance(v, collectionsABC.Mapping):
                    self.append(cls(**v))
                else:
                    # assume elements are ordered same as arguments of `cls`
                    self.append(cls(*v))

        def custom_get(self, k):
            raise NotImplementedError(
                "Serializable container is a list, not dictionary, does not use key access"
            )

        def custom_set(self, k, v):
            raise NotImplementedError(
                "Serializable container is a list, not dictionary, does not use key access"
            )

        # Defining the MutableSequence interface, to enable all list methods
        # like append, bracket access, etc.
        def check(self, v):
            if not isinstance(v, self._content_type):
                raise TypeError("{} not a valid type: {}".format(v, self._content_type))

        def __len__(self):
            return len(self.list)

        def __getitem__(self, i):
            return self.list[i]

        def __delitem__(self, i):
            del self.list[i]

        def __setitem__(self, i, v):
            self.check(v)
            self.list[i] = v

        def insert(self, i, v):
            self.check(v)
            self.list.insert(i, v)

        # def __str__(self):
        #    return str(self.list)

        def __str__(self):  # prettier reporting
            ret = "<{}_list> {}".format(
                self.__class__.__name__, [str(el) for el in self]
            )
            return ret

        def __repr__(self):
            return super(myclass).__repr__() + " " + self.list.__repr__()

        """
        def __repr__(self):  # pretier reporting to get around manual decorator
            ret = "<{}.{}-wrapped {}.{}_list object at {}>".format(
                serialize.__module__,
                serialize.__name__,
                self.__class__.__module__,
                self.__class__.__name__,
                hex(id(self)),
            )
            return ret
        """

    myclass.__name__ = cls.__name__ + "_list"  # retain name
    myclass.__module__ = cls.__module__
    return myclass


class SerializableTypedList(collections.MutableSequence, Serializable):
    """Create custom list that type-checks its contents.
    
    TypedList, with special handling for use with Serializable classes. 
    For creating a quick list containers when creating a new class is overkill.
    
    Note:
        taken from https://stackoverflow.com/questions/3487434/overriding-append-method-after-inheriting-from-a-python-list
        
        This is just a reference.

        deserialization behavior:
            OVERWRITES internal list completely
    """

    def __init__(self, oktype, *args):
        # if not issubclass(oktype, Serializable):
        #    raise TypeError("list contents must be a Serializable Class")
        self.oktype = oktype
        self.list = list()
        self.extend(list(args))
        self._serial_vars = []  # not used since we extend list and override to_dict()
        self._extra_vars = None

    # === SERIALIZABLE INTERFACE ===
    def to_dict(self):
        return self.to_list()

    def to_list(self):
        if issubclass(self.oktype, Serializable):
            return [el.to_dict() for el in self]
        else:
            # return list(self)
            return [encode(el) for el in self]

    def from_dict(self, d):
        """default behavior is to override"""
        if isinstance(d, collectionsABC.Mapping):  # dictlike
            d_iterator = d.values()
        elif isinstance(d, collectionsABC.Iterable):
            d_iterator = d
        else:
            raise ValueError("Unsupported input data type {}".format(type(d)))
        self.from_list(d_iterator)

    def from_list(self, l):
        del self[:]  # clear out the list for over-riding
        for v in l:
            if isinstance(v, collectionsABC.Mapping):
                self.append(self.oktype(**v))
            else:
                # assume elements are ordered same as arguments of `cls`
                self.append(self.oktype(*v))

    def custom_get(self, k):
        raise NotImplementedError(
            "Serializable container is a list, not dictionary, does not use key access"
        )

    def custom_set(self, k, v):
        raise NotImplementedError(
            "Serializable container is a list, not dictionary, does not use key access"
        )

    def __eq__(self, obj):
        ret = isinstance(obj, self.__class__) and (self.oktype == obj.oktype)
        if not ret:
            return False
        else:
            if len(self) != len(obj):
                return False
            for ik, k in enumerate(self):
                if self[ik] != obj[ik]:
                    return False
            return True

    # === SEQUENCE INTERFACE ===
    def check(self, v):
        if not isinstance(v, self.oktype):
            raise TypeError("{} not a valid type: {}".format(v, self.oktype))

    def __len__(self):
        return len(self.list)

    def __getitem__(self, i):
        return self.list[i]

    def __delitem__(self, i):
        del self.list[i]

    def __setitem__(self, i, v):
        self.check(v)
        self.list[i] = v

    def insert(self, i, v):
        self.check(v)
        self.list.insert(i, v)

    def __str__(self):  # prettier reporting
        ret = "<{}:{}> {}".format(
            self.__class__.__name__, self.oktype.__name__, [str(el) for el in self]
        )
        return ret

    def __repr__(self):
        return (
            super(SerializableTypedDict).__repr__()
            + " for type {} ".format(self.oktype.__name__)
            + self.list.__repr__()
        )


class SerializableTypedDict(collections.MutableMapping, Serializable):
    """Create custom dict that type-checks its contents.

    For creating a quick container when creating a new class is overkill.

    Can be used as a light-weight schema in the vein of Marshmallow.
    If using as a schema, note that *nested* schemas can only be preserved 
    if one defines new classes defining the schemas to be nested. 
    (i.e. if a SerializableTypedDict itself contains a SerializableTypedDict 
    constructed on the fly, the nested SerializableTypedDict's schema won't be preserved).

    Currently, from_dict defaults to *strict mode* (`self._deserialization_strict=True`), 
    only updates entries if they already exist in the dictionary.

    The non-strict mode performs a standard dictionary update (i.e. add new keys).

    Can consider completely over-writing the dictionary upon loading, 
    akin to the implementation in SerializableTypedList.

    Todo:
        Have the SerializableTypedDict preserve its own schema? 
        Or some way to save and load the schema? (essentially, self.keys, self._types, self._has_many)
    
    Note:
        used https://stackoverflow.com/questions/3387691/how-to-perfectly-override-a-dict/47361653#47361653 as reference
    """

    # === MUTABLE MAPPING INTERFACE ===
    def __init__(self, *args, **kwargs):
        self._store = OrderedDict()
        self._types = OrderedDict()
        self._has_many = OrderedDict()
        temp_dict = OrderedDict(*args, **kwargs)
        self.update(temp_dict)  # use the free update to set keys

        self._serial_vars = []  # not used
        self._extra_vars = None
        self._deserialization_strict = True

    def __getitem__(self, key):
        return self._store[(key)]

    def __setitem__(self, key, value):
        if key in self._store:
            # if key exists, then the appropriate type should already be logged.
            self.check(key, value)
        else:
            if isinstance(value, (TypedList, SerializableTypedList)):
                self.add_entry_type(key, value.oktype, has_many=True)

            else:
                self.add_entry_type(key, type(value), has_many=False)

        if isinstance(value, (TypedList, SerializableTypedList)):
            self._store[key] = value
        else:
            self._store[key] = self._types[key](value)

    def __delitem__(self, key):
        del self._store[(key)]

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        return len(self._store)

    def __repr__(self):
        # return super(self.__class__,self).__repr__() + " " + self.list.__repr__()
        ret = "{}({})".format(type(self).__name__, self._store.__repr__())
        return ret

    def __eq__(self, other):
        return Serializable.__eq__(self, other)

    # === HELPERS ===
    def check(self, key, v):
        """check if `v` is of the right type for key `key`
        only works if key in this dictionary!
        I.e. behavior when key not in dictionary should be handled elsewhere, e.g. `__setitem__()`
    
        Args:
            v (instance): should not be a class
        """
        if key not in self._store:  # new key-value pair, assess what type it is
            raise KeyError("key {} not found".format(key))
        elif self._has_many[key] == True:
            if not isinstance(v, SerializableTypedList):
                # if given value is not even a typed list
                raise TypeError(
                    "value for key {} must be a SerializableTypedList for type {}".format(
                        key, self._types[key]
                    )
                )
            else:  # must check if the types agree
                if v.oktype != self._types[key]:
                    raise TypeError(
                        "{} contains type {} that is not the valid type ({}) for key {}".format(
                            v, v.oktype, self._types[key], key
                        )
                    )
        else:
            if isinstance(v, SerializableTypedDict) and issubclass(
                self._types[key], SerializableTypedDict
            ):
                vprint(
                    "Caution! Schema for SerializableTypedDicts may not be saved properly unless the schema for the SerializableTypedDict is defined somewhere in a class."
                )
            elif isinstance(v, SerializableTypedList):
                # not used for now;  This would end up being a generic, un-type-checked SerializableTypedList
                # do not attempt to correct the schema here
                # TODO:
                #   handle this case more cleanly if it arises
                pass
            elif not isinstance(v, (self._types[key], type(None))):
                raise TypeError(
                    "{} not a valid type ({}) for key {}".format(
                        v, self._types[key], key
                    )
                )

        return True

    def add_entry_type(self, key, typ, has_many=None):
        """initialize an appropriate typed entry
        Args:
            has_many (bool): used only if not None

        Notes:
            For use before setting any values, i.e. to define the schema for this dict.
            
            Note that `typ` should be a class, i.e. it can't be a SerializableTypedList instantiated with an oktype. 
         
            I.e. if one intuitively wants to do: 
                add_entry_type(key,SerializableTypedList(typ))  

            one should actually do:

                add_entry_type(key,typ,has_many=True)
        
            because SerializableTypedList(typ) is an instance of SerializableTypedList, and is not a proper type.
            
            The alternative is to use: `mydict['key'] = SerializableTypedList(typ)`
        """
        # For handling the "intuitive but wrong" behavior described in the Notes
        if not inspect.isclass(typ):
            if (
                isinstance(typ, (TypedList, SerializableTypedList)) and has_many is None
            ):  # The only exception
                typ = typ.oktype
                has_many = True
            else:
                raise TypeError(
                    "template entry types can only be added for classes, not an instance like {}".format(
                        typ
                    )
                )

        if has_many:
            self._types[key] = typ
            self._has_many[key] = True
            self._store[key] = SerializableTypedList(typ)
        else:
            self._types[key] = typ
            self._has_many[key] = False
            self._store[key] = None  # default to initialize

    def copy_schema(self, other_dict):
        """Clear this dictionary's schema and copy from another SerializableTypedDict

        Args:
            other_dict (SerializableTypedDict): source schema that we want to copy
        """
        self._store.clear()
        self._types.clear()
        self._has_many.clear()

        for key in other_dict:
            self.add_entry_type(key, other_dict._types[key], other_dict._has_many[key])

    # === SERIALIZABLE INTERFACE ===
    def to_dict(self):
        """always send out a dict
        
        Essentially, this dict is the *state* of the system!"""
        od = OrderedDict()
        for k in self._store:
            od[k] = self.custom_get(k)  # custom handle (int,float,str,bool,list,dict)?
        return od

    def from_dict(self, d):
        """set state/variables from dict
        
        Note:
            deserialization defaults to strict mode. "lax" mode behaves like a normal dictionary update.

        Todo:
            add a mode that is for overwriting entire dict and schema
        """
        # print(self)

        if isinstance(d, collectionsABC.Mapping):  # dictlike
            for (k, v) in d.items():
                if k in self._store:
                    self.custom_set(k, v)
                else:
                    if self._deserialization_strict == True:
                        # Strict, only update values that are in the pre-defined schema (i.e. what is currently in the dict)
                        # This is like the standard Seralizable behavior with `_serial_vars`
                        if self._extra_vars is None:
                            self._extra_vars = OrderedDict()
                            self._extra_vars[k] = v
                    else:
                        # Non-strict mode: add keys that are not present, i.e. the default dictionary update behavior.
                        self.custom_set(k, v)

        elif isinstance(d, collectionsABC.Iterable):  # list-like
            raise ValueError("SerializableTypedDict can not load from iterables")
        else:
            raise ValueError(
                "Unsupported input data type {}, can't serialize".format(type(d))
            )

    def custom_get(self, k):
        if isinstance(self._store[k], Serializable):
            return self._store[k].to_dict()
        else:
            return encode(self._store[k])  # retrieve the getter/variable of choice

    def custom_set(self, k, v):
        # if serial_vars[k] points to a Serializable object,
        # then need to recursively update that object's values by calling from_dict!!
        # note that custom_set here does not type check, because we're deserializing/decoding and the input
        # won't have the right type if it's a non-standard json object type.
        if isinstance(self._store[k], Serializable):
            self._store[k].from_dict(v)
        else:
            self._store[k] = self._types[k](v)


# ===== Unused ======
# def serialize(tracked_args):


class serial_tracking(type):
    """Alternative approach using metaclass to inject default class variables into a class

    Can also use to change the __str__ and __repr__ of a class.
    AH, e.g. I can use the metaclass in a decorator!

    Try later, in case it is cleaner than using decorators. Maybe it's really not that much more handy than simply defining a Serializable interface, and manually inheriting Serializable!
    """

    def __new__(cls, clsname, bases, clsdict):
        c = super(serial_tracking).__new__(cls, clsname, bases, clsdict)
        c._serial_vars = []
        c._extra_vars = None
        # too clunky... will require classes to implement/define __clsname__
        # c.__clsname__ = ""
        # .__clsmodule__ = ""
        return c


class TypedList(collections.MutableSequence):
    """Create custom list that type-checks its contents.
    
    Note:
        taken from https://stackoverflow.com/questions/3487434/overriding-append-method-after-inheriting-from-a-python-list
        This is just a reference; we don't actually use this.

        One potential way to use this is:
        ```
        mylist = TypedList( ASerializableClass ) #returns an *object*
        mylist.to_dict = ...        #i.e. attach methods for the Serializable interface after the fact
        mylist.from_dict = ...      #OR, simply modify this class to contain the Serializable interface

        Current approach is to use the serialize_list() method above, which returns a *class*
        that can have its interface further modified.
        ```
    """

    def __init__(self, oktypes, *args):
        self.oktypes = oktypes
        self.list = list()
        self.extend(list(args))

    def check(self, v):
        if not isinstance(v, self.oktypes):
            raise TypeError("{} not a valid type: {}".format(v, self.oktypes))

    def __len__(self):
        return len(self.list)

    def __getitem__(self, i):
        return self.list[i]

    def __delitem__(self, i):
        del self.list[i]

    def __setitem__(self, i, v):
        self.check(v)
        self.list[i] = v

    def insert(self, i, v):
        self.check(v)
        self.list.insert(i, v)

    def __str__(self):
        return str(self.list)

    def __repr__(self):
        return super(TypedList).__repr__() + " " + self.list.__repr__()

