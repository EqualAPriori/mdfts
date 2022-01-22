""" Helper tools to create objects that automatically serialize and deserialize

By serializing to dictionary, can trivially read/write with (e.g. json).
(Provided that the dictionary only has basic data types, we don't do type checking here)

Note:
    I envision simply needing to have classes subclass Serializable...
    or use decorators to track which variables need serialization/deseriazliation

Todo:
    - decorators to track variables that can be saved to and read from dictionaries.
    - different use-methods:
        - subclass Serializable, and define _serial_vars and _extra_vars accordingly
        - use metaclass as aid to ensure _serial_vars is indeed defined
        - decorator that takes _serial_vars as argument
"""
from __future__ import absolute_import, division, print_function

from collections import OrderedDict
import collections

try:
    collectionsABC = collections.abc
except:
    collectionsABC = collections

import inspect


# def serialize(tracked_args):


class serial_tracking(type):
    """metaclass to inject default class variables into a class

    Can also use to change the __str__ and __repr__ of a class.
    AH, e.g. I can use the metaclass in a decorator!
    """

    def __new__(cls, clsname, bases, clsdict):
        c = super(serial_tracking).__new__(cls, clsname, bases, clsdict)
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


class SerializableTypedList(collections.MutableSequence, Serializable):
    """Create custom list that type-checks its contents.
    
    TypedList, restricted for use with Serializable classes. 
    For creating a quick list containers when creating a new class is overkill.
    
    Note:
        taken from https://stackoverflow.com/questions/3487434/overriding-append-method-after-inheriting-from-a-python-list
        
        This is just a reference.
    """

    def __init__(self, oktype, *args):
        if not issubclass(oktype, Serializable):
            raise TypeError("list contents must be a Serializable Class")
        self.oktype = oktype
        self.list = list()
        self.extend(list(args))
        self._serial_vars = []  # not used since we extend list and override to_dict()
        self._extra_vars = None

    # === SERIALIZABLE INTERFACE ===
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
        return super(SerializableTypedList).__repr__() + " " + self.list.__repr__()
