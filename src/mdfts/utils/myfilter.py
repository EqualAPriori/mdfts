"""
Simple utility container for filtering, matching, etc.

Main requirement: that the objects stored have an `__eq__` function!

High level: 
- the Filter objects allow for ad infinitum nesting of objects, e.g. by nesting more Filters. This is probably overkill.
- the FilterSet object requires a more restricted interface of tuples of tuples.
- syntax is Filter(pattern1,pattern2,pattern3)
    which internally gets processed into a tuple of tuples or sets
- use itertools.product and itertools.permutations to generate all patterns to be matched

=====
Functionality:
- equals
- contains/subset/overlapping patterns
    essentially need to check all permutations of all products and check for overlap
- match(item)
- filter(list or generator) --> return only items that match the filter
- ability to filter on particular attributes/functions
- ordered/unordered

Nice to have:
- serializable

Things that make using set harder:
- if unhashable (e.g. using custom == behavior instead of checking object identity), then can't use set/dict
Things that make serializable harder:
- nested tuples, etc.

Sketch:
1) Approach 1: all elements are hashable and the user is sure that hash(x) == hash(y) is sufficient for satisfactory value equality
   challenge: if the filter definition itself uses nested lists/tuples that aren't themselves hashable...

    a) see if two filters are equal

        set(self) == set(other)

2) Approach 2: general unhashable case, O(n^2) implementation

    a) see if two filters are equal, relies on *value equality* (e.g. how list.index and `a in b` behaves)

        def sequences_contain_same_items(a, b):
            for item in a:
                try:
                    i = b.index(item)
                except ValueError:
                    return False
                b = b[:i] + b[i+1:]
            return not b

3) General functions:
    b) see if filter contains an item --> not really needed?

        other in self

    c) match

        if len(other) != len(self):
            return False
        for pattern in product(self):
            for ix,x in other:
                if x not in pattern[ix]: break
                return True # only if all x match to their corresponding patterns
        return False
    
    d) filter

        return [ x for x in other if self.match(x) ]

Desired usage:

f = Filter(obj1,obj2,[obj3,obj4])   #3-body filter. e.g. same syntax as itertools.
f.match(obj1)                       #returns False b/c not an equivalent # keys in each pattern
f.match(obj3,obj2,obj1)             #returns True

f = Filter([obj1,obj3])             #1-body filter with 2 elements
f.match(obj1)                       #returns True
                                    

f = Filter([obj1,obj3],obj1,obj2)   #1-body filter with 2 elements
f.match(obj1,obj3,obj2)             #returns True

Filter([obj1]) == Filter(obj1)      #returns True
                                    #DON'T try to overload syntax too much beyond this, or else get amibiguities
                                    #Another alternative is that the input must be single argument, e.g. f = Filter([obj1,obj2,[obj3,obj4]])
                                    #but can add that later

f1 = Filter([obj1,obj2],[obj3,obj4])
f2 = Filter([obj3,obj4],[obj2,obj1])
f1 == f2                            #returns True


"""
from __future__ import absolute_import, division, print_function
import collections
from multiprocessing.sharedctypes import Value

from pytest import Item

# from multiprocessing.sharedctypes import Value

try:
    collectionsABC = collections.abc
except:
    collectionsABC = collections


from collections import OrderedDict
import itertools
import warnings

import sys

sys.path.insert(1, "/Users/kshen/code/mdfts/src")
from mdfts.utils import serial


class A(serial.Serializable):
    _serial_vars = ["x"]

    def __init__(self, x):
        self.x = x

    def __eq__(self, other):
        if type(self) == type(other):
            return self.x == other.x
        else:
            return self.x == other

    @classmethod
    def init_from_dict(cls, d, *args, **kwargs):
        print(cls)
        try:
            obj = cls(0, *args, **kwargs)
        except:
            print("Could not initialize {}".format(cls))

        # final update
        obj.from_dict(d)
        return obj


@serial.serialize("x")
class A_hashable:
    def __init__(self, x):
        self.x = x

    def __eq__(self, other):
        if type(self) == type(other):
            return self.x == other.x
        else:
            return self.x == other

    def __hash__(self):
        # return hash((self.__class__, self.x)) #if want to distinguish from `x``
        return hash(self.x)  # if want object to look the same as `x`


class Filter(collectionsABC.Sequence):
    def __init__(self, *args, ordered=False):
        """better to use tuples instead of lists (since tuples are hashable)"""
        self._ordered = ordered
        self._pattern = tuple(
            [Filter.unique_values(pattern) for pattern in args]
        )  # args gets read in as a tuple

    # === SEQUENCE INTERFACE ===
    # def __iter__(self, *args):
    #    return self._pattern.__iter__()
    def __getitem__(self, ii):
        return self._pattern[ii]

    def __len__(self):
        return len(self._pattern)

    def __add__(self, other):
        if type(self) != type(other):
            raise ValueError("Can only add common Typed Filters together")
        else:
            return Filter(*(self._pattern + other._pattern))

    # === TYPEDFILTER INTERFACE ===
    def __eq__(self, other):
        if not isinstance(other, Filter):
            if isinstance(other, collectionsABC.Iterable):
                other = Filter(*other)
            else:
                other = Filter(other)

        if len(self) != len(other):
            return False

        if self._ordered:
            for ix, item1 in enumerate(self):
                item2 = other[ix]
                if not eq_order_insensitive(item1, item2):
                    return False
            return True

        for item1 in self:
            pattern_found = False
            for ind, item2 in enumerate(other):
                try:  # using set optimization
                    if set(item1) == set(item2):  # Filter._pattern should give tuples
                        other = other[:ind] + other[ind + 1 :]
                        pattern_found = True
                        break
                except TypeError:
                    if eq_order_insensitive(item1, item2):
                        other = other[:ind] + other[ind + 1 :]
                        pattern_found = True  # found a match, pop the element out and start searching for next item
                        break
            if not pattern_found:
                return False  # if made it through w/out finding a match for item1

        return not other  # if other now empty, then self==other

    def __repr__(self):
        # return "{}{}".format(self.__class__, self._pattern)
        return "Filter{}".format(self._pattern)

    # === REST of functionality ===
    def match(self, *args):
        """See if a set of values matches the filter's pattern
        Current syntax expects # arguments = # elements in the filter's pattern.
        This is to match the constructor.
        An alternative syntax is to take a *single* iterable containing the patterns.

        Returns:
            bool: whether or not match was found
        """
        if len(args) != len(self):
            return False
        # iterate over permutations of the filter's pattern until we get one
        # that matches the order the arguments are input in
        if self._ordered:
            permutations = [self._pattern]
        else:
            permutations = itertools.permutations(self._pattern)

        for pattern in permutations:
            print("pattern {}".format(pattern))
            for ix, x in enumerate(args):
                if isinstance(pattern[ix], collectionsABC.Iterable):
                    if x not in pattern[ix]:
                        break  # try next pattern
                else:
                    if x != pattern[ix]:
                        break  # try next pattern
                return True  # only if all x match to their corresponding patterns
        return False

    def __contains__(self, pattern):
        return self.match(*pattern)

    def check_overlap(self, other):
        """Returns common matching subpatterns.
        The trick is to check all permutations of the patterns for completeness,
        but when outputting to not put repeat patterns.

        Args:
            other (Filter-like): e.g. filter._pattern also works

        Returns:
            list: common matching subpatterns
        """
        if len(other) != len(self):
            return False
        matches = []

        for subpattern in itertools.product(*self._pattern):
            # print("SUBPATTERN {}".format(subpattern))
            if self._ordered == True:
                if isinstance(other, Filter):
                    permutation_iterator = [other._pattern]
                else:
                    permutation_iterator = [other]
            elif isinstance(other, Filter):
                permutation_iterator = itertools.permutations(other._pattern)
            else:
                permutation_iterator = itertools.permutations(other)
            for permutation in permutation_iterator:
                for other_subpattern in itertools.product(*permutation):
                    # print("other {}".format(other_subpattern))
                    if subpattern == other_subpattern:
                        if subpattern not in matches:
                            matches.append(subpattern)
        return matches

        # i.e. check if two patterns might have overlapping "members"
        # maybe can check Cartesian product

    def filter(self, iterator):
        for el in iterator:
            if self.match(el):
                yield el

    @staticmethod
    def unique_values(x):
        """Takes an iterable or an instance and returns unique *values*, as determined by __eq__

        Args:
            x (any): thing to extract unique values from

        Returns:
            tuple or instance: returns instance/object
        """
        unique = []
        if isinstance(x, collectionsABC.Iterable) and not isinstance(x, str):
            for el in x:
                if isinstance(el, collectionsABC.Iterable) and not isinstance(
                    el, (Filter, str)
                ):
                    el1 = Filter(el)
                else:
                    el1 = el
                # print(el1, unique)
                if el1 not in unique:
                    unique.append(el1)
            unique = tuple(unique)  # but this still preserves order...
            return unique
        else:
            return (x,)  # `x` prints nicer, but (x,) is iter-ready for itertools


class TypedFilter(collectionsABC.Sequence, serial.Serializable):
    """(Serializable) Typed Filter
    Potentially can use hash optimizations, if the `oktype` is hashable.
    Should be immutable in usage, also good for hashing.

    Simplified implementation: don't need ad-infinitum nesting of Filters.
    Simply as a tuple of tuples (of objects) is enough.
    """

    def __init__(self, *args, oktype=None, ordered=False):
        self._oktype = oktype
        self._ordered = ordered
        self._pattern = tuple(
            [self.unique_values(pattern) for pattern in args]
        )  # args gets read in as a tuple
        self._serial_vars = []  # not used since self._patterns stores the contents
        self._extra_vars = None

    # === SEQUENCE INTERFACE ===
    def __getitem__(self, ii):
        # leads to infinite recursion because __eq__ tries to cast contents into Filter,
        # which triggers another equality check and can never unpack contents
        # return TypedFilter(*self._pattern[ii], oktype=self._oktype)
        return self._pattern[ii]

    def __len__(self):
        return len(self._pattern)

    # === TypedFilter INTERFACE ===
    # TODO: build in ordered checking
    def __eq__(self, other):
        if not isinstance(other, (TypedFilter, Filter)):
            # Then turn into a filter, which pre-processes uniqueness
            if isinstance(other, collectionsABC.Iterable):
                other = TypedFilter(*other, oktype=self._oktype)
            else:  # i.e. a "1-body" filter --> provides equality checking against tuples
                other = TypedFilter(other, oktype=self._oktype)
        else:
            if self._oktype != other._oktype:
                return False

        if len(self) != len(other):
            return False

        if self._ordered:
            for ix, item1 in enumerate(self):
                item2 = other[ix]
                if item1 != item2:
                    return False
            return True

        for item1 in self:
            subpattern_found = False
            for ind, item2 in enumerate(other):
                # print(item1, item2)
                if eq_order_insensitive(item1, item2):
                    other = other[:ind] + other[ind + 1 :]
                    subpattern_found = True  # found a match, pop the element out and start searching for next item
            if not subpattern_found:
                return False  # if made it through w/out finding a match for item1

        return not other  # if other now empty,f1 then self==other

    def __add__(self, other):
        if type(self) != type(other):
            raise ValueError("Can only add common Typed Filters together")
        elif self._oktype != other._oktype:
            warnings.warn("Dissimilar Filter Types, using untyped Filter instead")
            return Filter(*(self._pattern + other._pattern))
        else:
            return TypedFilter(*(self._pattern + other._pattern), oktype=self._oktype)

    def __repr__(self):
        # return "{}{}".format(self.__class__, self._pattern)
        return "Filter{}{}".format(self._oktype, self._pattern)

    # === SERIALIZABLE INTERFACE ===

    # === Other Functionality ===
    def __contains__(self, pattern):
        return self.match(*pattern)

    def match(self, *args):
        """See if a set of values matches the filter's pattern
        Current syntax expects # arguments = # elements in the filter's pattern.
        This is to match the constructor.
        An alternative syntax is to take a *single* iterable containing the patterns.

        Returns:
            [type]: [description]
        """
        if len(args) != len(self):
            return False
        # iterate over permutations of the filter's pattern until we get one
        # that matches the order the arguments are input in
        if self._ordered:
            permutations = [self._pattern]
        else:
            permutations = itertools.permutations(self._pattern)

        for pattern in permutations:
            print("pattern {}".format(pattern))
            for ix, x in enumerate(args):
                if isinstance(pattern[ix], collectionsABC.Iterable):
                    if x not in pattern[ix]:
                        break  # try next pattern
                else:
                    if x != pattern[ix]:
                        break  # try next pattern
                return True  # only if all x match to their corresponding patterns
        return False

    def check_overlap(self, other):
        """Returns common matching subpatterns.
        The trick is to check all permutations of the patterns for completeness,
        but when outputting to not put repeat patterns.

        Args:
            other (Filter-like): e.g. filter._pattern also works

        Returns:
            list: common matching subpatterns
        """
        if len(other) != len(self):
            return False
        matches = []

        for subpattern in itertools.product(*self._pattern):
            # print("SUBPATTERN {}".format(subpattern))
            if self._ordered == True:
                if isinstance(other, Filter):
                    permutation_iterator = [other._pattern]
                else:
                    permutation_iterator = [other]
            elif isinstance(other, Filter):
                permutation_iterator = itertools.permutations(other._pattern)
            else:
                permutation_iterator = itertools.permutations(other)
            for permutation in permutation_iterator:
                for other_subpattern in itertools.product(*permutation):
                    # print("other {}".format(other_subpattern))
                    if subpattern == other_subpattern:
                        if subpattern not in matches:
                            matches.append(subpattern)
        return matches

        # i.e. check if two patterns might have overlapping "members"
        # maybe can check Cartesian product

    def filter(self, iterator):
        for el in iterator:
            if self.match(el):
                yield el

    def unique_values(self, x):
        """Takes an iterable or an instance and returns unique *values*, as determined by __eq__

        Args:
            x (any): thing to extract unique values from

        Returns:
            tuple or instance: returns instance/object
        """
        unique = []
        if isinstance(x, collectionsABC.Iterable) and not isinstance(x, str):
            for el in x:
                print(el, x)
                if isinstance(el, collectionsABC.Iterable) and not isinstance(
                    el, (TypedFilter, str)
                ):
                    el1 = TypedFilter(self._oktype, el)
                elif isinstance(el, Filter):
                    el1 = TypedFilter(self._oktype, *el._pattern)
                else:
                    el1 = el
                if not isinstance(el1, self._oktype) or (
                    isinstance(el1, TypedFilter) and el1._oktype != self._oktype
                ):
                    raise TypeError(
                        "This typed filter requires type {}, received {}".format(
                            self._oktype, el1
                        )
                    )
                if el1 not in unique:
                    unique.append(el1)
            unique = tuple(unique)  # but this still preserves order...
            return unique
        else:
            return (x,)  # `x` prints nicer, but (x,) is iter-ready for itertools


class FilterSet0(collectionsABC.Sequence):
    """A FilterSet that only works with hashable (and ideally immutable) types.

    Essentially, processes arguments into a tuple of (frozen) sets.

    Only nests 1 deep.
    """

    def __init__(self, *args, ordered=False):
        self._pattern = self.process_args(args)
        self._ordered = ordered
        self._patterns = list(itertools.product(*self._pattern))
        self._serial_vars = [
            "_pattern",
            "_ordered",
        ]
        self._extra_vars = None

    def __eq__(self, other):
        if self._ordered:
            return self._pattern == other._pattern
        else:
            for p in itertools.permutations(other._pattern):
                if p == self._pattern:
                    return True
            return False

    def __hash__(self):
        return hash((self._pattern, self.ordered))

    def __repr__(self):
        return "FilterSet {}".format(self._pattern)

    def __getitem__(self, ii):
        return self._pattern[ii]

    def __len__(self):
        return len(self._pattern)

    def __add__(self, other):
        if type(self) != type(other):
            raise ValueError("Can only add FilterSets together")
        else:
            return FilterSet(*(self._pattern + other._pattern), oktype=self._oktype)

    def __contains__(self, pattern):
        return self.match(*pattern)

    def match(self, *args, mode=0):
        """See if a set of values matches the filter's pattern
        Current syntax expects # arguments = # elements in the filter's pattern.
        This is to match the constructor.
        An alternative syntax is to take a *single* iterable containing the patterns.

        Returns:
            bool: whether or not match is found
        """
        return match(self._pattern, args, self._ordered, mode=0)

    def check_overlap(self, other):
        """Returns common matching subpatterns.
        The trick is to check all permutations of the patterns for completeness,
        but when outputting to not put repeat patterns.

        Args:
            other (Filter-like): e.g. filter._pattern also works

        Returns:
            list: common matching subpatterns
        """
        if len(other) != len(self):
            return False
        matches = []

        for subpattern in itertools.product(*self._pattern):
            # print("SUBPATTERN {}".format(subpattern))
            if self._ordered == True:
                if isinstance(other, Filter):
                    permutation_iterator = [other._pattern]
                else:
                    permutation_iterator = [other]
            elif isinstance(other, Filter):
                permutation_iterator = itertools.permutations(other._pattern)
            else:
                permutation_iterator = itertools.permutations(other)
            for permutation in permutation_iterator:
                for other_subpattern in itertools.product(*permutation):
                    # print("other {}".format(other_subpattern))
                    if subpattern == other_subpattern:
                        if subpattern not in matches:
                            matches.append(subpattern)
        return matches

    def process_args(self, args):
        # Only nest 1-deep
        if not isinstance(args, collectionsABC.Iterable):  # single element
            return (frozenset((args,)),)
        else:  # This is the typical case, as *args gets processed into at least (arg1,)
            pattern = []
            for arg in args:
                if isinstance(arg, collectionsABC.Iterable):
                    pattern.append(frozenset(arg))
                else:
                    pattern.append(frozenset([arg]))
            return tuple(pattern)


class FilterSet(collectionsABC.Sequence):
    """A TypedFilterSet that only preferably with hashable (and ideally immutable) types.

    Essentially, processes arguments into a tuple of (frozen) sets.

    Only nests 1 deep.
    """

    def __init__(self, *args, ordered=False):
        self._pattern, self._hashable = process_pattern(args)
        self._ordered = ordered
        self._oktype = self.infer_type()

        self._patterns = list(itertools.product(*self._pattern))
        self._serial_vars = [
            "_pattern",
            "_ordered",
        ]
        self._extra_vars = None

    # === SEQUENCE INTERFACE AND DUNDER METHODS ===
    def __eq__(self, other):
        if not isinstance(other, FilterSet):
            if not isinstance(other, collectionsABC.Iterable):
                other = (other,)
            else:
                other = FilterSet(*other, ordered=self._ordered)
        if len(self) != len(other) or self._ordered != other._ordered:
            print("fail1")
            return False

        if self._hashable:  # use set equality for order-insensitive subpattern matching
            if self._ordered:
                return self._pattern == other._pattern
            else:
                for p in itertools.permutations(other._pattern):
                    if p == self._pattern:
                        return True
                print("fail2")
                return False
        else:  # need custom function for order-insensitive matching of subpatterns
            if self._ordered:
                for item1, item2 in zip(self, other):
                    if not eq_order_insensitive(item1, item2):
                        print("fail3")
                        return False
                return True
            else:
                # custom handling of order-insensitive tuple of order-insensitive tuples
                return eq_order_insensitive(self, other)

    def __hash__(self):
        return hash((self._pattern, self.ordered))

    def __repr__(self):
        return "FilterSet {}".format(self._pattern)

    def __getitem__(self, ii):
        return self._pattern[ii]

    def __len__(self):
        return len(self._pattern)

    def __add__(self, other):
        if type(self) != type(other):
            raise ValueError("Can only add FilterSets together")
        else:
            return FilterSet(*(self._pattern + other._pattern), oktype=self._oktype)
        """
        elif self._oktype != other._oktype:
            warnings.warn("Dissimilar Filter Types, using untyped Filter instead")
            return Filter(*(self._pattern + other._pattern))
        else:
            return TypedSet(*(self._pattern + other._pattern), oktype=self._oktype)
        """

    def __contains__(self, pattern):
        return self.match(*pattern)

    # === CUSTOM Filter METHODS ===
    def filter(self, iterator, pre_process=None, post_process=None):
        for el in iterator:
            if pre_process is not None:
                target = pre_process(el)  # e.g. extract an attribute only
            else:
                target = el

            if self.match(target):
                if post_process(el):  # e.g. additional checks, like bonding
                    yield el

    def match(self, *args, mode=0):
        """See if a set of values matches the filter's pattern

        Current syntax matches constructor syntax, e.g. like numpy.random.rand(2,2).
        Alternative syntax is to take a *single* iterable, e.g. like numpy.random.random([2,2])

        Returns:
            bool: whether or not match is found
        """
        return match(self._pattern, args, self._ordered, mode=0)

    def check_overlap(self, other):
        """Returns common matching subpatterns.
        The trick is to check all permutations of the patterns for completeness,
        but when outputting to not put repeat patterns.

        Args:
            other (Filter-like): e.g. filter._pattern also works

        Returns:
            list: common matching subpatterns
        """
        if len(other) != len(self):
            return False
        matches = []

        for subpattern in itertools.product(*self._pattern):
            # print("SUBPATTERN {}".format(subpattern))
            if self._ordered == True:
                if isinstance(other, Filter):
                    permutation_iterator = [other._pattern]
                else:
                    permutation_iterator = [other]
            elif isinstance(other, Filter):
                permutation_iterator = itertools.permutations(other._pattern)
            else:
                permutation_iterator = itertools.permutations(other)
            for permutation in permutation_iterator:
                for other_subpattern in itertools.product(*permutation):
                    # print("other {}".format(other_subpattern))
                    if subpattern == other_subpattern:
                        if subpattern not in matches:
                            matches.append(subpattern)
        return matches

    def infer_type(self):
        seen_types = []
        if len(self._pattern) == 0:
            return None
        for sub_pattern in self._pattern:
            # self._pattern should already be a tuple of tuples
            for el in sub_pattern:
                if type(el) not in seen_types:
                    seen_types.append(type(el))
        if all([t == seen_types[0] for t in seen_types]):
            return seen_types[0]
        else:
            return None


class SerializableFilterSet(FilterSet, serial.Serializable):
    def __init__(self, *args, ordered=False):
        super(SerializableFilterSet, self).__init__(*args, ordered=ordered)
        if self._oktype is None:
            if len(self._pattern) > 0:
                raise TypeError(
                    "SerializableFilter must have everything be a single type!"
                )
            else:
                warnings.warn("Empty filter, still need to specify a type")
        _serial_vars = [
            "_pattern",
            "_ordered",
            "_oktype",
        ]

    # === SERIALIZABLE INTERFACE ===
    def custom_get(self, k):
        """_summary_

        Args:
            k (str): dictionary key

        Raises:
            KeyError: if not a standard key of SerializableFilterSets

        Returns:
            depends: on the field requested

        the _pattern part likely needs over-riding for more complex objects
        """
        if k in ["_ordered", "ordered"]:
            return self._ordered
        elif k in ["_oktype", "type", "oktype"]:
            return self._oktype
        elif k in ["_pattern", "pattern"]:
            res = []
            for sub_pattern in self._pattern:
                subres = []
                for sp in sub_pattern:
                    if isinstance(sp, serial.Serializable):
                        subres.append(sp.to_dict())
                    else:
                        subres.append(sp)
                res.append(subres)
            return res
        else:
            raise KeyError("Unknown key type for a SerializableFilterSet")

    def custom_set(self, k, v, decoder=None):
        """_summary_

        Args:
            k (str): key
            v (depends): value to set. should be in *decoded* format.

        Note:
            do the type coercing here.
            Maybe... need another code in the schema who knows the correct type to coerce things into!
        """
        if k in ["_ordered", "ordered"]:
            if type(v) == bool:
                self._ordered = v
            else:
                raise ValueError("value for key `ordered` must be `bool`")
        elif k in ["_oktype", "type", "oktype"]:
            self._oktype = oktype
            # really should not allow setting of oktype this way, dangerous
        elif k in ["_pattern", "pattern"]:
            res = []
            pattern, _ = process_pattern(v)
            print(pattern)
            for sub_pattern in pattern:
                # pattern should be list of lists!
                subres = []
                for sp in sub_pattern:
                    if issubclass(self._oktype, serial.Serializable):
                        subres.append(self._oktype.init_from_dict(sp))
                    else:
                        subres.append(serial.decode(sp, self._oktype))
                res.append(tuple(subres))
            self._pattern, self._hashable = process_pattern(res)
        else:
            raise KeyError("Unknown key type for a SerializableFilterSet")

        # self._pattern, self._hashable = process_pattern(args)
        # self._oktype = self.infer_type()

    # @classmethod
    # def init_from_dict(cls, d, *args, **kwargs):
    #    pass


# === Utility Functions Performing Some Set-like Operations and Pattern Expansion, for Un-hashable objects
def match(master_pattern, target_pattern, ordered=False, mode=0):
    """Matching patterns represented by tuples

    Args:
        master_pattern (tuple): each element can itself be an iterable with several options that need to be searched over
        target_pattern (tuple): should be non-iterable in each element
        ordered (bool, optional): whether or not the match should be ordered. Defaults to False.
        mode (int, optional): which algorithm for doing matching. Defaults to 0.

    Returns:
        bool: whether target_pattern was found in master_pattern

    Note:
        patterns should be n-tuples where each element is itself an iterable, even for 1-element,
        e.g. ((x,),(1,2)) is a 2-body filter with 'x' in the first slot and 1 or 2 in the second slot.
        The code will also handle (x,(1,2)) for convenience, but the above should be preferred for clarity.
    """
    if len(target_pattern) != len(master_pattern):
        return False
    # iterate over permutations of the filter's pattern until we get one
    # that matches the order the arguments are input in

    if mode == 0:  # naive way
        if ordered:
            permutations = [master_pattern]
        else:
            permutations = itertools.permutations(master_pattern)

        for pattern in permutations:
            matched = []
            for ix, x in enumerate(target_pattern):
                if isinstance(pattern[ix], collectionsABC.Iterable):
                    if isinstance(x, collectionsABC.Iterable) and len(list(x)) == 1:
                        if x[0] not in pattern[ix]:
                            matched.append(False)
                            break
                    elif x not in pattern[ix]:
                        matched.append(False)
                        break  # try next pattern
                else:  # in case the pattern is given as a single element instead of the less ambiguous (x,) format
                    if x != pattern[ix]:
                        matched.append(False)
                        break  # try next pattern

                matched.append(True)

            if all(matched):
                return True  # only if all x match to their corresponding patterns
    elif mode == 1:  # itertools way, only works if everything is given in (1,) format
        # costly, preferably patterns are pre-computed
        patterns = list(itertools.product(master_pattern))
        if ordered:
            return target_pattern in patterns
        else:
            for p in itertools.permutations(target_pattern):
                if p in patterns:  # try next pattern
                    return True

    return False


def uniqify(x):
    """Takes an iterable or an instance and returns unique *values*, as determined by __eq__
    this implementation does *not* nest and recurse through.

    This is for un-hashables that don't use set.

    Args:
        x (iterable): thing to extract unique values from

    Returns:
        tuple: returns tuple of processed values
    """
    if isinstance(x, collectionsABC.Iterable) and not isinstance(x, str):
        unique = []
        for el in x:
            if el not in unique:
                unique.append(el)
        return tuple(unique)
        # or, for fun, a 2-liner with list comprehension:
        # unique = []
        # [x for x in l if x not in used and (unique.append(x) or True)]
    else:
        return (x,)  # `x` prints nicer, but (x,) is iter-ready for itertools


def process_pattern(args):
    """Takes a list of arguments and turns it into a pattern that is a tuple of tuples (or frozen sets)

    Consider alternative syntax: take in unpacked sequence of arguments (i.e. `def process_pattern(*args)`)

    Returns:
        tuple: of tuples/frozen sets
        bool: whether everything was hashable
    """
    # Only nest 1-deep!
    if isinstance(args, str) or not isinstance(
        args, collectionsABC.Iterable
    ):  # single element
        if isinstance(args, collectionsABC.Hashable):
            hashable = True
            pattern = (frozenset((args,)),)
        else:
            hashable = False
            pattern = ((args,),)
    else:  # This is the typical case, as *args gets processed into at least (arg1,)
        preprocessed = [
            tuple(x) if isinstance(x, collectionsABC.Iterable) else x for x in args
        ]
        print(preprocessed)

        hashable = all([isinstance(x, collectionsABC.Hashable) for x in preprocessed])

        pattern = []
        if hashable:
            try:
                for arg in preprocessed:
                    if isinstance(arg, collectionsABC.Iterable):
                        pattern.append(frozenset(arg))
                    else:
                        pattern.append(frozenset([arg]))
            except:
                hashable = False
        if not hashable:
            for arg in preprocessed:
                if isinstance(arg, collectionsABC.Iterable):
                    pattern.append(uniqify(arg))
                else:
                    pattern.append(uniqify(arg))

    return tuple(pattern), hashable


def eq_order_insensitive(x, y):
    """Checks two iterables to see if they have the same contents

    In recursively order-insensitive form.

    If all elements are unique (i.e. no repeats), this is like testing set equality by value
    that works for unhasable objects that have cutom __eq__().

    The main difference is that 1-element iterables are treated as if they were just their content,
    i.e. 1 == (1,) by this function.

    Args:
        x (tuple, non-iterable):
        y (tuple, non-iterable):
    """
    if isinstance(x, collectionsABC.Iterable) and isinstance(
        y, collectionsABC.Iterable
    ):
        if (
            isinstance(x, collectionsABC.Collection)
            and isinstance(y, collectionsABC.Collection)
            and len(x) != len(y)
        ):
            return False
        for item1 in x:
            pattern_found = False
            for ind, item2 in enumerate(y):
                if eq_order_insensitive(item1, item2):
                    y = y[:ind] + y[ind + 1 :]
                    pattern_found = True
                    break  # found match, popped it, search for next item1
            if not pattern_found:  # short circuit optimization; not necessary
                break  # this is only hit if no matches are found and the loop didn't break earlier
            """ old implementation that didn't do recursion:
            try:
                i = y.index(item1)
            except ValueError:
                return False
            y = y[:i] + y[i + 1 :]
            """
        return not y
    else:
        if isinstance(x, collections.Iterable) and len(list(x)) == 1:
            x = x[0]
        if isinstance(y, collections.Iterable) and len(list(y)) == 1:
            y = y[0]
        return x == y


if __name__ == "__main__":
    f1 = Filter(A(1), A(2))
    f2 = Filter(1, 2)
    print(f1 == f2)  # True

    f3 = TypedFilter(A(1), A(2), oktype=A)
    f4 = TypedFilter(1, 2, oktype=int)
    f5 = TypedFilter(A(2), A(1), oktype=A)
    print(f3 == f5)  # True

    fs1 = FilterSet(1, 2, 3, 4, (4, 4, 5))
    fs2 = FilterSet(A_hashable(1), 2, 3, 4, (4, 4, 5))
    print(fs1 == fs2)  # True

    import time

    start = time.time()
    for ii in range(1000):
        fs1.match(1, 5, 3, 4, 2)
    print("took {}s to do 1000x matches".format(time.time() - start))

    sfs = SerializableFilterSet(A(1), A(2))
    d = sfs.to_dict()
    print(d)
    sfs.from_dict(d)
    sfs.match(1, 2)  # True
    sfs.match((1,), (2,))  # True
    sfs.match(1, 3)  # False

    # to make a new SerializableFilterSet:
    sfs1 = SerializableFilterSet()
    sfs1._oktype = int
    sfs1.custom_set("pattern", [1])
    sfs1.custom_set("_pattern", [[1]])  # equivalent
    sfs1.custom_set("_pattern", [(1,)])  # equivalent

    sfs2 = SerializableFilterSet()
    sfs2._oktype = A
    sfs2.from_dict(d)

    sfs == sfs2  # True
    sfs2 == sfs  # True
