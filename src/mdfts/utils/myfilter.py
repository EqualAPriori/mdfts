"""
Simple utility container for filtering, matching, etc.

Main requirement: that the objects stored have an `__eq__` function!

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


class A:
    def __init__(self, x):
        self.x = x

    def __eq__(self, other):
        if type(self) == type(other):
            return self.x == other.x
        else:
            return self.x == other


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
                if item1 != item2:
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
                    if _eq_unordered_seq(item1, item2):
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
        self.oktype = oktype
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
        # return TypedFilter(*self._pattern[ii], oktype=self.oktype)
        return self._pattern[ii]

    def __len__(self):
        return len(self._pattern)

    # === TypedFilter INTERFACE ===
    # TODO: build in ordered checking
    def __eq__(self, other):
        if not isinstance(other, (TypedFilter, Filter)):
            # Then turn into a filter, which pre-processes uniqueness
            if isinstance(other, collectionsABC.Iterable):
                other = TypedFilter(*other, oktype=self.oktype)
            else:  # i.e. a "1-body" filter --> provides equality checking against tuples
                other = TypedFilter(other, oktype=self.oktype)
        else:
            if self.oktype != other.oktype:
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
                if _eq_unordered_seq(item1, item2):
                    other = other[:ind] + other[ind + 1 :]
                    subpattern_found = True  # found a match, pop the element out and start searching for next item
            if not subpattern_found:
                return False  # if made it through w/out finding a match for item1

        return not other  # if other now empty,f1 then self==other

    def __add__(self, other):
        if type(self) != type(other):
            raise ValueError("Can only add common Typed Filters together")
        elif self.oktype != other.oktype:
            warnings.warn("Dissimilar Filter Types, using untyped Filter instead")
            return Filter(*(self._pattern + other._pattern))
        else:
            return TypedFilter(*(self._pattern + other._pattern), oktype=self.oktype)

    def __repr__(self):
        # return "{}{}".format(self.__class__, self._pattern)
        return "Filter{}{}".format(self.oktype, self._pattern)

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
                    el1 = TypedFilter(self.oktype, el)
                elif isinstance(el, Filter):
                    el1 = TypedFilter(self.oktype, *el._pattern)
                else:
                    el1 = el
                if not isinstance(el1, self.oktype) or (
                    isinstance(el1, TypedFilter) and el1.oktype != self.oktype
                ):
                    raise TypeError(
                        "This typed filter requires type {}, received {}".format(
                            self.oktype, el1
                        )
                    )
                if el1 not in unique:
                    unique.append(el1)
            unique = tuple(unique)  # but this still preserves order...
            return unique
        else:
            return (x,)  # `x` prints nicer, but (x,) is iter-ready for itertools


def _eq_unordered_seq(x, y):
    """for pattern-by-pattern comparison in filter

    Args:
        x (tuple, non-iterable):
        y (tuple, non-iterable):
    """
    if isinstance(x, collectionsABC.Iterable) and isinstance(
        y, collectionsABC.Iterable
    ):
        # print(x, y)
        for item in x:
            try:
                i = y.index(item)
            except ValueError:
                return False
            y = y[:i] + y[i + 1 :]
        return not y
    else:
        if isinstance(x, collections.Iterable) and len(list(x)) == 1:
            x = x[0]
        if isinstance(y, collections.Iterable) and len(list(y)) == 1:
            y = y[0]
        return x == y


class FilterSet:
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
        return hash(self._pattern)

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
            return FilterSet(*(self._pattern + other._pattern), oktype=self.oktype)
        """
        elif self.oktype != other.oktype:
            warnings.warn("Dissimilar Filter Types, using untyped Filter instead")
            return Filter(*(self._pattern + other._pattern))
        else:
            return TypedSet(*(self._pattern + other._pattern), oktype=self.oktype)
        """

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

        """ alternative implementation, simpler!
        for p in itertools.permutations(args):
            if p in self._patterns:  # try next pattern
                return True
        """

        for pattern in permutations:
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


'''
from mdfts.utils import serial

class Filter(serial.Serializable):
    """1-element filter
    Require that the things to be filtered are hashable!!!    
    """

    _serial_vars = ["components"]

    def __init__(self, *args):

        self.components = args  # tuple of things
        n_els = len(args)
        # verify that everything is eq-able
        for ind1, el1 in enumerate(args):
            try:
                el1 == el1
            except NotImplementedError:
                raise ValueError(
                    "object {} ({}) does not have __eq__ functionality".format(
                        ind1, el1
                    )
                )

            for ind2 in range(ind1 + 1, len(args)):
                el2 = args[ind2]
                try:
                    el1 == el2
                except NotImplementedError:
                    raise ValueError(
                        "obj {} ({}) and {} ({}) can not be __eq__ compared".format(
                            ind1, el1, ind2, el2
                        )
                    )

    def __eq__(self, other):
        # order shouldn't matter
        pass


class FilterN(serial.Serializable):
    """multi-element filter, simply provides shorthand for using Filter."""

    _serial_vars = serial.SerializableTypedList(Filter)
'''
if __name__ == "__main__":
    f1 = Filter(A(1), A(2))
    f2 = Filter(1, 2)
    print(f1 == f2)  # True

    f3 = TypedFilter(A(1), A(2), oktype=A)
    f4 = TypedFilter(1, 2, oktype=int)
    f5 = TypedFilter(A(2), A(1), oktype=A)
    f3 == f5  # True

    fs1 = FilterSet(1, 2, 3, 4, (4, 4, 5))
    fs2 = FilterSet(A_hashable(1), 2, 3, 4, (4, 4, 5))
    fs1 == fs2  # True

    import time

    start = time.time()
    for ii in range(1000):
        fs1.match(1, 5, 3, 4, 2)
    print("took {}s".format(time.time() - start))
