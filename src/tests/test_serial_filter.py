from __future__ import absolute_import, division, print_function

import unittest
from collections import OrderedDict

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


@serial.serialize(["x", "y"])
class B:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y


class TestFilter(unittest.TestCase):
    def setUp(self):
        pass

    def run_filters(self, filterclass, obj1, obj2, obj3, obj4):
        """basic tests to be repeated for different object types"""
        # 3-body filter. e.g. same syntax as itertools.
        f = filterclass(obj1, obj2, [obj3, obj4])

        self.assertFalse(f.match(obj1))  # returns False b/c unequal # keys in pattern
        self.assertTrue(f.match(obj3, obj2, obj1))  # order-insensitive by default
        self.assertTrue((obj3, obj2, obj1) in f)  # alternative syntax for match

        f = filterclass([obj1, obj3])  # 1-body filter with 2 elements
        self.assertTrue(f.match(obj1))

        f = filterclass([obj1, obj3], obj1, obj2)  # 1-body filter with 2 elements
        self.assertTrue(f.match(obj1, obj3, obj2))

        self.assertTrue(f.match(obj2, (obj3,), obj1))
        f.match(obj2, (obj3,), obj2)
        self.assertFalse(f.match(obj2, (obj3,), obj2))

    def test_FilterInt(self):
        self.run_filters(serial.FilterSet, 1, 2, 3, 4)

    def test_FilterStr(self):
        self.run_filters(serial.FilterSet, "a", "b", "cd", "qwer")

    def test_FilterObj(self):
        self.run_filters(serial.FilterSet, B(1, 2), B(2, 3), B(4, 5), B(5, 6))

    # Next need to test a serializable (and typed) version of the filter

    # test FilterSet
    def test_FilterSet(self):
        self.run_filters(serial.FilterSet, 1, 2, 3, 4)

        self.assertTrue(serial.eq_order_insensitive(((2, 1), (1, 2)), ((2, 1), (1, 2))))
        self.assertTrue(
            serial.eq_order_insensitive(((2, 1, 3), (1, 2)), ((2, 1), (1, 2, 3)))
        )

    def test_TypedFilterSet(self):
        self.run_filters(serial.SerializableFilterSet, 1, 2, 3, 4)

        sfs0 = serial.SerializableFilterSet(A(1), A(2))
        d = sfs0.to_dict()
        sfs0.from_dict(d)
        self.assertTrue(sfs0.match(1, 2))
        self.assertTrue(sfs0.match((1,), (2,)))
        self.assertFalse(sfs0.match(1, 3))

        # to make a new SerializableFilterSet:
        sfs1 = serial.SerializableFilterSet()
        sfs1._oktype = int
        sfs2 = serial.SerializableFilterSet()
        sfs2._oktype = int
        sfs1.custom_set("pattern", [1])
        sfs2.custom_set("_pattern", [[1]])  # equivalent
        self.assertTrue(sfs1 == sfs2)
        sfs1.custom_set("_pattern", [(1,)])  # equivalent
        self.assertTrue(sfs1 == sfs2)

        sfs3 = serial.SerializableFilterSet()
        sfs3._oktype = A
        sfs3.from_dict(d)
        self.assertTrue(sfs0 == sfs3)
        self.assertTrue(sfs3 == sfs0)


if __name__ == "__main__":
    unittest.main()
