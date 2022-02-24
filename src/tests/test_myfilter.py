from __future__ import absolute_import, division, print_function

import unittest
from collections import OrderedDict

from mdfts.utils import myfilter, serial


@serial.serialize(["x", "y"])
class A:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y


class TestFilter(unittest.TestCase):
    def setUp(self):
        pass

    def run_1el_filters(self, o1, o2):
        # make filters
        f1 = myfilter.Filter(o1)
        f2 = myfilter.Filter(o1)
        f3 = myfilter.Filter(o2)
        f4 = myfilter.Filter(None)

        # check if filters equal
        self.assertEqual(f1, f1)
        self.assertEqual(f1, f2)
        self.assertFals(f1, f3)

        # check if tuple input works
        self.assertTrue(f1.equivalent(f2))
        self.assertTrue(f1.equivalent(f2))
        self.assertFalse(f1.equivalent(f3))
        self.assertTrue(f1.equivalent(o1))
        self.assertTrue(f1.equivalent([o1]))
        # check

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
        print("THIS TEST")
        f.match(obj2, (obj3,), obj2)
        self.assertFalse(f.match(obj2, (obj3,), obj2))

    def test_FilterBasic(self):
        # Test basic syntax and processing
        f1 = myfilter.Filter((1, 2, (3, 4, 5, (6, 7))), 2, (3, 3, 4), 6)
        expected_output = "((1, 2, Filter((3, 4, 5, Filter((6, 7),)),)), 2, (3, 4), 6)"
        self.assertTrue(f1._pattern.__str__(), expected_output)

        # Test that the output allows for round-trip initialization of a filter
        f2 = myfilter.Filter(
            (
                1,
                2,
                myfilter.Filter(
                    (
                        3,
                        4,
                        5,
                        myfilter.Filter(
                            (6, 7),
                        ),
                    ),
                ),
            ),
            2,
            (3, 4),
            6,
        )
        self.assertEqual(f1._pattern, f2._pattern)
        self.assertEqual(f1, f2)

        # Test that order doesn't matter
        f3 = myfilter.Filter(2, (3, 3, 4), 6, (1, 2, (3, 4, 5, (6, 7))))
        self.assertEqual(f1, f3)

        f4 = myfilter.Filter(1, 2)
        f5 = myfilter.Filter(2, 1)
        self.assertEqual(f4, f5)
        f6 = myfilter.Filter(2, (1,))
        self.assertEqual(f5, f6)

        # Test matching
        f7 = myfilter.Filter((1, 2, 3), (3, 4))
        f8 = myfilter.Filter((3, 4), 2)
        expected_overlap = [(2, 3), (2, 4)]
        self.assertEqual(f7.check_overlap(f8), expected_overlap)
        f8 = myfilter.Filter((3, 4), 4)
        expected_overlap = [(3, 4)]
        self.assertEqual(f7.check_overlap(f8), expected_overlap)

    def test_FilterInt(self):
        self.run_filters(myfilter.Filter, 1, 2, 3, 4)

    def test_FilterStr(self):
        self.run_filters(myfilter.Filter, "a", "b", "cd", "qwer")

    def test_FilterObj(self):
        self.run_filters(myfilter.Filter, A(1, 2), A(2, 3), A(4, 5), A(5, 6))

    # Next need to test a serializable (and typed) version of the filter

    # test FilterSet
    def test_FilterSet(self):
        self.run_filters(myfilter.FilterSet0, 1, 2, 3, 4)

        self.assertTrue(
            myfilter.eq_order_insensitive(((2, 1), (1, 2)), ((2, 1), (1, 2)))
        )
        self.assertTrue(
            myfilter.eq_order_insensitive(((2, 1, 3), (1, 2)), ((2, 1), (1, 2, 3)))
        )
        self.run_filters(myfilter.FilterSet, 1, 2, 3, 4)

    def test_TypedFilterSet(self):
        self.run_filters(myfilter.SerializableFilterSet, 1, 2, 3, 4)

        sfs0 = myfilter.SerializableFilterSet(myfilter.A(1), myfilter.A(2))
        d = sfs0.to_dict()
        sfs0.from_dict(d)
        self.assertTrue(sfs0.match(1, 2))
        self.assertTrue(sfs0.match((1,), (2,)))
        self.assertFalse(sfs0.match(1, 3))

        # to make a new SerializableFilterSet:
        sfs1 = myfilter.SerializableFilterSet()
        sfs1._oktype = int
        sfs2 = myfilter.SerializableFilterSet()
        sfs2._oktype = int
        sfs1.custom_set("pattern", [1])
        sfs2.custom_set("_pattern", [[1]])  # equivalent
        self.assertTrue(sfs1 == sfs2)
        sfs1.custom_set("_pattern", [(1,)])  # equivalent
        self.assertTrue(sfs1 == sfs2)

        sfs3 = myfilter.SerializableFilterSet()
        sfs3._oktype = myfilter.A
        sfs3.from_dict(d)
        self.assertTrue(sfs0 == sfs3)
        self.assertTrue(sfs3 == sfs0)


if __name__ == "__main__":
    unittest.main()
