"""
Todo:
    - verify that to_dict() is json-saveable, or at least tests with encoders to test json-saveability with more complicated data types.
"""
from __future__ import absolute_import, division, print_function

import unittest
from collections import OrderedDict

from mdfts.utils import serial


class TestSerialize(unittest.TestCase):
    def setUp(self):
        @serial.serialize(["x", "y"])
        class Srlz:
            def __init__(self, x, y):
                self.x = x
                self.y = y

            @classmethod
            def init_from_dict(cls, d, *args, **kwargs):
                try:
                    # obj = object.__new__(cls)
                    obj = cls(None, None, *args, **kwargs)
                except:
                    print("Could not initialize {}".format(cls))

                # final update
                obj.from_dict(d)
                return obj

        self.Srlz = Srlz

    def test_Serializable(self):
        # test that it tracks the variable indicated in the decorator
        s1 = self.Srlz(5, True)
        self.assertEqual(s1._serial_vars, ["x", "y"])

        # test that __eq__ works
        s2 = self.Srlz(6, False)
        self.assertNotEqual(s1, s2)

        # test that the to_dict() works
        data = s1.to_dict()
        self.assertEqual(data, OrderedDict([("x", 5), ("y", True)]))

        # test that from_dict()
        s2.from_dict(data)
        self.assertEqual(s1, s2)

    def test_SerializableTypedList(self):
        # test that the oktype works for non-serializable types, and type checking is working
        sl1 = serial.SerializableTypedList(int)
        with self.assertRaises(TypeError):
            sl1.append(1.0)

        # test that the oktype works for Serializable types, and type checking is working
        sl2 = serial.SerializableTypedList(self.Srlz)
        with self.assertRaises(TypeError):
            sl2.append(1.0)

        # check that list operations and to_dict() work
        sl2.append(self.Srlz(5, True))
        sl2.extend([self.Srlz(6, False), self.Srlz(7, True)])
        self.assertEqual(
            sl2.to_dict(),
            [
                OrderedDict([("x", 5), ("y", True)]),
                OrderedDict([("x", 6), ("y", False)]),
                OrderedDict([("x", 7), ("y", True)]),
            ],
        )

        # check that list operations and from_dict() work
        sl3 = serial.SerializableTypedList(self.Srlz)
        sl3.from_dict(sl2.to_dict())
        self.assertEqual(sl2.to_dict(), sl3.to_dict())

    def test_SerializbleTypedDict(self):
        # === Test different ways of initializing, and the has_many flag
        sd1 = serial.SerializableTypedDict()
        sd1.add_entry_type("sl1", int, has_many=False)
        sd1.add_entry_type("sl2", int, has_many=False)
        sd1.add_entry_type("sl3", int, has_many=True)
        sd1["sl1"] = 1
        sd1["sl2"] = 2
        sd1["sl3"] = serial.SerializableTypedList(int, *[7, 8, 9])

        sd2 = serial.SerializableTypedDict([("sl1", 1)])
        sd2["sl2"] = 2
        sd2.add_entry_type("sl3", int, has_many=True)
        sd2["sl3"].extend([7, 8, 9])
        # alternatively:
        # sd2["sl3"] = serial.SerializableTypedList(int)
        # sd2["sl3"].extend([7,8,9])
        # or, directly without first adding the entry type:
        # sd2["sl3"] = serial.SerializableTypedList(int, *[7, 8, 9])

        self.assertEqual(sd1, sd2)
        self.assertEqual(sd1._types, sd2._types)
        self.assertEqual(sd1._has_many, sd2._has_many)

        # === test that to_dict and from_dict work even for more complex objects
        sd3 = serial.SerializableTypedDict()
        sd3.add_entry_type("t1", self.Srlz, has_many=True)
        sd4 = serial.SerializableTypedDict()
        sd4.add_entry_type("t1", self.Srlz, has_many=True)

        sd3["t1"].extend([self.Srlz(5, True), self.Srlz(6, False), self.Srlz(7, True)])
        sd4.from_dict(sd3.to_dict())
        self.assertEqual(sd3, sd4)

        expected_data = OrderedDict(
            [
                (
                    "t1",
                    [
                        OrderedDict([("x", 5), ("y", True)]),
                        OrderedDict([("x", 6), ("y", False)]),
                        OrderedDict([("x", 7), ("y", True)]),
                    ],
                )
            ]
        )
        self.assertEqual(sd4.to_dict(), expected_data)

        # === test that from_dict does not add previously undefined keys to the schema
        expected_data["t2"] = 50
        sd4.from_dict(expected_data)

        self.assertEqual(sd3._types, sd4._types)
        self.assertEqual(len(sd3), len(sd4))
        self.assertNotEqual(sd4.to_dict(), expected_data)

        # === Test that update() throws error b/c it doesn't know how to cast basic types (e.g. dict) into more complicated types (in this case, the expected self.Srlz class).
        with self.assertRaises(TypeError):
            sd4.update(expected_data)

        # === Test schema copying
        sd5 = serial.SerializableTypedDict()
        # first validate that strict from_dict doesn't update any keys
        sd5.from_dict(sd4.to_dict())
        self.assertEqual(len(sd5), 0)
        # now copy schema and see that from_dict now works
        sd5.copy_schema(sd4)
        sd5.from_dict(sd4.to_dict())
        self.assertEqual(sd4, sd5)


if __name__ == "__main__":
    unittest.main()
