from __future__ import absolute_import, division, print_function

import unittest
from collections import OrderedDict

from mdfts.utils import topology


class TestTopology(unittest.TestCase):
    def setUp(self):
        pass

    def test_segment(self):
        s = topology.Segment(["A", "B", ("A", 3)])
        self.assertEqual(s.sequence(), ("A", "B", "A", "A", "A"))
        self.assertEqual(s.sequence_compact(), (("A", 1), ("B", 1), ("A", 3)))
        self.assertEqual(s.sequence_compactest(), ("A", "B", ("A", 3)))
        self.assertEqual(
            s.to_dict(),
            OrderedDict(
                [
                    ("statistics", "DGC"),
                    ("block_species", ["A", "B", "A"]),
                    ("n_per_block", [1, 1, 3]),
                ]
            ),
        )
        s1 = topology.Segment(["A", "A", "A"])
        self.assertEqual(s1.block_species, ["A", "A", "A"])
        self.assertEqual(s1.n_per_block, [1, 1, 1])
        s1.to_compact()
        self.assertEqual(s1.block_species, ["A"])
        self.assertEqual(s1.n_per_block, [3])

    def test_chain_and_top(self):
        FT = topology.FTSTopology()
        c = FT.add_chaintype([("A", 10)])
        a1 = FT.add_armtype([("B", 3), ("C", 2)])
        FT.add_graft(c, a1, [1, 2, 4, 8], 0, multiplicity=2)
        # ^^ is equivalent to: FT.add_graft(c,a1,[1,1,2,2,4,4,8,8],0,multiplicity = 1)
        a2 = FT.add_armtype([("A", 2)])
        FT.add_graft(a1, a2, [0, 2, 4])
        # FT.visualize()
        # topology.plt.show()
        self.assertTrue(FT.isvalid())  # True, i.e. no infinite recursion chains
        FT.add_graft(a1, a2, [0, 2, 4])
        self.assertTrue(FT.isvalid())  # True, i.e. no infinite recursion chains
        FT.add_graft(a2, a1, [0], -1)
        self.assertFalse(FT.isvalid())

        # Check Serialization
        f = topology.FTSTopology.init_from_dict(FT.to_dict())
        self.assertEqual(f.to_dict(), FT.to_dict())

        # Check alternative equivalent definitions
        FT1 = topology.FTSTopology()
        c = FT1.add_chaintype([("A", 5), ("A", 5)])
        a1 = FT1.add_armtype([("B", 3), ("C", 2)])
        FT1.add_graft(c, a1, [1, 1, 2, 2, 4, 4, 8, 8], 0, multiplicity=1)
        a2 = FT1.add_armtype(["A", "A"])
        FT1.add_graft(a1, a2, [0, 2, 4])

        # Shorthand, should be same as above!
        FT2 = topology.FTSTopology()
        FT2.add_path(
            [("Aaaa", 10)],
            [
                (("Bbbb", 3), ("Cccc", 2)),
                [1, 1, 2, 2, 4, 4, 8, 8],
                [[("Aaaa", 2)], [0, 2, 4]],
            ],
            mode=1,
        )

        # Shorthand, with multiple graft types
        FT3 = topology.FTSTopology()
        FT3.add_path(
            [("Aaaa", 10)],
            [  # graft type 1
                (("Bbbb", 3), ("Cccc", 2)),
                [1, 1, 2, 2, 4, 4, 8, 8],
                [[("Aaaa", 2)], [0, 2, 4]],  # subgrafts off of graft type 1!
            ],
            [(("Dddd", 3)), [5, 5, 4, 3, 2]],  # graft type 2
            mode=1,  # make new arm_types as we go
        )

        # Test Enumeration
        FT3 = topology.FTSTopology()
        u = FT3.add_path(
            [("Aaaa", 10)],
            [
                (("Bbbb", 3), ("Cccc", 2)),
                [1, 1, 2, 2, 4, 4, 8, 8],
                [[("Aaaa", 2)], [0, 2, 4]],
            ],
            mode=1,
            as_chain_name="Ch",
        )

        FT4 = topology.FTSTopology()
        graft_def = FT3.get_grafts(0)
        u = FT4.add_path(FT3.arm_types[0], *graft_def, mode=1, as_chain_name="Ch")

        self.assertNotEqual(FT4.to_dict(), FT3.to_dict())

        FT5 = FT4.fully_enumerate()
        self.assertEqual(FT5.to_dict(), FT4.to_dict())

        # Read write
        FT5.save()
        FT6 = topology.FTSTopology()
        FT6.load()
        self.assertEqual(FT5.to_dict(), FT6.to_dict())


if __name__ == "__main__":
    unittest.main()
