from __future__ import absolute_import, division, print_function

import unittest
from collections import OrderedDict

from mdfts.utils import topology


def make_sample_2gen():
    ft = topology.FTSTopology()
    u = ft.add_path(
        [("Aaaa", 10)],
        [
            (("Bbbb", 3), ("Cccc", 2)),
            [1, 1, 8, 8],
            [[("Aaaa", 2)], [0, 2, 4]],
        ],
        mode=1,
        as_chain_name="Ch",
    )
    return ft


sample_2gen_verbose = (
    [
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Bbbb",
        "Bbbb",
        "Bbbb",
        "Cccc",
        "Cccc",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Bbbb",
        "Bbbb",
        "Bbbb",
        "Cccc",
        "Cccc",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Bbbb",
        "Bbbb",
        "Bbbb",
        "Cccc",
        "Cccc",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Bbbb",
        "Bbbb",
        "Bbbb",
        "Cccc",
        "Cccc",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
        "Aaaa",
    ],
    [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (1, 10),
        (10, 11),
        (11, 12),
        (12, 13),
        (13, 14),
        (10, 15),
        (15, 16),
        (12, 17),
        (17, 18),
        (14, 19),
        (19, 20),
        (1, 21),
        (21, 22),
        (22, 23),
        (23, 24),
        (24, 25),
        (21, 26),
        (26, 27),
        (23, 28),
        (28, 29),
        (25, 30),
        (30, 31),
        (8, 32),
        (32, 33),
        (33, 34),
        (34, 35),
        (35, 36),
        (32, 37),
        (37, 38),
        (34, 39),
        (39, 40),
        (36, 41),
        (41, 42),
        (8, 43),
        (43, 44),
        (44, 45),
        (45, 46),
        (46, 47),
        (43, 48),
        (48, 49),
        (45, 50),
        (50, 51),
        (47, 52),
        (52, 53),
    ],
)


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
        FT3 = make_sample_2gen()

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

    def test_expansion(self):
        # Test expansion to beads
        FT = make_sample_2gen()
        self.assertEqual(FT.expand_graft_to_beads(0), sample_2gen_verbose)

        FT2 = FT.fully_enumerate()
        self.assertEqual(FT.expand_to_beads(), FT2.expand_to_beads())


if __name__ == "__main__":
    unittest.main()
