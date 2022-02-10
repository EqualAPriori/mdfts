from __future__ import absolute_import, division, print_function

import unittest

import mdfts.forcefield as ff


def create_test_forcefield_1():
    f = ff.ForceField()
    f.add_bead_type(ff.BeadType("A", 2.0, -1.0))
    f.add_bead_type(ff.BeadType("B", 5.0, 1.0))
    return f


class TestForceField(unittest.TestCase):

    def test_forcefield_1_bead_types(self):

        # create forcefield
        f = create_test_forcefield_1()

        # check that attributes of bead types are correct
        bead_type_A = f.get_bead_type("A")
        self.assertEqual("A", bead_type_A.name)
        self.assertEqual(2.0, bead_type_A.smear_length)
        self.assertEqual(-1.0, bead_type_A.charge)
        bead_type_B = f.get_bead_type("B")
        self.assertEqual("B", bead_type_B.name)
        self.assertEqual(5.0, bead_type_B.smear_length)
        self.assertEqual(1.0, bead_type_B.charge)

        # check reordering
        self.assertEqual(["A", "B"], list(f.bead_names))
        f.reorder_bead_types(["B", "A"])
        self.assertEqual(["B", "A"], list(f.bead_names))

    def test_forcefield_1_reordering(self):

        # create forcefield
        f = create_test_forcefield_1()

        # check reordering
        self.assertEqual(["A", "B"], list(f.bead_names))
        f.reorder_bead_types(["B", "A"])
        self.assertEqual(["B", "A"], list(f.bead_names))

        # ensure that correct error is thrown
        with self.assertRaises(ValueError):
            f.reorder_bead_types(["C"])

    def test_forcefield_1_serialization(self):

        # create forcefield
        f = create_test_forcefield_1()

        # check that forcefield serializes correctly
        f_dict = f.to_dict()
        self.assertEqual(['_kT', '_bead_types', '_potentials'], list(f_dict.keys()))
        self.assertEqual(1.0, f_dict['_kT'])
        bead_type_A = f_dict['_bead_types'][0]
        self.assertEqual(bead_type_A['name'], "A")
        self.assertEqual(bead_type_A['smear_length'], 2.0)
        self.assertEqual(bead_type_A['charge'], -1.0)
        bead_type_B = f_dict['_bead_types'][1]
        self.assertEqual(bead_type_B['name'], "B")
        self.assertEqual(bead_type_B['smear_length'], 5.0)
        self.assertEqual(bead_type_B['charge'], 1.0)
        self.assertEqual([], list(f_dict['_potentials']))

        # check that we can initialize a new forcefield using serialized dictionary
        g = ff.ForceField()
        g.from_dict(f_dict)
        g_dict = g.to_dict()
        self.assertEqual(['_kT', '_bead_types', '_potentials'], list(g_dict.keys()))
        self.assertEqual(1.0, g_dict['_kT'])
        bead_type_A = g_dict['_bead_types'][0]
        self.assertEqual(bead_type_A['name'], "A")
        self.assertEqual(bead_type_A['smear_length'], 2.0)
        self.assertEqual(bead_type_A['charge'], -1.0)
        bead_type_B = g_dict['_bead_types'][1]
        self.assertEqual(bead_type_B['name'], "B")
        self.assertEqual(bead_type_B['smear_length'], 5.0)
        self.assertEqual(bead_type_B['charge'], 1.0)
        self.assertEqual([], list(g_dict['_potentials']))


if __name__ == '__main__':
    unittest.main()
