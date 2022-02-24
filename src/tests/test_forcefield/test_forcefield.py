from __future__ import absolute_import, division, print_function

import unittest

import os

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

    def test_import_from_sim(self):
        kT = 313.15 * 8.314462618 / 1000
        ff_path = os.path.join("sim_forcefields/pba_313.15K_ff.dat")
        f = ff.load_from_sim_ff(ff_path, kT=kT)
        gaussian_D4_D4 = f.get_potential(ff.Gaussian, 'D4', 'D4')
        self.assertAlmostEqual(2.84325299791, gaussian_D4_D4.excl_vol.value)
        self.assertAlmostEqual(4.4444e-01, gaussian_D4_D4.Kappa.value)
        self.assertAlmostEqual(3.9391e-01 / kT, gaussian_D4_D4.B.value)
        gaussian_Bpba_D4 = f.get_potential("Gaussian", "D4", "Bpba")
        self.assertAlmostEqual(2.1099736395319293, gaussian_Bpba_D4.excl_vol.value)
        self.assertAlmostEqual(4.4444e-01, gaussian_Bpba_D4.Kappa.value)
        self.assertAlmostEqual(2.9232e-01 / kT, gaussian_Bpba_D4.B.value)
        gaussian_Bpba_Bpba = f.get_potential("Gaussian", "Bpba", "Bpba")
        self.assertAlmostEqual(1.3695484344717714, gaussian_Bpba_Bpba.excl_vol.value)
        self.assertAlmostEqual(4.4444e-01, gaussian_Bpba_Bpba.Kappa.value)
        self.assertAlmostEqual(1.8974e-01 / kT, gaussian_Bpba_Bpba.B.value)
        harmonic_bond_Bpba_Bpba = f.get_potential(ff.HarmonicBond, "Bpba", "Bpba")
        self.assertAlmostEqual(3.8718e-01, harmonic_bond_Bpba_Bpba.r0.value)
        self.assertAlmostEqual(1.7903e+03 / kT, harmonic_bond_Bpba_Bpba.K.value)


if __name__ == '__main__':
    unittest.main()
