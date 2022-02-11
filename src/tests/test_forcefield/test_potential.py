from __future__ import absolute_import, division, print_function

import unittest

import math

import mdfts.forcefield as ff


class TestPotential(unittest.TestCase):

    def test_compare_bead_types(self):

        # create BeadTypes
        bead_type_A = ff.BeadType("A", 1.0, 1.0)
        bead_type_B = ff.BeadType("B", 2.0, -1.0)

        # create potentials and compare BeadTypes_PairPotential
        potential1 = ff.potential._PairPotential(bead_type_A, bead_type_B)
        potential2 = ff.potential._PairPotential(bead_type_B, bead_type_A)
        self.assertTrue(potential1.compare_bead_types(potential2))

    def test_Gaussian(self):

        # create BeadTypes
        bead_type_A = ff.BeadType("A", 1.0, 1.0)
        bead_type_B = ff.BeadType("B", 2.0, -1.0)

        # create Gaussian potential
        gaussian = ff.Gaussian(bead_type_A, bead_type_B)

        # check values of parameters
        self.assertEqual(0.0, gaussian.excl_vol.value)
        self.assertEqual(float, gaussian.excl_vol.type)
        self.assertEqual(False, gaussian.excl_vol.fixed)
        self.assertEqual(0.0, gaussian.B.value)
        self.assertEqual(float, gaussian.B.type)
        self.assertEqual(False, gaussian.B.fixed)
        self.assertEqual(0.1, gaussian.Kappa.value)
        self.assertEqual(float, gaussian.Kappa.type)
        self.assertEqual(True, gaussian.Kappa.fixed)

        # trying to set value of Kappa should raise error
        with self.assertRaises(NotImplementedError):
            gaussian.Kappa.value = 1.0

        # setting value of excl_vol should change B
        gaussian.excl_vol.value = 1.0
        self.assertEqual(1.0, gaussian.excl_vol.value)
        self.assertEqual(gaussian.excl_vol.value * (gaussian.Kappa.value / math.pi)**1.5, gaussian.B.value)

        # test serialization
        g_dict = gaussian.to_dict()
        g1 = ff.Gaussian(bead_type_A, bead_type_B)
        g1.from_dict(g_dict)
        self.assertEqual(g_dict, g1.to_dict())


if __name__ == '__main__':
    unittest.main()
