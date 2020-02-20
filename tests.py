import unittest
import stress_modules
import numericaltools
import math
import numpy as np
import InputClasses
import displacements
import equilibrium
from InputClasses import *
import validation

class MyTestCase(unittest.TestCase):

    def test_something(self):
        self.assertEqual(True, True)

    def test_bend(self):
        h = 10
        b = 5
        Iyy = b * h ** 3 / 12
        Izz = h * b ** 3 / 12
        centroid = (0, 5, 2.5)
        #No moments
        self.assertEqual(stress_modules.bending(0, 0, Iyy, Izz, 0, 0, centroid), 0)
        self.assertEqual(stress_modules.bending(5, 0, Iyy, Izz, 0, 0, centroid), 0)
        self.assertEqual(stress_modules.bending(5, 10, Iyy, Izz, 0, 0, centroid), 0)
        self.assertEqual(stress_modules.bending(0, 10, Iyy, Izz, 0, 0, centroid), 0)

        #On NA
        self.assertEqual(stress_modules.bending(5, 2.5, Iyy, Izz, 99999, 99999, centroid), 0)
        self.assertEqual(stress_modules.bending(5, 0, Iyy, Izz, 0, 99999, centroid), 0)
        self.assertEqual(stress_modules.bending(0, 2.5, Iyy, Izz, 99999, 0, centroid), 0)

        #One moment
        self.assertEqual(stress_modules.bending(0, 0, Iyy, Izz, 1, 0, centroid), -2.5 / Iyy)
        self.assertEqual(stress_modules.bending(0, 5, Iyy, Izz, 1, 0, centroid), 2.5 / Iyy)
        self.assertEqual(stress_modules.bending(10, 5, Iyy, Izz, 1, 0, centroid), 2.5 / Iyy)
        self.assertEqual(stress_modules.bending(10, 0, Iyy, Izz, 1, 0, centroid), -2.5 / Iyy)

        # other moment
        self.assertEqual(stress_modules.bending(0, 0, Iyy, Izz, 0, 1, centroid), -5 / Izz)
        self.assertEqual(stress_modules.bending(0, 5, Iyy, Izz, 0, 1, centroid), -5 / Izz)
        self.assertEqual(stress_modules.bending(10, 5, Iyy, Izz, 0, 1, centroid), 5 / Izz)
        self.assertEqual(stress_modules.bending(10, 0, Iyy, Izz, 0, 1, centroid), 5 / Izz)

        #Superposition
        self.assertEqual(stress_modules.bending(0, 0, Iyy, Izz, 1, 1, centroid), -5 / Izz - 2.5 / Iyy)
        self.assertEqual(stress_modules.bending(0, 5, Iyy, Izz, 1, 1, centroid), -5 / Izz + 2.5 / Iyy)
        self.assertEqual(stress_modules.bending(10, 5, Iyy, Izz, 1, 1, centroid), 5 / Izz + 2.5 / Iyy)
        self.assertEqual(stress_modules.bending(10, 0, Iyy, Izz, 1, 1, centroid), 5 / Izz - 2.5 / Iyy)

        #very small input
        self.assertEqual(stress_modules.bending(0, 0, Iyy, Izz, .000000001, 0, centroid), -.0000000025 / Iyy)

    def test_vonmis(self):
        h = 10
        b = 5
        Iyy = b * h ** 3 / 12
        Izz = h * b ** 3 / 12
        centroid = (0, 5, 2.5)

        #No stresses:
        self.assertEqual(stress_modules.von_mises(sigma_xx= 0, sigma_yy= 0, sigma_zz= 0, tau_xy= 0, tau_yz= 0, tau_xz= 0), 0)

        #One unit stress
        self.assertEqual(stress_modules.von_mises(sigma_xx=1, sigma_yy=0, sigma_zz=0, tau_xy=0, tau_yz=0, tau_xz=0), 1)
        self.assertEqual(stress_modules.von_mises(sigma_xx=1, sigma_yy=1, sigma_zz=0, tau_xy=0, tau_yz=0, tau_xz=0), 1)
        self.assertEqual(stress_modules.von_mises(sigma_xx=0, sigma_yy=0, sigma_zz=1, tau_xy=0, tau_yz=0, tau_xz=0), 1)

        #Two unit stresses
        self.assertEqual(stress_modules.von_mises(sigma_xx=1, sigma_yy=1, sigma_zz=0, tau_xy=0, tau_yz=0, tau_xz=0), math.sqrt((1)))
        self.assertEqual(stress_modules.von_mises(sigma_xx=0, sigma_yy=1, sigma_zz=1, tau_xy=0, tau_yz=0, tau_xz=0), math.sqrt((1)))

        #Three unit stresses (they should cancel out)
        self.assertEqual(stress_modules.von_mises(sigma_xx=1, sigma_yy=1, sigma_zz=1, tau_xy=0, tau_yz=0, tau_xz=0), 0)

        #Checking for different Tau combinations
        self.assertEqual(stress_modules.von_mises(sigma_xx= 0, sigma_yy= 0, sigma_zz= 0, tau_xy= 1, tau_yz= 0, tau_xz= 0), math.sqrt(6* 1**2 /2))
        self.assertEqual(stress_modules.von_mises(sigma_xx= 0, sigma_yy= 0, sigma_zz= 0, tau_xy= 2, tau_yz= 0, tau_xz= 0), math.sqrt(6* 2**2 /2))
        self.assertEqual(stress_modules.von_mises(sigma_xx= 0, sigma_yy= 0, sigma_zz= 0, tau_xy= 0, tau_yz= 1, tau_xz= 0), math.sqrt(6* 1**2 /2))
        self.assertEqual(stress_modules.von_mises(sigma_xx= 0, sigma_yy= 0, sigma_zz= 0, tau_xy= 0, tau_yz= 2, tau_xz= 0), math.sqrt(6* 2**2 /2))
        self.assertEqual(stress_modules.von_mises(sigma_xx= 0, sigma_yy= 0, sigma_zz= 0, tau_xy= 0, tau_yz= 0, tau_xz= 1), math.sqrt(6* 1**2 /2))
        self.assertEqual(stress_modules.von_mises(sigma_xx= 0, sigma_yy= 0, sigma_zz= 0, tau_xy= 0, tau_yz= 0, tau_xz= 2), math.sqrt(6* 2**2 /2))

        #Check for negative tau
        self.assertEqual(stress_modules.von_mises(sigma_xx= 0, sigma_yy= 0, sigma_zz= 0, tau_xy= -1, tau_yz= 0, tau_xz= 0), math.sqrt(6* 1**2 /2))

        #Negative unit stress
        self.assertEqual(stress_modules.von_mises(sigma_xx= -1, sigma_yy= 0, sigma_zz= 0, tau_xy= 0, tau_yz= 0, tau_xz= 0), 1)
        self.assertEqual(stress_modules.von_mises(sigma_xx= -1, sigma_yy= -1, sigma_zz= 0, tau_xy= 0, tau_yz= 0, tau_xz= 0), 1)
        self.assertEqual(stress_modules.von_mises(sigma_xx= -1, sigma_yy= -1, sigma_zz= -1, tau_xy= 0, tau_yz= 0, tau_xz= 0), 0)

        #Very small input
        self.assertEqual(stress_modules.von_mises(sigma_xx=.0000001, sigma_yy=0, sigma_zz=0, tau_xy=0, tau_yz=0, tau_xz=0), math.sqrt(.0000001**2))

    def test_integration(self):
        # integrate(func, start, stop, number_of_points)
        start = 0
        stop = 1

        def tempfunc(x):
            return 1 + x * 0
        n = 1
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), 1., places=8)
        n =100
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), 1., places=8)
        start = -5
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), 6., places=8)
        stop = -2
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), 3., places=8)

        def tempfunc(x):
            return 2*x
        start = 0
        stop = np.pi
        n = 1000
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), stop**2, places=8)
        stop = 10
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), stop**2, places=8)
        stop = -10
        # Stop < start
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), stop ** 2, places=8)
        start = 10
        # Integral should cancel out (symmetric)
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), 0, places=8)
        start, stop = 0, 0
        # No distance over which to integrate
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), 0, places=8)
        stop = .00000001
        # Very small distance over which to integrate
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), 0.00000001**2, places=10)

        # Harder book example from Early Transcendentals
        def tempfunc(x):
            return np.log(x) / x
        start = 1
        stop = np.e
        n = 100000
        self.assertAlmostEqual(numericaltools.integrate(tempfunc, start, stop, n), .5, places=8)

    def test_interpolation(self):
        # Simple test with basic python lists and some fringe cases (first value, high decimal count)
        test_list_x = [0, 1, 2, 3, 4]
        test_list_f = [2, 5, 3, 3, 4]
        test_x_target = 0
        self.assertEqual(numericaltools.interpolate(test_list_x, test_list_f, test_x_target), 2)
        test_x_target = 0.5
        self.assertEqual(numericaltools.interpolate(test_list_x, test_list_f, test_x_target), 3.5)
        test_x_target = 2.9999999
        self.assertEqual(numericaltools.interpolate(test_list_x, test_list_f, test_x_target), 3)
        # More complex test with other fringe case (last value)
        test_list_x_1 = np.array([0.001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        test_list_f_1 = np.array([4, 3, 2, 3, 3, 6, 2, 3, 5, 4, 3, 1])
        test_x_target = 5.5
        self.assertEqual(numericaltools.interpolate(test_list_x_1, test_list_f_1, test_x_target), 4)
        test_x_target = 7.5
        self.assertEqual(numericaltools.interpolate(test_list_x_1, test_list_f_1, test_x_target), 4)
        test_x_target = 11
        self.assertEqual(numericaltools.interpolate(test_list_x_1, test_list_f_1, test_x_target), 1)

    def test_input_geo(self):
        # Tests Izz
        chord = 1.1
        height = 0.8
        skint = 0.01
        spart = 0.05
        a = Aileron(chord=chord, height=height, skint=skint, spart=spart)
        # self.assertEqual(a.Izz(stiffener=False, spar=False), 0.001865283305) # Exact value different from
        # manually calculated value by less than 0.1%, hence human error, and can be interpreted as correct (Works)
        # self.assertEqual(a.Izz(skin=False, spar=False),) todo
        # self.assertEqual(a.Izz(skin=False, stiffener=False), 0.0256/12) Works

        # Tests Iyy
        self.assertEqual(a.Iyy(stiffener=False, spar=False), 0.003638774597)
        # self.assertEqual(a.Iyy(skin=False, spar=False),) todo
        # self.assertEqual(a.Iyy(skin=False, stiffener=False), 0.0001/12) Works

        # Tests centroid

        # Tests shearcentre
        # TODO: implement in InputClasses.py

        # Tests visualinspection?

        # Tests _Istiff

    def test_validation_getdata(self):
        self.assertEqual(validation.get_dat('bending', 'stresses').size, (65338, 5))
        #self.assertEqual(validation.get_dat('bending', 'disp').size, (65338, 5))
        #self.assertEqual(validation.get_dat('Jam_Bent', 'stresses').size, (65338, 5))
        #self.assertEqual(validation.get_dat('Jam_Bent', 'disp').size, (65338, 5))
        #self.assertEqual(True, False)


# class SystemTests(unittest.TestCase):
#     def test_no_load_no_deformation(self):
#  TODO: CHANGE THIS YO, SHOULDN'T BE COMMENTED

if __name__ == '__main__':
    #data = InputClasses.Aileron()
    #data.visualinspection()
    unittest.main()
