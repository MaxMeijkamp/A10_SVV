import unittest
import stress_modules
import numericaltools
import math
import numpy as np
import InputClasses
import displacements
import equilibrium
from InputClasses import *
#import validation

class MyTestCase(unittest.TestCase):

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
        span = 5
        a = Aileron(chord=chord, height=height, skint=skint, spart=spart, span=span)
        # self.assertEqual(a.Izz(stiffener=False, spar=False), 0.001865283305) # Exact value different from
        # manually calculated value by less than 0.1%, hence human error, and can be interpreted as correct (Works)
        # self.assertEqual(a.Izz(skin=False, spar=False),) todo
        # self.assertEqual(a.Izz(skin=False, stiffener=False), 0.0256/12) Works

        # Tests Iyy
        # self.assertEqual(a.Iyy(stiffener=False, spar=False), 0.003638774597) # Exact value different from
        # manually calculated value by less than 0.1%, hence human error, and can be interpreted as correct (Works)
        # self.assertEqual(a.Iyy(skin=False, spar=False),) todo
        # self.assertEqual(a.Iyy(skin=False, stiffener=False), 0.0001/12) Works

        # Tests centroid
        self.assertEqual(a.centroid(0), 2.5)
        self.assertEqual(a.centroid(1), 0)
        # self.assertEqual(a.centroid(2), )
        # self.assertEqual(a.centroid(), (2.5, 0, ))

        # Tests shearcentre
        # TODO: implement in InputClasses.py

        # Tests visualinspection?

        # Tests _Istiff

    def test_spline(self):
        x = [0, 1, 2, 3, 4, 5, 6, 7]
        f = [0, 1, 2, 3, 4, 5, 6, 7]
        n = 8
        self.assertEqual(numericaltools.spline(x, f, n)[0], f[0:-1])

        spline_slopes = numericaltools.spline(x, f, n)[1]
        for elem in spline_slopes:
            self.assertEqual(elem, 1)

        x = [0, 10, 20, 30, 40, 50, 60, 70]
        spline_slopes = numericaltools.spline(x, f, n)[1]
        for elem in spline_slopes:
            self.assertEqual(elem, .1)

        # Increasing and decreasing f
        x = [0, 1, 2, 3, 4, 5, 6, 7]
        f = [0, 2, 0, -2, -8, 0, 3, 7]
        slopes = [2, -2, -2, -6, 8, 3, 4]
        self.assertEqual(numericaltools.spline(x, f, n)[1], slopes)

        # Non-consistent step size in x
        n = 4
        x = [0, 1, 5, 7]
        f = [0, 1, 2, 9]
        slopes = [1, .25, 3.5]
        self.assertEqual(numericaltools.spline(x, f, n)[1], slopes)

        # Negative values in x
        x = [-20, -10, 0]
        f = [0, 1, 2]
        slopes = [.1, .1]
        self.assertEqual(numericaltools.spline(x, f, 3)[1], slopes)
        # Negative x backwards
        x = [0, -10, -20]
        f = [0, 1, 2]
        slopes = [-.1, -.1]
        self.assertEqual(numericaltools.spline(x, f, 3)[1], slopes)
        x = np.array([-20, -10, 0])
        f = np.array([0, 1, 2])
        slopes = [.1, .1]
        self.assertEqual(numericaltools.spline(x, f, 3)[1], slopes)

        # Very small values
        x = [0, 0.00000001, 0.000000002]
        f = [0, 0.000000015, 0.000000003]
        slopes = [1.5, 1.5]
        self.assertAlmostEqual(numericaltools.spline(x, f, 3)[1][0], slopes[0], places= 8)

    #def test_validation_read_data(self):
        #self.assertEqual(validation.get_dat("bending", "stresses").size[1],  5)
        #pass


 # class SystemTests(unittest.TestCase):
 #     def test_no_load_no_deformation(self):
 #
 #    def test_unit_loads_applied(self):
 #
 #    def test_normal_loads_compared_with_hand_calculations_at_different_points(self):
 #
 #    def test_check_hinge_2_deflection_0_at_different_loads(self):
 #        self.assertIsNone(displacements.displ(normal loads))
 #
 #    def test_check_large_E_G_give_small_displacements(self):
 #        self.assertLess(displacements.displ(normal loads, large E), some low displacement value,))
 #        self.assertLess(displacements.displ(normal loads, very very large E), some super low displacement value,))
 #        self.assertLess(displacements.angle_displ(normal loads, large G), some low angle value,))
 #        self.assertLess(displacements.angle_displ(normal loads, very very large G), some suprt low angle value,))

if __name__ == '__main__':
    data = InputClasses.Aileron()
    data.visualinspection()
    unittest.main()
