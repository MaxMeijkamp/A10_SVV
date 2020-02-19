import unittest
import stress_modules
import numericaltools
import math
import numpy as np

#simple case: rectangle with y = 10, x = 5
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

    def test_integration(self):
        #integrate(func, start, stop, number_of_points)
        start = 0
        stop = 1
        def tempfunc(x):
            return 1. + x * 0.
        n = 1
        self.assertAlmostEqual(numericaltools.integrateP(tempfunc, start, stop, n), 1., places= 8)
        n =100
        self.assertAlmostEqual(numericaltools.integrateP(tempfunc, start, stop, n), 1., places=8)
        start = -5
        self.assertAlmostEqual(numericaltools.integrateP(tempfunc, start, stop, n), 6., places=8)
        stop = -2
        self.assertAlmostEqual(numericaltools.integrateP(tempfunc, start, stop, n), 3., places=8)

        def tempfunc(x):
            return 2*x
        start = 0
        stop = np.pi
        n = 1000
        self.assertAlmostEqual(numericaltools.integrateP(tempfunc, start, stop, n), stop**2, places=8)
        stop = 10
        self.assertAlmostEqual(numericaltools.integrateP(tempfunc, start, stop, n), stop**2, places=8)
        stop = -10
        self.assertAlmostEqual(numericaltools.integrateP(tempfunc, start, stop, n), stop ** 2, places=8)
        start = 10
        self.assertAlmostEqual(numericaltools.integrateP(tempfunc, start, stop, n), 0, places=8)


if __name__ == '__main__':
    unittest.main()
