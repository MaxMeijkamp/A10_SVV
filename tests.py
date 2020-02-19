import unittest
import stress_modules

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

        from numericaltools import integrate
        import unittest

        def func(x):
            y = 2 * x + 1
            return y

        print(integrate(func, 0, 6, 3))

        def setUp(self):
            pass

        def test_integrate(self):
            self.assertEqual(integrate(func, 0, 6, 3), 42)

        # Other things to test: Negative number of steps, 0 steps, stepsize negative, reaction to different start and stops,

if __name__ == '__main__':
    unittest.main()
