from numericaltools import integrate
import unittest

def func(x):
    y = 2*x + 1
    return y

print(integrate(func, 0, 6, 3))

class SimpleTest(unittest.TestCase):

    def setUp(self):
        pass
    def test_integrate(self):
        self.assertEqual( integrate(func, 0, 6, 3), 42)

#Other things to test: Negative number of steps, 0 steps, stepsize negative, reaction to different start and stops,

if __name__ == '__main__':
    unittest.main()