import InputClasses
from A10_SVV_VerificationModel import main, Energy
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import sys
import unittest
import InputClasses
mainfile = main
aileron = InputClasses.Aileron()


class MyTestCase(unittest.TestCase):
#    def test_Izz_Iyy(self):
#        self.assertAlmostEqual(mainfile.Izz, aileron.Izz(), 9)
#        self.assertAlmostEqual(mainfile.Iyy, aileron.Iyy(), 9)

    def test_centroid(self):
#        self.assertAlmostEqual(mainfile.yc, aileron.centroid(1), 9)
        self.assertAlmostEqual(mainfile.zc, aileron.centroid(2), 9)

if __name__ == '__main__':
#    print(mainfile.Izz, mainfile.Iyy)
    print(aileron.centroid())
#    print(mainfile.stcoord,mainfile.totarea)
#    print(aileron._stiffcoord(7)[0]-aileron.radius)
    unittest.main()
