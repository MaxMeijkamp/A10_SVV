import InputClasses
from A10_SVV_VerificationModel import main, Energy
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import sys
import unittest
import InputClasses
import numericaltools
import numpy as np
mainfile = main
aileron = InputClasses.Aileron()


class MyTestCase(unittest.TestCase):
    def test_Izz_Iyy(self):
        self.assertAlmostEqual(mainfile.Izz, aileron.Izz(), 9)
        self.assertAlmostEqual(mainfile.Iyy, aileron.Iyy(), 9)

    def test_centroid(self):
        self.assertAlmostEqual(mainfile.yc, aileron.centroid(1), 9)
        self.assertAlmostEqual(mainfile.zc+aileron.radius, aileron.centroid(2), 9)

    def test_stiffloc(self):
        st_list = aileron.stiffLoc()
        for stiffener in st_list:
            stiffener = np.array(stiffener)
            stiffener[0] += aileron.radius
        starray = np.array(st_list)
        starray, stmainfile = np.sort(starray, axis = 0), np.sort(mainfile.stcoord, axis = 0)
        for i in range(len(starray)):
            self.assertEqual(starray[i], stmainfile[i])

if __name__ == '__main__':
#    print(mainfile.Izz, mainfile.Iyy)
#    print(main.d1, main.e1)
    print(mainfile.stcoord)
    print("aaa",np.array(aileron.stiffLoc(3)))
    unittest.main()

