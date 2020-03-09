import InputClasses
import main, Energy
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import sys
import unittest
import InputClasses
import numericaltools
import numpy as np
aileron = InputClasses.Aileron()
#app_loads = InputClasses.AppliedLoads()
import equilibrium
mainfile = main

class MyTestCase(unittest.TestCase):
    def test_Izz_Iyy(self):
        self.assertAlmostEqual(mainfile.Izz, aileron.Izz(), 9)
        self.assertAlmostEqual(mainfile.Iyy, aileron.Iyy(), 9)

    def test_centroid(self):
        self.assertAlmostEqual(mainfile.yc, aileron.centroid(1), 9)
        self.assertAlmostEqual(mainfile.zc+aileron.radius, aileron.centroid(2), 9)

    def test_shear_centre(self):
        self.assertAlmostEqual(mainfile.zsc, aileron.shearcentre(self)[0], 9)

    def test_stiffloc(self):
        st_list = aileron.stiffLoc()
        starray = np.array(st_list)
        starray, stmainfile = np.sort(starray, axis = 0), np.sort(mainfile.stcoord, axis = 0)
        for i in range(len(starray)):
            self.assertAlmostEqual(starray[i][0]-aileron.radius, stmainfile[i][0], 9)
            self.assertAlmostEqual(starray[i][1], stmainfile[i][1], 9)

   # def test_stiff_s_pos(self):

class SystemTests(unittest.TestCase):
   def test_no_load_no_deformation(self):
       self.assertIsNone(equilibrium.calc_deflection_y(loads=0, u = np.nan))
       self.assertIsNone(equilibrium.calc_deflection_z(loads=0, w = np.nan))

   def test_unit_loads_applied(self):
       self.assertAlmostEqual(equilibrium.calc_deflection_y(loads=1, u = np.nan))
       self.assertAlmostEqual(equilibrium.calc_deflection_z(loads=1, w = np.nan))

   def test_normal_loads_compared_with_hand_calculations_at_different_points(self):
       self.assertAlmostEqual(equilibrium.calc_deflection_y(loads=np.nan, u=np.nan))
       self.assertAlmostEqual(equilibrium.calc_deflection_y(loads=np.nan, u=np.nan))
       self.assertAlmostEqual(equilibrium.calc_deflection_y(loads=np.nan, u=np.nan))

       self.assertAlmostEqual(equilibrium.calc_deflection_z(loads=np.nan, w=np.nan))
       self.assertAlmostEqual(equilibrium.calc_deflection_z(loads=np.nan, w=np.nan))
       self.assertAlmostEqual(equilibrium.calc_deflection_z(loads=np.nan, w=np.nan))

   def test_check_hinge_2_deflection_0_at_different_loads(self):
        self.assertIsNone(displacements.displ(some load))
        self.assertIsNone(displacements.displ(some load))
        self.assertIsNone(displacements.displ(some load))
        self.assertIsNone(displacements.displ(some load))

   def test_check_large_E_G_give_small_displacements(self):
       Elarge =
       Evlarge =
       Glarge =
       Gvlarge =

        self.assertLess(displacements.displ(normal loads, large E), some low displacement value,))
        self.assertLess(displacements.displ(normal loads, very very large E), some super low displacement value,))

        self.assertLess(displacements.angle_displ(normal loads, large G), some low angle value,))
        self.assertLess(displacements.angle_displ(normal loads, very very large G), some suprt low angle value,))


if __name__ == '__main__':
#    print(mainfile.Izz, mainfile.Iyy)
#    print(main.d1, main.e1)
#    print(mainfile.stcoord)
#    print("aaa",np.array(aileron.stiffLoc(3)))
    unittest.main()

