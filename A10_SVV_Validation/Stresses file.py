#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 13:34:55 2020

@author: mustafawahid
"""


import numpy as np
import matplotlib.pyplot as plt

"Von Mises Bending Values"


a = np.array([0.21416849999999998, 0.2491405, 0.275432, 0.3018725, 0.3160955, 0.34132949999999995, 0.35193850000000004, 0.250922, 0.21329199999999998, 0.18317, 0.15722, 0.13265300000000002, 0.1086245, 0.08616805, 0.06534845, 0.050042249999999996, 0.044589050000000005, 0.0533709, 0.07051815, 0.09196105, 0.11489450000000001, 0.138847, 0.163308, 0.1882005, 0.21488000000000002, 0.2438735, 0.2960785, 0.06448155, 0.01781195, 0.039232050000000004, 0.08005465, 0.11964149999999998, 0.1628055, 0.2084945, 0.2857475, 0.10112080000000001, 0.108149, 0.1171475, 0.125912, 0.135342, 0.144457, 0.15379500000000002, 0.1629555, 0.1721955, 0.18135649999999998, 0.19052200000000002, 0.199684, 0.208844, 0.218179, 0.228055, 0.2387625, 0.2541775, 0.27699, 0.3617365, 0.09582925, 0.09399555, 0.1731305, 0.1321405, 0.09707725, 0.07126735000000001, 0.06569354999999999, 0.08091860000000001])

b = np.array([(1198.5, -13.2644758, 100.75370025000001), (1198.5, 13.2644758, 100.75370025000001), (1198.5, -100.75370025000001, 13.2644758), (1198.5, -100.01399229999998, -12.1875), (1198.5, -89.6875, 0.0), (1198.5, 89.6875, 0.0), (1198.5, 100.0139923, -12.1875), (1198.5, 100.75370025000001, 13.2644758), (1198.5, 90.120285, -60.69078825), (1198.5, 95.04197690000001, -36.5625), (1198.5, -2.4356840876332786, -490.559204), (1198.5, 2.4356840823667216, -490.559204), (1198.5, -90.120285, -60.69078825), (1198.5, -95.0419769, -36.5625), (1198.5, -38.8894758, 93.88750075), (1198.5, -61.86422350000001, 80.623024), (1198.5, -80.623024, 61.8642235), (1198.5, -93.88750075, 38.8894758), (1198.5, 93.88750075, 38.8894758), (1198.5, 80.623024, 61.8642235), (1198.5, 61.86422350000001, 80.623024), (1198.5, 38.8894758, 93.88750075), (1198.5, 7.307050824999999, -466.67762749999997), (1198.5, 12.178417465, -442.79605100000003), (1198.5, 17.049784199999998, -418.9144745), (1198.5, 21.921150675, -395.03289800000005), (1198.5, 26.792517175, -371.1513215), (1198.5, 31.663883675, -347.269745), (1198.5, 36.535251599999995, -323.388153), (1198.5, 41.4066181, -299.5065765), (1198.5, 46.277984599999996, -275.625), (1198.5, 51.149353000000005, -251.74341600000002), (1198.5, 56.020719500000006, -227.8618395), (1198.5, 60.892086000000006, -203.9802625), (1198.5, 65.7634525, -180.098686), (1198.5, 70.634819, -156.2171095), (1198.5, 75.5061874, -132.335529), (1198.5, 80.377553925, -108.45394885), (1198.5, 85.248918525, -84.5723686), (1198.5, -7.307050825, -466.6776275), (1198.5, -12.178417465, -442.79605100000003), (1198.5, -17.049784199999998, -418.9144745), (1198.5, -21.921150675, -395.03289800000005), (1198.5, -26.792517175, -371.1513215), (1198.5, -31.663883675, -347.269745), (1198.5, -36.535251599999995, -323.388153), (1198.5, -41.4066181, -299.5065765), (1198.5, -46.277984599999996, -275.625), (1198.5, -51.149353000000005, -251.74341600000002), (1198.5, -56.0207195, -227.8618395), (1198.5, -60.892086000000006, -203.9802625), (1198.5, -65.7634525, -180.09868600000001), (1198.5, -70.634819, -156.2171095), (1198.5, -75.5061874, -132.335529), (1198.5, -80.377553925, -108.45394885), (1198.5, -85.24891852500001, -84.5723686), (1198.5, 64.0625, 0.0), (1198.5, 38.4375, 0.0), (1198.5, 12.8125, 0.0), (1198.5, -12.8125, 0.0), (1198.5, -38.4375, 0.0), (1198.5, -64.0625, 0.0)])


"Von Mises Jam Bent"
c = np.array([0.2205195, 0.2593535, 0.283264, 0.30415749999999997, 0.3117065, 0.3013935, 0.1744705, 0.128961, 0.09977849999999999, 0.09882455000000001, 0.116076, 0.13419599999999998, 0.13574000000000003, 0.1725805, 0.184198, 0.1960305, 0.203196, 0.2097135, 0.215689, 0.2216505, 0.22752450000000002, 0.233744, 0.2398195, 0.24662499999999998, 0.25303549999999997, 0.260871, 0.267673, 0.277269, 0.284365, 0.2976005, 0.306301, 0.33610249999999997, 0.4141685, 0.0598689, 0.08053975, 0.0902868, 0.1142785, 0.1371515, 0.167458, 0.1925595, 0.2295005, 0.24917850000000002, 0.2697815, 0.370124, 0.244205, 0.203584, 0.1751405, 0.148788, 0.125356, 0.10117380000000001, 0.07997099999999999, 0.060407050000000004, 0.0494782, 0.049096100000000004, 0.06206475, 0.0807905, 0.10327615, 0.126703, 0.1513065, 0.1760825, 0.2015835, 0.22553])


d = np.array([(1223.5, -13.2644758, 100.75370025000001), (1223.5, 13.2644758, 100.75370025000001), (1223.5, -100.75370025000001, 13.2644758), (1223.5, -89.6875, 0.0), (1223.5, -100.01399229999998, -12.1875), (1223.5, 100.75370025000001, 13.2644758), (1223.5, 100.0139923, -12.1875), (1223.5, 89.6875, 0.0), (1223.5, 95.04197690000001, -36.5625), (1223.5, 90.120285, -60.69078825), (1223.5, 2.4356840823667216, -490.559204), (1223.5, -2.4356840876332786, -490.559204), (1223.5, -95.0419769, -36.5625), (1223.5, -90.120285, -60.69078825), (1223.5, -38.8894758, 93.88750075), (1223.5, -61.86422350000001, 80.623024), (1223.5, -80.623024, 61.8642235), (1223.5, -93.88750075, 38.8894758), (1223.5, 38.8894758, 93.88750075), (1223.5, 61.86422350000001, 80.623024), (1223.5, 80.623024, 61.8642235), (1223.5, 93.88750075, 38.8894758), (1223.5, 85.248918525, -84.5723686), (1223.5, 80.377553925, -108.45394885), (1223.5, 75.5061874, -132.335529), (1223.5, 70.634819, -156.2171095), (1223.5, 65.7634525, -180.098686), (1223.5, 60.892086000000006, -203.9802625), (1223.5, 56.020719500000006, -227.8618395), (1223.5, 51.149353000000005, -251.74341600000002), (1223.5, 46.277984599999996, -275.625), (1223.5, 41.4066181, -299.5065765), (1223.5, 36.535251599999995, -323.388153), (1223.5, 31.663883675, -347.269745), (1223.5, 26.792517175, -371.1513215), (1223.5, 21.921150675, -395.03289800000005), (1223.5, 17.049784199999998, -418.9144745), (1223.5, 12.178417465, -442.79605100000003), (1223.5, 7.307050825000001, -466.67762749999997), (1223.5, 64.0625, 0.0), (1223.5, 38.4375, 0.0), (1223.5, 12.8125, 0.0), (1223.5, -12.8125, 0.0), (1223.5, -38.4375, 0.0), (1223.5, -64.0625, 0.0), (1223.5, -7.307050825, -466.6776275), (1223.5, -12.178417465, -442.79605100000003), (1223.5, -17.049784199999998, -418.9144745), (1223.5, -21.921150675, -395.03289800000005), (1223.5, -26.792517175, -371.1513215), (1223.5, -31.663883675, -347.269745), (1223.5, -36.535251599999995, -323.388153), (1223.5, -41.4066181, -299.5065765), (1223.5, -46.277984599999996, -275.625), (1223.5, -51.149353000000005, -251.74341600000002), (1223.5, -56.0207195, -227.8618395), (1223.5, -60.892086000000006, -203.9802625), (1223.5, -65.7634525, -180.09868600000001), (1223.5, -70.634819, -156.2171095), (1223.5, -75.5061874, -132.335529), (1223.5, -80.377553925, -108.45394885), (1223.5, -85.248918525, -84.5723686)])


"Von Mises Jam Straight"
e= np.array([0.06539385, 0.07119075, 0.07738275, 0.07686955000000001, 0.07526469999999999, 0.08571315, 0.10670775, 0.05450015, 0.0432357, 0.0425478, 0.04093, 0.042533749999999995, 0.0418722, 0.04252915, 0.04214295, 0.04231235, 0.04214455, 0.0420385, 0.0420496, 0.0416792, 0.04199005, 0.0414802, 0.0429669, 0.04358185, 0.049747349999999996, 0.054993150000000005, 0.0898794, 0.121699, 0.0969829, 0.08619175000000001, 0.0835436, 0.08084155, 0.07815659999999999, 0.0823637, 0.09829175000000001, 0.139935, 0.141038, 0.13901449999999999, 0.13895400000000002, 0.137111, 0.1365785, 0.1351735, 0.1348645, 0.1339035, 0.134104, 0.1333885, 0.134242, 0.1332615, 0.1346365, 0.1319665, 0.133083, 0.126809, 0.13232, 0.14360299999999998, 0.152804, 0.145417, 0.0629212, 0.06779550000000001, 0.0802894, 0.1048646, 0.1336, 0.178384])

f = np.array([(1073.5, -13.2644758, 100.75370025000001), (1073.5, 13.2644758, 100.75370025000001), (1073.5, -100.75370025000001, 13.2644758), (1073.5, -100.01399229999998, -12.1875), (1073.5, -89.6875, 0.0), (1073.5, -90.120285, -60.69078825), (1073.5, -95.0419769, -36.5625), (1073.5, -2.43568468, -490.559204), (1073.5, 2.43568468, -490.559204), (1073.5, 89.6875, 0.0), (1073.5, 100.0139923, -12.1875), (1073.5, 100.75370025000001, 13.2644758), (1073.5, 90.120285, -60.69078825), (1073.5, 95.04197690000001, -36.5625), (1073.5, -38.8894758, 93.88750075), (1073.5, -61.86422350000001, 80.623024), (1073.5, -80.623024, 61.8642235), (1073.5, -93.88750075, 38.8894758), (1073.5, -7.307051179999999, -466.6776275), (1073.5, -12.1784177, -442.79605100000003), (1073.5, -17.049784199999998, -418.9144745), (1073.5, -21.9211502, -395.03289800000005), (1073.5, -26.7925167, -371.1513215), (1073.5, -31.6638832, -347.269745), (1073.5, -36.535251599999995, -323.388153), (1073.5, -41.4066181, -299.5065765), (1073.5, -46.277984599999996, -275.625), (1073.5, -51.149353000000005, -251.74341600000002), (1073.5, -56.0207195, -227.8618395), (1073.5, -60.892086000000006, -203.9802625), (1073.5, -65.7634525, -180.09868600000001), (1073.5, -70.634819, -156.2171095), (1073.5, -75.5061874, -132.335529), (1073.5, -80.37755200000001, -108.45394885), (1073.5, -85.2489166, -84.5723686), (1073.5, 64.0625, 0.0), (1073.5, 38.4375, 0.0), (1073.5, 12.8125, 0.0), (1073.5, -12.8125, 0.0), (1073.5, -38.4375, 0.0), (1073.5, -64.0625, 0.0), (1073.5, 85.2489166, -84.5723686), (1073.5, 80.37755200000001, -108.45394885), (1073.5, 75.5061874, -132.335529), (1073.5, 70.634819, -156.2171095), (1073.5, 65.7634525, -180.098686), (1073.5, 60.892086000000006, -203.9802625), (1073.5, 56.020719500000006, -227.8618395), (1073.5, 51.149353000000005, -251.74341600000002), (1073.5, 46.277984599999996, -275.625), (1073.5, 41.4066181, -299.5065765), (1073.5, 36.535251599999995, -323.388153), (1073.5, 31.6638832, -347.269745), (1073.5, 26.7925167, -371.1513215), (1073.5, 21.9211502, -395.03289800000005), (1073.5, 17.049784199999998, -418.9144745), (1073.5, 12.1784177, -442.79605100000003), (1073.5, 7.30705118, -466.67762749999997), (1073.5, 38.8894758, 93.88750075), (1073.5, 61.86422350000001, 80.623024), (1073.5, 80.623024, 61.8642235), (1073.5, 93.88750075, 38.8894758)])

"S12 Bending"
g = np.array([-0.015370600000000002, -0.00477873, 0.00422869, 0.01932675, 0.045613, 0.07755395000000001, 0.057743699999999995, 0.06750805, 0.0503835, 0.0468539, 0.0393374, 0.035602499999999995, 0.031486, 0.02933375, 0.027210449999999997, 0.0263428, 0.02558935, 0.0257477, 0.02612105, 0.02736385, 0.02898165, 0.0317567, 0.035006300000000004, 0.040458, 0.0451005, 0.05637895, 0.05557395, 0.0120726, 0.009377475, 0.009013675, 0.0120644, 0.018055349999999998, 0.02451205, 0.0414831, 0.053885249999999996, -0.02462655, -0.02450015, -0.025235800000000003, -0.025427400000000003, -0.0261085, -0.02658345, -0.02726155, -0.0280546, -0.028932350000000003, -0.0303122, -0.031720700000000004, -0.03417665, -0.036587049999999996, -0.040898699999999996, -0.0446834, -0.052432, -0.055738449999999995, -0.0732546, -0.0627711, -0.023401150000000003, -0.02320725, -0.02470515, -0.0337603, -0.04017695, -0.040517449999999997, -0.031202149999999998, -0.0185325])


h = np.array([(1198.5, -13.2644758, 100.75370025000001), (1198.5, 13.2644758, 100.75370025000001), (1198.5, -100.75370025000001, 13.2644758), (1198.5, -100.01399229999998, -12.1875), (1198.5, -89.6875, 0.0), (1198.5, 89.6875, 0.0), (1198.5, 100.0139923, -12.1875), (1198.5, 100.75370025000001, 13.2644758), (1198.5, 90.120285, -60.69078825), (1198.5, 95.04197690000001, -36.5625), (1198.5, -2.4356840876332786, -490.559204), (1198.5, 2.4356840823667216, -490.559204), (1198.5, -90.120285, -60.69078825), (1198.5, -95.0419769, -36.5625), (1198.5, -38.8894758, 93.88750075), (1198.5, -61.86422350000001, 80.623024), (1198.5, -80.623024, 61.8642235), (1198.5, -93.88750075, 38.8894758), (1198.5, 93.88750075, 38.8894758), (1198.5, 80.623024, 61.8642235), (1198.5, 61.86422350000001, 80.623024), (1198.5, 38.8894758, 93.88750075), (1198.5, 7.307050824999999, -466.67762749999997), (1198.5, 12.178417465, -442.79605100000003), (1198.5, 17.049784199999998, -418.9144745), (1198.5, 21.921150675, -395.03289800000005), (1198.5, 26.792517175, -371.1513215), (1198.5, 31.663883675, -347.269745), (1198.5, 36.535251599999995, -323.388153), (1198.5, 41.4066181, -299.5065765), (1198.5, 46.277984599999996, -275.625), (1198.5, 51.149353000000005, -251.74341600000002), (1198.5, 56.020719500000006, -227.8618395), (1198.5, 60.892086000000006, -203.9802625), (1198.5, 65.7634525, -180.098686), (1198.5, 70.634819, -156.2171095), (1198.5, 75.5061874, -132.335529), (1198.5, 80.377553925, -108.45394885), (1198.5, 85.248918525, -84.5723686), (1198.5, -7.307050825, -466.6776275), (1198.5, -12.178417465, -442.79605100000003), (1198.5, -17.049784199999998, -418.9144745), (1198.5, -21.921150675, -395.03289800000005), (1198.5, -26.792517175, -371.1513215), (1198.5, -31.663883675, -347.269745), (1198.5, -36.535251599999995, -323.388153), (1198.5, -41.4066181, -299.5065765), (1198.5, -46.277984599999996, -275.625), (1198.5, -51.149353000000005, -251.74341600000002), (1198.5, -56.0207195, -227.8618395), (1198.5, -60.892086000000006, -203.9802625), (1198.5, -65.7634525, -180.09868600000001), (1198.5, -70.634819, -156.2171095), (1198.5, -75.5061874, -132.335529), (1198.5, -80.377553925, -108.45394885), (1198.5, -85.24891852500001, -84.5723686), (1198.5, 64.0625, 0.0), (1198.5, 38.4375, 0.0), (1198.5, 12.8125, 0.0), (1198.5, -12.8125, 0.0), (1198.5, -38.4375, 0.0), (1198.5, -64.0625, 0.0)])

"S12 Jam Bent"
i = np.array([-0.0261377, -0.04176935, -0.04988365, -0.058439649999999996, -0.07066205, -0.08341575000000001, -0.008253845, 0.0128026, 0.0343729, 0.0542435, 0.066194, 0.07533315, 0.0716843, 0.08522859999999999, 0.08845185, 0.09227725, 0.093831, 0.09486575, 0.09546615, 0.09599845, 0.09632665, 0.09693915, 0.0972645, 0.0983719, 0.09888345, 0.101183, 0.1020945, 0.10681099999999999, 0.1079295, 0.116976, 0.1148565, 0.1326395, 0.108152, 0.0294833, 0.0456245, 0.04819735, 0.0496933, 0.04254795, 0.0347268, 0.0145533, -0.010985249999999998, -0.0634196, -0.07194255, -0.0419843, -0.061001, -0.04419175, -0.0453312, -0.03648685, -0.03534215, -0.03068375, -0.0299264, -0.02775305, -0.0275773, -0.02676095, -0.0271674, -0.0271811, -0.028470950000000002, -0.029297049999999998, -0.032412250000000004, -0.034653500000000004, -0.042093950000000005, -0.048099])


j = np.array([(1223.5, -13.2644758, 100.75370025000001), (1223.5, 13.2644758, 100.75370025000001), (1223.5, -100.75370025000001, 13.2644758), (1223.5, -89.6875, 0.0), (1223.5, -100.01399229999998, -12.1875), (1223.5, 100.75370025000001, 13.2644758), (1223.5, 100.0139923, -12.1875), (1223.5, 89.6875, 0.0), (1223.5, 95.04197690000001, -36.5625), (1223.5, 90.120285, -60.69078825), (1223.5, 2.4356840823667216, -490.559204), (1223.5, -2.4356840876332786, -490.559204), (1223.5, -95.0419769, -36.5625), (1223.5, -90.120285, -60.69078825), (1223.5, -38.8894758, 93.88750075), (1223.5, -61.86422350000001, 80.623024), (1223.5, -80.623024, 61.8642235), (1223.5, -93.88750075, 38.8894758), (1223.5, 38.8894758, 93.88750075), (1223.5, 61.86422350000001, 80.623024), (1223.5, 80.623024, 61.8642235), (1223.5, 93.88750075, 38.8894758), (1223.5, 85.248918525, -84.5723686), (1223.5, 80.377553925, -108.45394885), (1223.5, 75.5061874, -132.335529), (1223.5, 70.634819, -156.2171095), (1223.5, 65.7634525, -180.098686), (1223.5, 60.892086000000006, -203.9802625), (1223.5, 56.020719500000006, -227.8618395), (1223.5, 51.149353000000005, -251.74341600000002), (1223.5, 46.277984599999996, -275.625), (1223.5, 41.4066181, -299.5065765), (1223.5, 36.535251599999995, -323.388153), (1223.5, 31.663883675, -347.269745), (1223.5, 26.792517175, -371.1513215), (1223.5, 21.921150675, -395.03289800000005), (1223.5, 17.049784199999998, -418.9144745), (1223.5, 12.178417465, -442.79605100000003), (1223.5, 7.307050825000001, -466.67762749999997), (1223.5, 64.0625, 0.0), (1223.5, 38.4375, 0.0), (1223.5, 12.8125, 0.0), (1223.5, -12.8125, 0.0), (1223.5, -38.4375, 0.0), (1223.5, -64.0625, 0.0), (1223.5, -7.307050825, -466.6776275), (1223.5, -12.178417465, -442.79605100000003), (1223.5, -17.049784199999998, -418.9144745), (1223.5, -21.921150675, -395.03289800000005), (1223.5, -26.792517175, -371.1513215), (1223.5, -31.663883675, -347.269745), (1223.5, -36.535251599999995, -323.388153), (1223.5, -41.4066181, -299.5065765), (1223.5, -46.277984599999996, -275.625), (1223.5, -51.149353000000005, -251.74341600000002), (1223.5, -56.0207195, -227.8618395), (1223.5, -60.892086000000006, -203.9802625), (1223.5, -65.7634525, -180.09868600000001), (1223.5, -70.634819, -156.2171095), (1223.5, -75.5061874, -132.335529), (1223.5, -80.377553925, -108.45394885), (1223.5, -85.248918525, -84.5723686)])

"S12 Jam Straight"
k = np.array([-0.032919000000000004, -0.0407639, -0.0446572, -0.042213, -0.03975835, -0.0268, 0.011308470000000001, -0.0010979524999999999, -0.00370194, -0.01014239, -0.01122045, -0.01424225, -0.01429165, -0.015524300000000001, -0.01544155, -0.01597495, -0.0160039, -0.016248699999999998, -0.01647725, -0.01661605, -0.017249300000000002, -0.0175304, -0.019185149999999998, -0.0202999, -0.024094499999999998, -0.026242349999999998, -0.0354155, 0.0304811, 0.042537649999999996, 0.0469542, 0.05032185, 0.0509129, 0.048735, 0.0459464, 0.034888050000000004, 0.07274215, 0.07744175, 0.07919435, 0.0812789, 0.0817262, 0.0826041, 0.08263465, 0.08306474999999999, 0.0828817, 0.0831562, 0.0826769, 0.0828821, 0.0817205, 0.08166225, 0.0787205, 0.07742505, 0.07083075, 0.06796745, 0.05195155, 0.0574656, 0.06899849999999999, -0.02097865, -0.005625514999999999, 0.01115375, 0.03032495, 0.0456373, 0.06223545])

l = np.array([(1198.5, -13.2644758, 100.75370025000001), (1198.5, 13.2644758, 100.75370025000001), (1198.5, -100.75370025000001, 13.2644758), (1198.5, -100.01399229999998, -12.1875), (1198.5, -89.6875, 0.0), (1198.5, 89.6875, 0.0), (1198.5, 100.0139923, -12.1875), (1198.5, 100.75370025000001, 13.2644758), (1198.5, 90.120285, -60.69078825), (1198.5, 95.04197690000001, -36.5625), (1198.5, -2.4356840876332786, -490.559204), (1198.5, 2.4356840823667216, -490.559204), (1198.5, -90.120285, -60.69078825), (1198.5, -95.0419769, -36.5625), (1198.5, -38.8894758, 93.88750075), (1198.5, -61.86422350000001, 80.623024), (1198.5, -80.623024, 61.8642235), (1198.5, -93.88750075, 38.8894758), (1198.5, 93.88750075, 38.8894758), (1198.5, 80.623024, 61.8642235), (1198.5, 61.86422350000001, 80.623024), (1198.5, 38.8894758, 93.88750075), (1198.5, 7.307050824999999, -466.67762749999997), (1198.5, 12.178417465, -442.79605100000003), (1198.5, 17.049784199999998, -418.9144745), (1198.5, 21.921150675, -395.03289800000005), (1198.5, 26.792517175, -371.1513215), (1198.5, 31.663883675, -347.269745), (1198.5, 36.535251599999995, -323.388153), (1198.5, 41.4066181, -299.5065765), (1198.5, 46.277984599999996, -275.625), (1198.5, 51.149353000000005, -251.74341600000002), (1198.5, 56.020719500000006, -227.8618395), (1198.5, 60.892086000000006, -203.9802625), (1198.5, 65.7634525, -180.098686), (1198.5, 70.634819, -156.2171095), (1198.5, 75.5061874, -132.335529), (1198.5, 80.377553925, -108.45394885), (1198.5, 85.248918525, -84.5723686), (1198.5, -7.307050825, -466.6776275), (1198.5, -12.178417465, -442.79605100000003), (1198.5, -17.049784199999998, -418.9144745), (1198.5, -21.921150675, -395.03289800000005), (1198.5, -26.792517175, -371.1513215), (1198.5, -31.663883675, -347.269745), (1198.5, -36.535251599999995, -323.388153), (1198.5, -41.4066181, -299.5065765), (1198.5, -46.277984599999996, -275.625), (1198.5, -51.149353000000005, -251.74341600000002), (1198.5, -56.0207195, -227.8618395), (1198.5, -60.892086000000006, -203.9802625), (1198.5, -65.7634525, -180.09868600000001), (1198.5, -70.634819, -156.2171095), (1198.5, -75.5061874, -132.335529), (1198.5, -80.377553925, -108.45394885), (1198.5, -85.24891852500001, -84.5723686), (1198.5, 64.0625, 0.0), (1198.5, 38.4375, 0.0), (1198.5, 12.8125, 0.0), (1198.5, -12.8125, 0.0), (1198.5, -38.4375, 0.0), (1198.5, -64.0625, 0.0)])

vmb=(0.3617365, 5190.0) #von mises bending
vmjb=(0.4141685, 2391.0) #von mises Jam Bent
vmjs=(0.178384, 5369.0) #von mises Jam Straight
s12b=(0.07755395000000001, 384.0) #s12 bending
s12jb=(0.1326395, 2390.0) #s12 Jam Bent
s12js=(0.0831562, 5181.0) #s12 Jam Straight



def plotcrossection(values,coordinates,plotname):

    ylist=[]
    zlist=[]
    for l in coordinates:
        ylist.append(l[1])
        zlist.append(l[2])

    plt.scatter(zlist,ylist,c=values)
    plt.title(plotname)
    plt.xlabel('z-location')
    plt.ylabel('y-location')
    plt.colorbar()
    plt.show()
    return

def plotcrossection(values,coordinates,plotname):

    ylist=[]
    zlist=[]
    for l in coordinates:
        ylist.append(l[1])
        zlist.append(l[2])

    plt.scatter(zlist,ylist,c=values)
    plt.title(plotname)
    plt.xlabel('z-location')
    plt.ylabel('y-location')
    plt.colorbar()
    plt.show()
    return


plotcrossection(a,b,'Von Mises Bending')
#plotcrossection(c,phi_1,'Von Mises Jam Bent' )
#plotcrossection(e,f,'Von Mises Jam Straight')
#plotcrossection(g,h, 'S12 Bending')
#plotcrossection(i,j,'S12 Jam Bent')
#plotcrossection(k,l,'S12 Jam Straight')




















