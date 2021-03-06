from validation import *

#a = get_dat('bending', 'stresses')
x , y = get_twist('Jam_Bent')
x2 , y2 = get_twist('Jam_Straight')
x3 , y3 = get_twist('bending')

plt.scatter(x, y)
plt.scatter(x2, y2, label = 'Jam Straight')
plt.scatter(x3, y3, label = 'Bending')
plt.scatter(x,(y2+y3), label = 'Bending + Jam Straight')
plt.title('Rate of twist ')
plt.xlabel('x - location')
plt.ylabel('Angle [rad]')
plt.legend()

plt.show()

