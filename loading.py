import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from geometry import Dataset

def makeGrid(minX, minZ, maxX, maxZ, numx, numz):
    stepX = (maxX - minX) / (numx-1)
    stepZ = (maxZ - minZ) / (numz-1)
    return np.mgrid[minX:maxX+stepX*0.1:stepX, minZ:maxZ+stepZ*0.1:stepZ]


if __name__ == "__main__":
    data = np.genfromtxt("aerodata.csv", delimiter=",")
    data[0,0] = 0.034398
    size = data.size
    a = Dataset()

    x, z = makeGrid(a.minx, a.minz, a.maxx, a.maxz, 81, 41)
    fig = plt.figure()
    frame1 = fig.add_subplot(2,1,1,projection='3d')
    frame1.plot_surface(x, z, data, cmap=cm.coolwarm, rstride=1, cstride=1, linewidth=1, antialiased=True)
    xs = np.reshape(x, size)
    zs = np.reshape(z, size)
    ys = np.reshape(data, size)
    frame2 = fig.add_subplot(2,1,2,projection='3d')
    frame2.plot_trisurf(xs, zs, ys, cmap=cm.coolwarm, linewidth=1, antialiased=True)
    plt.show()



