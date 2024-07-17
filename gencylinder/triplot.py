#!/usr/bin/python

# This function plots the the triangulated surface obaine from the simulations
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from matplotlib import cm

# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import




#print "This is the name of the script: ", sys.argv[0]
#print "Number of arguments: ", len(sys.argv)
#print "The arguments are: " , str(sys.argv)

MESH_PATH = "./outfiles/mesh.dat"

r0 = np.loadtxt("./outfiles/init_coords.dat", delimiter="\t") # Coordinates 1
tri = np.loadtxt(MESH_PATH, delimiter="\t") # Triangles

x0=r0[:,0]
y0=r0[:,1]
z0=r0[:,2]



# Calculates the coordination number
N_tri = len(tri)
Np = len(r0)
cn = np.zeros(Np)
for t in range(0, N_tri):
    v1 = int(tri[t,0])
    v2 = int(tri[t,1])
    v3 = int(tri[t,2])
    cn[v1] = cn[v1]+1
    cn[v2] = cn[v2]+1
    cn[v3] = cn[v3]+1
    
# Five-fold and seven fold defects
n5=0
n7=0
for i in range(0,Np):
    if cn[i] == 5:
        n5=n5+1
    if cn[i] == 7:
        n7=n7+1

print "Five-fold defects: " + str(n5)
print "Seven-fold defects: " + str(n7)



fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.set_xlim(-1.2,1.2)
#ax.set_ylim(-1.2,1.2)
#ax.set_zlim(-1.2,1.2)

#cmap=cm.coolwarm
#plt.cm.Spectral
#ax.plot_trisurf(x1, y1, z1, triangles=tri1, cmap=cm.coolwarm)
#edgecolor=(0.5,0.5,0.5,0.7)
ax.plot_trisurf(x0, y0, z0, triangles=tri, color=(1,0.0,0.0,1.0), edgecolor=(0.0,0.0,0.0,0.7))

plt.axis('off')
plt.show()
plt.savefig("mesh.jpg", dpi=600)
plt.close()
