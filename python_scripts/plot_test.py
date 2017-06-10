import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math

def potential_ecc_dlat(lat,lon,radius,omega,e,time):
    # Convert lat to colat
    # lat = np.pi/2.0 - lat

    # If numpy array, convert lat and lon using meshgrid
    # lat, lon = np.meshgrid(lat,lon)

    factor = omega**2 * radius**2*e


    M = omega*time
    pot = -factor *0.75*np.sin(2*lat) * (3*np.cos(M) * (1.0 + np.cos(2*lon)) + 4*np.sin(M)*np.sin(2*lon))


    return pot

# Creating a Triangulation without specifying the triangles results in the
# Delaunay triangulation of the points.

x = np.loadtxt("lons.txt",delimiter=',')
y = np.loadtxt("lats.txt",delimiter=',')
z = np.loadtxt("z.txt",delimiter=',')

# Create the Triangulation; no triangles so Delaunay triangulation created.
triang = tri.Triangulation(x, y)

fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2)

data = potential_ecc_dlat(y,x,2575.5e3,4.559e-6,0.0288, 0.0)


ax1.set_aspect('equal')
ax2.set_aspect('equal')
c1 = ax1.tricontourf(triang, z, cmap=plt.cm.Spectral)
plt.colorbar(c1,ax=ax1,orientation='horizontal')

c2 = ax2.tricontourf(triang, data, cmap=plt.cm.Spectral)
plt.colorbar(c2,ax=ax2,orientation='horizontal')

ax1.set_xlim([0,np.pi*2.0])
ax1.set_ylim([-0.5*np.pi,0.5*np.pi])

ax2.set_xlim([0,np.pi*2.0])
ax2.set_ylim([-0.5*np.pi,0.5*np.pi])
ax = plt.gca()
# ax.set_xlim

plt.show()
