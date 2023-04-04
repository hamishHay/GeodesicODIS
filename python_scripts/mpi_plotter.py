import h5py 
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys

NP = 2
GL = 4
ITER = int(sys.argv[1])

f1 = h5py.File("input_files/grid_l"+str(GL)+".000.h5", 'r')

node_dset = f1["NODES"]

NODE_NUM1 = node_dset.attrs["NODE_NUM_NO_GHOSTS"]

lat1 = node_dset["LAT"][:NODE_NUM1]
lon1 = node_dset["LON"][:NODE_NUM1]

f1.close()

f2 = h5py.File("input_files/grid_l"+str(GL)+".001.h5", 'r')

node_dset = f2["NODES"]

NODE_NUM2 = node_dset.attrs["NODE_NUM_NO_GHOSTS"]

lat2 = node_dset["LAT"][:NODE_NUM2]
lon2 = node_dset["LON"][:NODE_NUM2]

f2.close()

df = h5py.File("DATA/data.000.h5", 'r')
eta1 = df["displacement"][ITER, :NODE_NUM1]
df.close()

df = h5py.File("DATA/data.001.h5", 'r')
eta2 = df["displacement"][ITER, :NODE_NUM2]
df.close()

lon = np.concatenate( (lon1, lon2) )
lat = np.concatenate( (lat1, lat2) )
eta = np.concatenate( (eta1, eta2) )

triang_eta = tri.Triangulation(lon, lat)

# fig, ax = plt.subplots()

# ax.tricontourf(triang, eta)

# plt.show()

f1 = h5py.File("input_files/grid_l"+str(GL)+".000.h5", 'r')

face_dset = f1["FACES"]

FACE_NUM1 = face_dset.attrs["FACE_NUM_NO_GHOSTS"]

lat1 = face_dset["LAT"][:FACE_NUM1]
lon1 = face_dset["LON"][:FACE_NUM1]

f1.close()

f2 = h5py.File("input_files/grid_l"+str(GL)+".001.h5", 'r')

face_dset = f2["FACES"]

FACE_NUM2 = face_dset.attrs["FACE_NUM_NO_GHOSTS"]

lat2 = face_dset["LAT"][:FACE_NUM2]
lon2 = face_dset["LON"][:FACE_NUM2]

f2.close()

df = h5py.File("DATA/data.000.h5", 'r')
u1 = df["east velocity"][ITER, :FACE_NUM1]
v1 = df["north velocity"][ITER, :FACE_NUM1]
vx1 = df["x velocity"][ITER, :NODE_NUM1]
vy1 = df["y velocity"][ITER, :NODE_NUM1]
vz1 = df["z velocity"][ITER, :NODE_NUM1]
df.close()

df = h5py.File("DATA/data.001.h5", 'r')
u2 = df["east velocity"][ITER, :FACE_NUM2]
v2 = df["north velocity"][ITER, :FACE_NUM2]
vx2 = df["x velocity"][ITER, :NODE_NUM2]
vy2 = df["y velocity"][ITER, :NODE_NUM2]
vz2 = df["z velocity"][ITER, :NODE_NUM2]
df.close()

lon = np.concatenate( (lon1, lon2) )
lat = np.concatenate( (lat1, lat2) )
u = np.concatenate( (u1, u2) )
v = np.concatenate( (v1, v2) )
vx = np.concatenate( (vx1, vx2) )
vy = np.concatenate( (vy1, vy2) )
vz = np.concatenate( (vz1, vz2) )

triang_uv = tri.Triangulation(lon, lat)

fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(ncols=6)

ax1.tricontourf(triang_uv, u)
ax2.tricontourf(triang_uv, v)
ax3.tricontourf(triang_eta, eta)
ax4.tricontourf(triang_eta, vx)
ax5.tricontourf(triang_eta, vy)
ax6.tricontourf(triang_eta, vz)

# a1.set_aspect("equal")

plt.show()
