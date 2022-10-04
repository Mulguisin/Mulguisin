import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import time
import collections
import mgs_class

# -------- Data read --------------
f = h5.File("rand_gal_data_210218_sig10.hdf5",'r')

x1 = np.array(f['x_rsd'])
y1 = np.array(f['y_rsd'])
z1 = np.array(f['z_rsd'])
f.close()
print('Length of data : ',len(x1))

# -----------------------------------

# ------- Calculate MGS -----------------
Rcut = 10.0

MGS = mgs_class.mulguisin(Rcut,x1,y1,z1)

Nmgs, imgs, clg, clm, cng = MGS.get_mgs()

print('Total number of mgs : ' ,Nmgs)

# Make link information
"""
links = np.zeros((1,3))
for i in range(Nmgs):
	for j in range(cng[i]):
		id_mom = clm[i][j]
		id_chi = clg[i][j]
		links = np.append(links,[[id_mom,id_chi,i]],axis=0)
links = np.delete(links,0,0)
"""

# Make topological information
"""
sta = time()
Totlength = []
for i in range(Nmgs):
    Totlength.append(MGS.Get_TotLength(i))

Generations = []
for i in range(Nmgs):
	for j in range(cng[i]):
		Generations.append(MGS.Get_generation(i,j))

Length = []
for i in range(Nmgs):
	for j in range(cng[i]):
		Length.append(MGS.Get_Length(i,j))

nChild = []
for i in range(Nmgs):
	for j in range(cng[i]):
		nChild.append(MGS.Get_Child(i,j))

nDegree = []
for i in range(Nmgs):
	for j in range(cng[i]):
		nDegree.append(MGS.Get_Degree(i,j))

nBranch = []
for i in range(Nmgs):
	nBranch.append(MGS.Get_Branch(i))

OpenAngle = []
for i in range(Nmgs):
	for j in range(cng[i]):
		openangle = MGS.Get_OpenAngle(i,j)
		OpenAngle.append(openangle)

PolarAngle = []
for i in range(Nmgs):
	for j in range(cng[i]):
		PolarAngle.append(MGS.Get_LinkPloarAngle(i,j))

PolarAngle_origin = []
for i in range(Nmgs):
	for j in range(cng[i]):
		PolarAngle_origin.append(MGS.Get_PolarAngle(i,j))

end = time()
print('Time for calculating topological information = ', end - sta)

Totlength = np.asarray(Totlength)
Length = np.asarray(Length)
Generations = np.asarray(Generations)
nChild = np.asarray(nChild)
nDegree = np.asarray(nDegree)
nBranch = np.asarray(nBranch)
OpenAngle = np.asarray(OpenAngle)
PolarAngle = np.asarray(PolarAngle)
PolarAngle_origin = np.asarray(PolarAngle_origin)

mask_open = np.where(OpenAngle>=0)
mask_polar = np.where(PolarAngle>=0)
mask_polar_origin = np.where(PolarAngle_origin>=0)

fig, ax = plt.subplots(2,4,figsize=(16,8))
ax[0,0].hist(Totlength,histtype='step')
ax[0,0].set_xlabel('Total length of each cluster')
ax[0,1].hist(Length,bins=100,histtype='step')
ax[0,1].set_xlabel('Total length of one generation')
ax[0,2].hist(Generations,histtype='step')
ax[0,2].set_xlabel('Generation')
ax[0,3].hist(nChild,bins=6,histtype='step')
ax[0,3].set_xlabel('Number of children for each galaxy')
ax[1,0].hist(nDegree,histtype='step')
ax[1,0].set_xlabel('Degree')
ax[1,1].hist(nBranch,bins=35,histtype='step')
ax[1,1].set_xlabel('Number of branch')
ax[1,2].hist(OpenAngle[mask_open]*(180./np.pi),bins=20,histtype='step')
ax[1,2].set_xlabel('Opening Angle')
ax[1,3].hist(PolarAngle[mask_polar]*(180./np.pi),bins=30,histtype='step')
ax[1,3].set_xlabel('Polar Angle')

plt.tight_layout()
plt.show()
"""
