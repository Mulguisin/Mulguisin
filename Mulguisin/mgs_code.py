from numba import njit
from numba.typed import List
from scipy.spatial import cKDTree
import numpy as np

@njit
def mgs3d(Rcut,x1,y1,z1,isort):
	Nmgs=0
	imgs = List()
	cng = List()
	clg = List()
	clm = List()
	#cll = List()
	for i in range(len(x1)):
		ic = isort[i]
		Icc=-1
		Iwg=-1
		Dmin=999

		for ig in range(0,Nmgs):
			Dist=999
			for icl in range(0,cng[ig]):
				itmp=clg[ig][icl]
				d2 = (x1[ic]-x1[itmp])**2+(y1[ic]-y1[itmp])**2+(z1[ic]-z1[itmp])**2
				if (d2<Dist): 
					idmin=itmp
					Dist=d2
			if (Dist < Dmin):
				Iwg=ig
				Icc=idmin
				Dmin=Dist
		if (Dmin < Rcut*Rcut):
			imgs.append(Iwg)
			clg[Iwg].append(ic)
			clm[Iwg].append(Icc)
			cng[Iwg]+=1
		else:
			imgs.append(Nmgs)	
			clg.append([ic])
			clm.append([ic])
			cng.append(1)
			Nmgs+=1
	return Nmgs,imgs,clg,clm,cng

@njit
def mgs2d(Rcut,x1,y1,isort):
	Nmgs=0
	imgs = List()
	cng = List()
	clg = List()
	clm = List()
	#cll = List()
	for i in range(len(x1)):
		ic = isort[i]
		Icc=-1
		Iwg=-1
		Dmin=999

		for ig in range(0,Nmgs):
			Dist=999
			for icl in range(0,cng[ig]):
				itmp=clg[ig][icl]
				d2 = (x1[ic]-x1[itmp])**2+(y1[ic]-y1[itmp])**2
				if (d2<Dist): 
					idmin=itmp
					Dist=d2
			if (Dist < Dmin):
				Iwg=ig
				Icc=idmin
				Dmin=Dist
		if (Dmin < Rcut*Rcut):
			imgs.append(Iwg)
			clg[Iwg].append(ic)
			clm[Iwg].append(Icc)
			cng[Iwg]+=1
		else:
			imgs.append(Nmgs)	
			clg.append([ic])
			clm.append([ic])
			cng.append(1)
			Nmgs+=1
	return Nmgs,imgs,clg,clm,cng

def mgs3d_tree(Rcut, x1, y1, z1, isort):
    """KD-tree based version of :func:`mgs3d`."""
    pos = np.vstack((x1, y1, z1)).T
    tree = cKDTree(pos)
    rank = np.empty(len(isort), dtype=np.int64)
    for i, idx in enumerate(isort):
        rank[idx] = i
    Nmgs = 0
    imgs = []
    cng = []
    clg = []
    clm = []
    group = np.full(len(x1), -1, dtype=np.int64)
    for i, ic in enumerate(isort):
        neigh = tree.query_ball_point(pos[ic], r=Rcut)
        best = -1
        dmin = None
        for nb in neigh:
            if nb == ic or rank[nb] >= i:
                continue
            d2 = ((x1[ic]-x1[nb])**2 + (y1[ic]-y1[nb])**2 + (z1[ic]-z1[nb])**2)
            if dmin is None or d2 < dmin:
                dmin = d2
                best = nb
        if best != -1:
            ig = group[best]
            imgs.append(ig)
            clg[ig].append(ic)
            clm[ig].append(best)
            cng[ig] += 1
            group[ic] = ig
        else:
            imgs.append(Nmgs)
            clg.append([ic])
            clm.append([ic])
            cng.append(1)
            group[ic] = Nmgs
            Nmgs += 1
    return Nmgs, imgs, clg, clm, cng


def mgs2d_tree(Rcut, x1, y1, isort):
    """KD-tree based version of :func:`mgs2d`."""
    pos = np.vstack((x1, y1)).T
    tree = cKDTree(pos)
    rank = np.empty(len(isort), dtype=np.int64)
    for i, idx in enumerate(isort):
        rank[idx] = i
    Nmgs = 0
    imgs = []
    cng = []
    clg = []
    clm = []
    group = np.full(len(x1), -1, dtype=np.int64)
    for i, ic in enumerate(isort):
        neigh = tree.query_ball_point(pos[ic], r=Rcut)
        best = -1
        dmin = None
        for nb in neigh:
            if nb == ic or rank[nb] >= i:
                continue
            d2 = ((x1[ic]-x1[nb])**2 + (y1[ic]-y1[nb])**2)
            if dmin is None or d2 < dmin:
                dmin = d2
                best = nb
        if best != -1:
            ig = group[best]
            imgs.append(ig)
            clg[ig].append(ic)
            clm[ig].append(best)
            cng[ig] += 1
            group[ic] = ig
        else:
            imgs.append(Nmgs)
            clg.append([ic])
            clm.append([ic])
            cng.append(1)
            group[ic] = Nmgs
            Nmgs += 1
    return Nmgs, imgs, clg, clm, cng

