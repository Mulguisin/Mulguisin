from numba import njit
from numba.typed import List

@njit
def mgs(Rcut,x1,y1,z1,isort):
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
