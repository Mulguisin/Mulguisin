from time import time
import numpy as np
#import get_density
#import mgs_code

class mulguisin:
	def __init__(self, Rcut, x1, y1, z1):
		self.Rcut = Rcut
		self.x1 = x1
		self.y1 = y1
		self.z1 = z1
		self.Nmgs = None
		self.imgs = None
		self.cng = None
		self.clm = None
		self.clg = None
		self.isort = None

	def get_isort(self):
		positions = np.vstack((self.x1,self.y1,self.z1)).T
		print('Calculate Voronoi density')
		sta = time()
		den = get_density.voronoi_density(positions)
		end = time()
		print('Calculation is done. Time = ', end - sta)
		isort = np.flip(den.argsort())
		return isort

	def get_mgs(self):
		isort = self.get_isort()
		print('Calculate MGS')
		sta = time()
		self.Nmgs,self.imgs,self.clg,self.clm,self.cng = mgs_code.mgs(self.Rcut, self.x1, self.y1, self.z1, isort)
		end = time()
		print('Calculation is done. Time = ', end - sta)
		return self.Nmgs, self.imgs, self.clg, self.clm, self.cng

	def Get_TotLength(self,imgs=None):
		lsum = 0.
		for i in range(self.cng[imgs]):
			id_mom = self.clm[imgs][i]
			id_chi = self.clg[imgs][i]
			dist = (self.x1[id_mom]-self.x1[id_chi])**2 + (self.y1[id_mom]-self.y1[id_chi])**2 + (self.z1[id_mom]-self.z1[id_chi])**2
			lsum += np.sqrt(dist)
		return lsum

	def Get_Node(self,imgs=None,idgal=None):
		for i in range(self.cng[imgs]):
			if self.clg[imgs][i] == idgal:
				return self.clm[imgs][i]
	
	def Get_generation(self,imgs=None,idgal=None):
		ngen = 0
		id_chi = self.clg[imgs][idgal]
		id_mom = self.clg[imgs][0] # The origin of the generation, not mother galaxy
		while (id_chi != id_mom):
			id_node = self.Get_Node(imgs,id_chi)
			ngen+=1
			id_chi = id_node
		return ngen
	
	def Get_Length(self,imgs=None,idgal=None):
		lsum = 0.
		id_chi = self.clg[imgs][idgal]
		id_mom = self.clg[imgs][0]
		while (id_chi != id_mom):
			id_node = self.Get_Node(imgs,id_chi)
			dist = (self.x1[id_chi]-self.x1[id_node])**2 + (self.y1[id_chi]-self.y1[id_node])**2 + (self.z1[id_chi]-self.z1[id_node])**2
			lsum += np.sqrt(dist)
			id_chi = id_node
		return lsum
	
	def Get_Child(self,imgs=None,idgal=None):
		nchild = 0
		for i in range(idgal+1, self.cng[imgs]):
			if self.clm[imgs][i]==self.clg[imgs][idgal]: nchild+=1
		return nchild
	
	def Get_Degree(self,imgs=None,idgal=None):
		degree = 0
		id_mom = self.clg[imgs][0]
		for i in range(idgal+1,self.cng[imgs]):
			if self.clm[imgs][i] == self.clg[imgs][idgal]: degree+=1
		if self.clg[imgs][idgal] != id_mom: degree += 1
		return degree
	
	def Get_Branch(self,imgs=None):
		nbranch = 0
		for i in range(self.cng[imgs]):
			if self.Get_Child(imgs,i) == 0: nbranch += 1
		return nbranch
	
	def Get_OpenAngle(self,imgs=None,idgal=None):
		ang = -1
		id_chi = self.clg[imgs][idgal]
		id_mom = self.clg[imgs][0]
		if id_chi!=id_mom:
			id_node = self.Get_Node(imgs,id_chi)
			if id_node!=id_mom:
				id_axis = self.Get_Node(imgs,id_node)
				ang = self.Get_Angle(self.x1[id_chi],self.y1[id_chi],self.z1[id_chi],self.x1[id_node],self.y1[id_node],self.z1[id_node],self.x1[id_axis],self.y1[id_axis],self.z1[id_axis])
		return ang
	
	def Get_LinkPloarAngle(self,imgs=None,idgal=None):
		ang = -1
		id_chi = self.clg[imgs][idgal]
		id_mom = self.clg[imgs][0]
		if id_chi!=id_mom:
			id_node = self.Get_Node(imgs,id_chi)
			ang = self.Get_Angle(self.x1[id_chi],self.y1[id_chi],self.z1[id_chi],self.x1[id_node],self.y1[id_node],self.z1[id_node],0,0,0)
		return ang

	def Get_PolarAngle(self,imgs=None,idgal=None):
		ang = -1
		id_chi = self.clg[imgs][idgal]
		id_mom = self.clg[imgs][0]
		if id_chi!=id_mom:
			ang = self.Get_Angle(self.x1[id_chi],self.y1[id_chi],self.z1[id_chi],self.x1[id_mom],self.y1[id_mom],self.z1[id_mom],0,0,0)
		return ang
	
	def Get_Angle(self,x3=None,y3=None,z3=None,x2=None,y2=None,z2=None,x1=None,y1=None,z1=None):
		da = np.sqrt((x3-x2)**2 + (y3-y2)**2 + (z3-z2)**2)
		db = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
		dc = (x3-x1)**2 + (y3-y1)**2 + (z3-z1)**2
		cosalpha = (da**2 + db**2 - dc)/(2.*da*db)
		return np.arccos(cosalpha)
