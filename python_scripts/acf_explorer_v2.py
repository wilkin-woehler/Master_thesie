#author : wilkin woehler
#date : 11.05.2020

import h5py
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from scipy.stats import norm

class datenset:
	def __init__(self,path,N,rho_i,rho_f,eq,pr,seed):
		self.filename = path+"Out_{0}_{1}_{2}_eq{3}k_pr{4}k_{5}.h5".format(N,rho_i,rho_f,eq,pr,seed) 
		self.read_file = h5py.File(self.filename, "r")
		self.wdata = self.read_file['observables/watchtable']
	
	def acf(self):
		#extract time and largest cluster from watchtable
		start_offset=100
		self.largest_cluster=np.zeros(len(self.wdata)-start_offset)
		self.time=np.zeros(len(self.wdata)-start_offset)
		for i in range(start_offset,len(self.wdata)):
			self.largest_cluster[i-start_offset] = self.wdata[i][3]
			self.time[i-start_offset] = self.wdata[i][1]
		self.mean_lc = np.mean(self.largest_cluster)
		self.fluc_lc = self.largest_cluster - self.mean_lc
		N_acf=len(self.largest_cluster)
		self.tau = np.zeros(N_acf)
		
		for i in range(1,N_acf):
			delta_time=self.time[i:]-self.time[:-i]
			self.tau[i]=np.mean(delta_time)

		self.acf=np.zeros(N_acf)
		self.acf[0]=np.mean(self.fluc_lc**2)
		for i in range(1,N_acf):
			self.acf[i]=np.mean(self.fluc_lc[i:]*self.fluc_lc[:-i])


outfile_path="../outfiles/"
datensets=[]

rho_final=np.linspace(526,530,3,dtype=int)
seeds=np.linspace(10,12,3,dtype=int)

for rho_f in rho_final:
	for seed in seeds:
		datensets.append(datenset(outfile_path,1048576,450,rho_f,1,1,seed))
		#datensets[-1].acf()
		#print(rho_f,seed)

overview=False
if overview:
	fig,axe = plt.subplots(3,3)
	for i in range(3):
		for j in range(3):
			k=i*3+j
			dset=datensets[k]
			dset.acf()
			print("Calculated acf of {0} {1}".format(i,j))
			axe[i,0].scatter(dset.time, dset.largest_cluster,c='C{0}'.format(j),alpha=0.3)
			axe[i,0].plot(dset.time, dset.mean_lc*np.ones(len(dset.time)) ,c='C{0}'.format(j))
			axe[i,1].scatter(dset.time,dset.fluc_lc,c='C{0}'.format(j))
			axe[i,2].scatter(dset.tau,dset.acf,c='C{0}'.format(j))
	plt.show()
else:
	fig,axe = plt.subplots(3,sharex=True,sharey=True)
	ylim=0
	xlim=60
	for i in range(3):
		for j in range(3):
			k=i*3+j
			dset=datensets[k]
			dset.acf()
			print("Calculated acf of {0} {1}".format(rho_final[i],seeds[j]))
			axe[i].scatter(dset.tau,dset.acf,c='C{0}'.format(j), s=10)
			if max(dset.acf)>ylim:
				ylim = max(dset.acf)
		axe[i].text(0.5*xlim,70,"ACF for eta = {0:4.3f}".format(rho_final[i]/1000.0))
		axe[i].set_ylabel(r"ACF($\tau$)")
		axe[i].axhline(0, c='k', linewidth=1)
	axe[0].set_ylim(-ylim*0.1,ylim)
	axe[0].set_xlim(-5,xlim)
	axe[2].set_xlabel(r'$\tau$ ')
	plt.savefig("ACF.pdf")

