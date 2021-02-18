#author : wilkin woehler
#date : 11.05.2020

import h5py
import numpy as np
import matplotlib.pyplot as plt
from datetime import date


#read in h5(md) file from which data is supposed to be read
N=13500
rho=0.53
names=[]
seeds=np.linspace(1,4,4,dtype=int)


filename = 'outfiles/Out_32000_450_550_eq0k_pr1k_1.h5'
read_file=h5py.File(filename, "r")
l1=len(read_file['observables/Cluster_Distribution/value'])

maxc=200
pnt = np.zeros(shape=(maxc,l1,len(seeds)))
filenames=[]

for s in seeds:
	filenames.append('outfiles/Out_32000_450_550_eq0k_pr1k_{0}.h5'.format(s))

	read_file=h5py.File(filenames[s-1], "r")
#show all folders and subfolders, from which a specific one can be picked
#def print_names(name):
#    print(name)
#read_file.visit(print_names)

#show all folders and subfolders, from which a specific one can be picked
	def print_names(name,object_name):
	    print(name)
	    print(object_name)
	read_file.visititems(print_names)

#for i in range(3):
#	print(read_file['reset_sim/FEL_array/value'][i])
#print(read_file['particles/pr_big_clus/velocity/value'][0])

#for i in range(len(read_file['reset_sim/particles/velocity/step'])):
#	print(read_file['reset_sim/particles/velocity/step'][i])




#print(read_file['reset_sim/cells_first_part_id/value'][0])
#plt.hist(read_file['reset_sim/cells_first_part_id/value'][0])
#
#plt.savefig('test.png')

	for i in range(l1):
		print(read_file['observables/Cluster_Distribution/value'][i])
		for j in range(len(read_file['observables/Cluster_Distribution/value'][i][0])):
			try:
				k=read_file['observables/Cluster_Distribution/value'][i][0][j][0]
				if k < maxc:
					pnt[k,i,s]+=read_file['observables/Cluster_Distribution/value'][i][0][j][1]
			except:
				pass;

#for i in range(len(read_file['observables/watchtable'])):
#	print(read_file['observables/watchtable'][i])

pnt[0,:,:]=0
pnt=np.sum(pnt, axis=2)

pnt = np.log(pnt)
pnt[pnt==pnt[0,0]]=-1

fig = plt.figure()
im = plt.pcolormesh(np.linspace(-0.5,maxc+0.5,maxc+1),np.linspace(-0.5,l1+0.5,l1+1),pnt.T)
fig.colorbar(im)
plt.ylabel('t')
plt.xlabel('N')
plt.title('p(n,t) RAW')
#plt.xticks(phi,np.array(phi,dtype='str'))
#plt.yticks(N_part)
plt.show()
#plt.savefig("pnt.pdf")

print("-------------------")


