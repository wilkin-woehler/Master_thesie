#author : wilkin woehler
#date : 11.05.2020

import h5py
import numpy as np
import matplotlib.pyplot as plt
from datetime import date


#read in h5(md) file from which data is supposed to be read
N=13500
rho=0.53

filename = "../Out_{0}_{1}.h5".format(N,int(rho*1000))
#filename = "../Out_0.75_1.00_sm.h5"
read_file=h5py.File(filename, "r")

#show all folders and subfolders, from which a specific one can be picked
def print_names(name):
    print(name)
read_file.visit(print_names)

#show all folders and subfolders, from which a specific one can be picked
def print_names(name,object_name):
    print(name)
    print(object_name)
read_file.visititems(print_names)

print(read_file['parameters/all_parts'])


V_spheres=N*4/3*np.pi*1/2**3
V_box=V_spheres/rho
box_len=V_box**(1/3)
print("box length is : ", box_len)

data=read_file['particles/pr_all_parts/position/value'][()]
step_data=read_file['particles/pr_all_parts/position/step'][()]

tmp=read_file['particles/pr_all_parts/position/value'][()]
part_num=len(tmp[0])
fdate = date.today().strftime('%d_%m_%Y')

msd=np.zeros(len(data))
s_msd=np.zeros(len(data))
steps=np.zeros(len(data))
for i in range(len(data)): 
	r=data[i]-data[0]
	r_pbc= r- np.round( r/box_len)*box_len 
	r2= np.sum((r_pbc)**2,axis=1 )
	print(r2)
	msd[i] = np.mean(r2)
	s_msd[i] = np.std(r2)/np.sqrt(part_num)
	steps[i] = step_data[i]



plt.figure()
#plt.errorbar(steps/part_num,msd,s_msd)
plt.plot(steps/part_num,msd)
plt.ylabel('MSD')
plt.xlabel('steps/particle')
plt.savefig('Diffusion_{0}_{1}_{2}.pdf'.format(fdate,part_num,int(rho*1000)))


print(read_file['parameters/all_parts'])
print("-------------------")


