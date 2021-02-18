#author : wilkin woehler
#date : 11.05.2020
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import h5py

#read in h5(md) file from which data is supposed to be read
N=4000
rho=0.53 #rho_final in production phase

filename = "../Out_{0}_{1}.h5".format(N,int(rho*1000))
print('read in file: '+filename)
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


#open a text file to which the data can be written
#dump_file1=open("dump_position.txt","w+")
#dump_file2=open("dump_velocity.txt","w+")

#take out position data and write it in the dump_file
#data1=read_file['particles/eq_all_parts/position/value'][()]
#data2=read_file['particles/eq_all_parts/velocity/value'][()]

#for i in range(len(data1[0])):
#	dump_file1.write("{0} {1} {2}\n".format(data1[0][i][0],data1[0][i][1],data1[0][i][2]))
#	dump_file2.write("{0} {1} {2}\n".format(data2[0][i][0],data2[0][i][1],data2[0][i][2]))
	#print("Particle {0}".format(i))
	#print("p: {0} {1} {2}".format(data1[0][i][0],data1[0][i][1],data1[0][i][2]))
	#print("v: {0} {1} {2}".format(data2[0][i][0],data2[0][i][1],data2[0][i][2]))
#dump_file1.close()
#dump_file2.close()

print(read_file['parameters/all_parts'])
print("-------------------")

#obtain step numbers particle numbers and date
tmp=read_file['particles/eq_all_parts/position/step'][()]
eq_steps=tmp[-1]-tmp[0]

tmp=read_file['particles/co_all_parts/position/step'][()]
try:
	co_steps=tmp[-1]-tmp[0]
except:
	co_steps=0

tmp=read_file['particles/pr_all_parts/position/step'][()]
pr_steps=tmp[-1]-tmp[0]


tmp=read_file['particles/eq_all_parts/position/value'][()]
part_num=len(tmp[0])


V_spheres=N*4/3*np.pi*1/2**3
V_box=V_spheres/rho
box_len=V_box**(1/3)
print("box length is : ", box_len)


fdate = date.today().strftime('%d_%m_%Y')

#read in watchdata from measurements
data=read_file['observables/watchtable'][()]
steps=np.zeros(len(data))
solids=np.zeros(len(data))
msd=np.zeros(len(data))
s_msd=np.zeros(len(data))
for i in range(len(data)):
	steps[i]=data[i][0]
	solids[i]=data[i][1]
	msd[i]=data[i][3]
	s_msd[i]=data[i][4]
print(s_msd)

fig,axe = plt.subplots(2, sharex=True)
axe[0].plot(steps/part_num,solids)
axe[0].axvline(0)
axe[0].text(0,part_num*0.8,'{0:4.2e} Equilibration steps'.format(eq_steps))
axe[0].axvline((eq_steps)/part_num)
axe[0].text(eq_steps,part_num*0.6,'{0:4.2e} Compression steps'.format(co_steps))
axe[0].axvline((eq_steps+co_steps)/part_num)
axe[0].text(eq_steps+co_steps,part_num*0.4,'{0:4.2e} Production steps'.format(pr_steps))
axe[0].axvline((eq_steps+co_steps+pr_steps)/part_num)
axe[0].set_title('N={0}, phi_final = {1}, T= 0.75'.format(part_num,rho))
axe[0].set_ylabel('Number of solid particles')
axe[0].set_xlabel('step number per particle')
axe[0].set_ylim(0,part_num)

axe[1].errorbar(steps/part_num,msd,s_msd,fmt='.')
#axe[1].plot(steps/part_num,msd)
axe[1].set_ylabel('MSD (unwrapped)')
axe[1].set_xlabel('steps per particle')
axe[1].axvline(0)
axe[1].axvline((eq_steps)/part_num)
axe[1].axvline((eq_steps+co_steps)/part_num)
axe[1].axvline((eq_steps+co_steps+pr_steps)/part_num)

fig.savefig('overview_{0}_{1}_{2}.pdf'.format(fdate,part_num,int(rho*1000)))



