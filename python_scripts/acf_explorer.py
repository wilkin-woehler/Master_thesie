#author : wilkin woehler
#date : 11.05.2020

import h5py
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from scipy.stats import norm


#read in h5(md) file from which data is supposed to be read

filename = 'Out_1048576_450_530_eq1k_pr1k_10.h5'
read_file=h5py.File(filename, "r")


data = read_file['observables/watchtable']

largest_cluster=[]
time=[]
for i in range(100,len(data)):
	largest_cluster.append(data[i][3])
	time.append(data[i][1])
	#print(data[i][1])

#N_end=100001
#x=np.linspace(0,N_end,N_end+1)
#y=line(x,1,0)

#plt.figure()
#plt.scatter(x,y)
#plt.savefig('acf_source.png')

#N_acf=N_end-1

'''
s1=time.time()

acf=np.zeros(N_acf)
for i in range(N_acf):
	for j in range(len(y)-i):
		acf[i]+=y[j]*y[j+i]
	acf[i]=acf[i]/(len(y)-i)
e1=time.time()

acf_normed=acf/np.std(y)**2
'''

largest_cluster=np.array(largest_cluster)

base_mm2=25

mm2_a=base_mm2
mm_a=mm2_a*2
kernel_a=np.ones(mm_a)/mm_a
mlc_a=np.convolve(largest_cluster,kernel_a,mode='valid')


#daten=gib_mir
#halbe_fenster_laenge=25 # Fensterlänge 
#fenster_laenge=halbe_fenster_laenge*2
#kern=np.ones(fenster_laenge)/fenster_laenge
#gemittelte_daten=np.convolve(daten, kern, mode='valid')

#gemittelte_daten ~ daten[halbe_fenster_laenge:-halbe_fenster_laenge+1]
#fluktuationen = daten[halbe_fenster_laenge:-halbe_fenster_laenge+1] - gemittelte_daten


mm2_b=base_mm2*5
mm_b=base_mm2*10
kernel_b=norm.pdf(np.linspace(-mm2_b,mm2_b,mm_b),loc=0, scale=mm_b/20)
mlc_b=np.convolve(largest_cluster,kernel_b,mode='valid')


mlc_c=np.mean(largest_cluster)
#plt.scatter(np.linspace(-mm2_a,mm2_a,mm_a),kernel_a)
#plt.scatter(np.linspace(-mm2_b,mm2_b,mm_b),kernel_b)


largest_cluster_fluc_a=largest_cluster[mm2_a:-mm2_a+1]-mlc_a
largest_cluster_fluc_b=largest_cluster[mm2_b:-mm2_b+1]-mlc_b
largest_cluster_fluc_c=largest_cluster-mlc_c


N_acf_a=len(largest_cluster_fluc_a)-1
N_acf_b=len(largest_cluster_fluc_b)-1
N_acf_c=len(largest_cluster_fluc_c)-1

time_a= np.array(time)[mm2_a:-mm2_a+1]
time_b= np.array(time)[mm2_b:-mm2_b+1]
time_c= np.array(time)

tau_a= np.zeros(N_acf_a)
#stau= np.zeros(N_acf)
for i in range(1,N_acf_a):
	delta_time=time_a[i:]-time_a[:-i]
	tau_a[i]=np.mean(delta_time)
	#stau[i]=np.std(delta_time)


tau_b= np.zeros(N_acf_b)
#stau= np.zeros(N_acf)
for i in range(1,N_acf_b):
	delta_time=time_b[i:]-time_b[:-i]
	tau_b[i]=np.mean(delta_time)

tau_c= np.zeros(N_acf_c)
#stau= np.zeros(N_acf)
for i in range(1,N_acf_c):
	delta_time=time_c[i:]-time_c[:-i]
	tau_c[i]=np.mean(delta_time)

acf_a=np.zeros(N_acf_a)
acf_a[0]=np.mean(largest_cluster_fluc_a**2)
for i in range(1,N_acf_a):
	acf_a[i]=np.mean(largest_cluster_fluc_a[i:]*largest_cluster_fluc_a[:-i])


acf_b=np.zeros(N_acf_b)
acf_b[0]=np.mean(largest_cluster_fluc_b**2)
for i in range(1,N_acf_b):
	acf_b[i]=np.mean(largest_cluster_fluc_b[i:]*largest_cluster_fluc_b[:-i])

acf_c=np.zeros(N_acf_c)
acf_c[0]=np.mean(largest_cluster_fluc_c**2)
for i in range(1,N_acf_c):
	acf_c[i]=np.mean(largest_cluster_fluc_c[i:]*largest_cluster_fluc_c[:-i])


fig,axe = plt.subplots(3)
axe[0].scatter(time, largest_cluster,c='k',alpha=0.3)
axe[0].plot(time_c, mlc_c*np.ones(len(time_c)) ,c='k')
axe[0].plot(time_a, mlc_a,c='blue')
axe[0].plot(time_b, mlc_b,c='green')

axe[1].scatter(time_c,largest_cluster_fluc_c,c='k')
axe[1].scatter(time_a,largest_cluster_fluc_a,c='blue')
axe[1].scatter(time_b,largest_cluster_fluc_b,c='green')


axe[2].scatter(tau_a,acf_a,c='blue')
axe[2].scatter(tau_b,acf_b,c='green')
axe[2].scatter(tau_c,acf_c,c='k')

axe[2].axhline(0)
#axe[2].errorbar(np.arange(N_acf),tau,stau*1e5)
plt.show()
