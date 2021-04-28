import matplotlib.pyplot as plt
import scipy.optimize as so
import numpy as np



SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title




data=open("cluster_calculation_times.txt",'r')
time=[]
lc=[]
for line in data.readlines():
     #line=data.readline()
     #print(line)
     split_line=line.split(' ')
     if split_line[0]=='Cluster':
          #print(split_line)
          time.append(float(split_line[3]))
          lc.append(float(split_line[6][:-1]))

time=np.array(time[20:-40])
lc=np.array(lc[20:-40])

arg=np.polyfit(lc,time,2)
f=np.poly1d(arg)

x=np.linspace(0,300000,100)

plt.scatter(lc*1e-4,time,s=4,alpha=0.9, label ='Cluster finding routine')

plt.scatter(x*1e-4,f(x),s=2,alpha=0.6, label='Quadratic best fit')
plt.ylabel('CPU time in seconds')
plt.xlabel(r'Largest cluster size $N/10^4$')
#plt.title('q6q6 calculation costs')
#plt.show()
plt.legend()
plt.savefig('q6q6_calculation_time.pdf')
plt.savefig('../../../plots/q6q6_calculation_time.pdf')
