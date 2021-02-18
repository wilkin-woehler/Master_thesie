import matplotlib.pyplot as plt
import scipy.optimize as so
import numpy as np

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

time=time[20:-40]
lc=lc[20:-40]

arg=np.polyfit(lc,time,2)
f=np.poly1d(arg)

x=np.linspace(0,300000,100)

plt.scatter(lc,time,s=4,alpha=0.9)

plt.scatter(x,f(x),s=2,alpha=0.6)
plt.ylabel('seconds')
plt.xlabel('N largest cluster')
plt.title('q6q6 calculation costs')
#plt.show()
plt.savefig('q6q6_calculation_time.pdf')
