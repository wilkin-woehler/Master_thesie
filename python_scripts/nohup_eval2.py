import numpy as np
import matplotlib.pyplot as plt
import os as os
from numpy.fft import fft, fftfreq

text_file=open('nohup_nucleation_1m_low_density.txt','r')
#text_file=open('nohup_glass620_run_4000.txt','r')
#text_file2=open('nohup_glass620_run_4000_11.txt','r')
#text_file3=open('nohup_glass620_run.txt','r')

sf=[]
lc=[]
ls_sorted=[[],[],[]]
ti=[]
n1,n2= 6150,6820
t=0
for line in text_file.readlines():
	ls=line.split()
	#print(ls)
	if len(ls) > 3:
		if ls[0] =='finding' and ls[1] =='a' and ls[2] =='solid':
#			print((eval(ls[5])))
#			sf.append(eval(ls[5]))
			sf.append(int(ls[5][:4]))
			if int(ls[5][:4]) < n1:
				t=1
			elif int(ls[5][:4]) > n1 and int(ls[5][:4]) < n2:
				t=2
			elif int(ls[5][:4]) > n2:
				t=3

		if ls[0] =='and' and ls[1] =='a' and ls[2] =='largest':
			print(ls[6])
			lc.append(int(ls[6]))
			if t==1:
				ls_sorted[0].append(int(ls[6]))
			elif t==2:
				ls_sorted[1].append(int(ls[6]))
			elif t==3:
				ls_sorted[2].append(int(ls[6]))
		if ls[0] =='measuring' and ls[1] =='at' and ls[2] =='system':
			ti.append(float(ls[4]))


fig,axe = plt.subplots(4,figsize=(8,10),sharey=True)
#axe[0].scatter(np.zeros(len(sf)),sf)
#axe[0].hist(sf,bins=100)
#axe[0].axvline(n1)
#axe[0].axvline(n2)

mm=300
lc_mean=np.convolve(lc,np.ones(mm)/mm,mode='valid')
#lc_mean2=np.convolve(lc2,np.ones(mm)/mm,mode='valid')
#lc_mean3=np.convolve(lc3,np.ones(mm)/mm,mode='valid')
sf_mean=np.convolve(sf,np.ones(mm)/mm,mode='valid')*4000
#sf_mean2=np.convolve(sf2,np.ones(mm)/mm,mode='valid')*4000
#sf_mean3=np.convolve(sf3,np.ones(mm)/mm,mode='valid')*186624
#axe[0].set_title('N=4000, eta=0.62, seed=10/11')
#axe[2].set_title('N=186624, eta=0.62')
#axe[0].scatter(np.arange(len(lc)-mm+1), lc_mean, s=2, label='largest cluster')
#axe[0].scatter(np.arange(len(sf)-mm+1), sf_mean, s=2, label='solid fraction')
#axe[0].scatter(ti[:-mm+1], lc_mean, s=2, label='largest cluster')
#axe[0].scatter(ti[:-mm+1], sf_mean, s=2, label='solid fraction')



#axe[0].set_ylabel('#particles')
#axe[1].set_ylabel('#particles')
#axe[2].set_ylabel('#particles')
#axe[2].set_xlabel('time/delta_t')
#axe[1].scatter(np.arange(len(lc2)-mm+1),lc_mean2,s=2, label='largest cluster')
#axe[1].scatter(np.arange(len(sf2)-mm+1),sf_mean2,s=2, label='solid fraction')
#axe[1].scatter(ti2[:-mm+1],lc_mean2,s=2, label='largest cluster')
#axe[1].scatter(ti2[:-mm+1],sf_mean2,s=2, label='solid fraction')

#axe[2].scatter(ti3[:-mm+1],lc_mean3,s=2, label='largest cluster')
#axe[2].scatter(ti3[:-mm+1],sf_mean3,s=2, label='solid fraction')

#axe[0].legend()
#axe[1].legend()
#axe[2].legend()
axe[0].hist(ls_sorted[0],color='g',alpha=0.3,bins=np.arange(1,150))
axe[0].hist(ls_sorted[1],color='b',alpha=0.3,bins=np.arange(1,150))
axe[0].hist(ls_sorted[2],color='r',alpha=0.3,bins=np.arange(1,150))

for i in range(3):
	temp=np.convolve(ls_sorted[i],np.ones(mm)/mm,mode='valid')
	axe[1+i].scatter(np.arange(len(ls_sorted[i])),ls_sorted[i],c='C{0}'.format(i),label = 't={0}'.format(i),s=3, alpha = 0.5)
	axe[1+i].plot(np.arange(len(ls_sorted[i])-mm+1),temp,c='C{0}'.format(i),label = 't={0}'.format(i))

	#for j in range(3):	
	#	temp=np.convolve(ls_sorted[i][j::3],np.ones(mm)/mm,mode='valid')
	#	axe[1+i].scatter(np.arange(len(ls_sorted[i]))[j::3],ls_sorted[i][j::3],c='C{0}'.format(j),label = 't={0}'.format(j),s=3, alpha = 0.5)
	#	axe[1+i].plot(np.arange(len(ls_sorted[i]))[j:(-mm+1)*3:3],temp,c='C{0}'.format(j),label = 't={0}'.format(i))

	#axe[1+i].set_ylim(0,max(temp))
	#np.convolve(ls_sorted[i],np.ones(mm)/mm,mode='valid')
	#axe[1+i].legend()
	#F=fft(ls_sorted[i]).real
	#f=fftfreq(len(F),1)
	#f=f[1:int(len(f)/2)-mm+1]
	#F=np.convolve(F[1:int(len(F)/2)],np.ones(mm)/mm, mode='valid')
	
	#plt.figure()
	#plt.scatter(f,F)
	#plt.savefig('FFT_nucleation_{0}.pdf'.format(i))

fig.tight_layout()
fig.savefig('Nucleation_nohup_overview.png')
#print(np.max(ls_sorted[i]))
plt.show()


F=fft(ls_sorted[0]).real
f=fftfreq(len(F),1/len(F))
f=f[1:int(len(f)/2)]
F=F[1:int(len(F)/2)]

plt.figure()
plt.scatter(f,F)
plt.savefig('FFT_nucleation.pdf')


