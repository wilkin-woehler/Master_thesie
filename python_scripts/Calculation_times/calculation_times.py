import numpy as np
import matplotlib.pyplot as plt

#prelimnary data taking for deciding on aim
N=np.array([32,32,32,32,108,108,108,108,256,256,256,256,500,500,500,500,864,864,864,864])*1e3

T_OE=np.array([139,130,13,25,141,19,58,228,108,80,69,89,115,145,0,127,0,0,0,0])
STEPS_OE=np.array([578,558,46,92,142,20,60,240,46,34,28,36,24,30,1,28,1,1,1,1])

T_EQ=np.array([2716,2649,2736,2691,10238,9473,9560,9621,23300,23400,24100,24260,47500,47800,48200,44800,81400,81800,82600,84200])
STEPS_EQ=np.array([10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000])

t_nuc=np.array([8973,8670,850,1511,2262,4229,1060,3735,789,609,459,639,429,519,519,399,459,398,519,519])
t_eq=np.ones(20)*287.4

tau=t_nuc-t_eq




fig,axe=plt.subplots(4,figsize=(5,10))
axe[0].set_title('eta=0.536 on NEMO')
axe[0].set_ylabel('CPU_time/Step')
axe[0].scatter(N,T_OE/STEPS_OE,c='b',label='(OE)',s=15,alpha=0.5,marker='x')
axe[0].scatter(N,T_EQ/STEPS_EQ,c='g',label='(EQ)',s=15,alpha=0.5,marker='x')
axe[0].legend()

axe[1].set_ylabel('CPU_time/Event (mus)')
axe[1].scatter(N,T_OE/STEPS_OE/N*1e6,c='b',label='(OE)',s=15,alpha=0.5,marker='x')
axe[1].scatter(N,T_EQ/STEPS_EQ/N*1e6,c='g',label='(EQ)',s=15,alpha=0.5,marker='x')
#axe[1].set_ylim(0,1.1*max(max(T_EQ/STEPS_EQ/N),max(T_OE/STEPS_OE/N)))
axe[1].legend()

axe[2].set_ylabel('-log10((tau * N))')
axe[2].scatter(N,np.log10(1/(tau*N)))
#axe[2].set_ylim(0,max(np.log(1/(tau*N))))
axe[2].set_xlabel('N')

axe[3].set_ylabel('T_OE')
axe[3].scatter(N,T_OE)
#axe[3].set_ylim(0,)
axe[3].set_xlabel('N')




fig.tight_layout()
fig.savefig('calculation_times_to_nucleation.pdf')









#NEMO data taking with mostly done code
N2=np.array([64,70,80,90,100,110,120,130])**3*4

ts64=np.mean([9.65,9.66,11.02,11.03,11.05,11.01,10.99,11.09,11.13,10.91,11.15])
ts70=np.mean([14.5,13.5,14.5,14.73,14.52,14.23])
ts80=np.mean([20.67,21.29,21.58,22.04,21.60,21.48,22.27])
ts90=np.mean([35.63,36.00,36.03,36.17,36.18])
ts100=np.mean([45.29,45.98,46.25,46.76,45.94])
ts110=np.mean([53.95,52.58,52.38])
ts120=np.mean([73.81,73.99,74.14,71.97])
ts130=np.mean([93.05,94.93,93.72,93.57])

tm64=np.mean([7.966,8.00,7.910,7.982,7.817,7.87,7.6,7.43,7.43,7.42,7.53])
tm70=np.mean([10.22,10.36,10.22,9.94,10.14])
tm80=np.mean([15.89,16.34,15.33,15.45])
tm90=np.mean([24.73,24.55,25.08,24.75])
tm100=np.mean([31.66,33,31.9,32.1,31.8])
tm110=np.mean([40.52,41.98,40.99,39.6,38.76])
tm120=np.mean([53.53,52.91,52.69,52.35,54.45])
tm130=np.mean([69.35,66.97,65.92,66.97])

T_STEP2=np.array([ts64,ts70,ts80,ts90,ts100,ts110,ts120,ts130])
T_MEAS2=np.array([tm64,tm70,tm80,tm90,tm100,tm110,tm120,tm130])

p_step=np.polyfit(N2,T_STEP2,1)
p_meas=np.polyfit(N2,T_MEAS2,1)

step_linreg=np.poly1d(p_step)
meas_linreg=np.poly1d(p_meas)
print(p_step)
print(p_meas)

plt.figure()
plt.scatter(N2,T_STEP2,label="Step times = {0:4.2f} mus/step/particle ".format(p_step[0]*1e6))
plt.scatter(N2,T_MEAS2, label="Measure times = {0:4.2f} mus/meas/particle ".format(p_meas[0]*1e6))
plt.plot(N2,step_linreg(N2))
plt.plot(N2,meas_linreg(N2))
plt.legend()
plt.xlabel("N")
plt.ylabel("Time in seconds")
plt.savefig("calculation_times_measure_step.pdf")





N3=np.array([30,50,64])**3*4
base_file=np.array([3.2,20.2,48.8]) #mb
sim20file=np.array([512.324824,2378.084472,4983.284056]) #mb
simplus20file =np.array([512.580760,2378.382192,4984.703592]) #mb

comparison=[64**3*4,3600/8]

p_base=np.polyfit(N3,base_file,1)
p_reset=np.polyfit(N3,(sim20file-base_file)/20,1)

base_linreg=np.poly1d(p_base)
reset_linreg=np.poly1d(p_reset)


plt.figure()
plt.scatter(N3,base_file,label = "setup size = {0:4.0f} byte/snap/particle ".format(p_base[0]*1e6))
plt.scatter(N3,(sim20file- base_file)/20, label = "reset_sim size = {0:4.0f} byte/reset/particle ".format(p_reset[0]*1e6))
plt.plot(N3,base_linreg(N3))
plt.plot(N3,reset_linreg(N3))
#plt.scatter(N3,(simplus20file-sim20file)/20*1000*10, label = "measurement size *10,000")
#plt.scatter(comparison[0],comparison[1])
plt.axhline(0)
plt.xlabel("N")
plt.ylabel("File size in mb")
plt.legend()
plt.savefig("File_size.pdf")




