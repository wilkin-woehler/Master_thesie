#small tool to approximate wallltime of EDMD-Hard-Sphere-Simulations running on the NEMO-Cluster
import numpy as np

class Measurement():
	def __init__(self,N,eq_steps,pr_steps,meas_steps):
		self.N=N
		self.eq_steps=eq_steps
		self.pr_steps=pr_steps
		self.N_step=self.eq_steps+self.pr_steps
		self.meas_steps=meas_steps
		
		self.N_meas=self.N_step/self.meas_steps
		print("N = {0}, eq_steps = {1}, pr_steps = {2}, meas_steps = {3}".format(N,eq_steps,pr_steps,meas_steps))
		print("Simulation of T = {0:4.2f} deltat".format(self.N_step/67.6))

	def runtime(self):
		#working quite good (on a fast node it might be an overestimate)
		t_eval=7.64e-6*self.N
		t_step=10.51e-6*self.N
		T_eval=t_eval*self.N_meas/3600
		T_step=t_step*self.N_step/3600
		print("Measure time: {0:4.2f} h + step time: {1:4.2f} h = {2:4.2f}".format(T_eval,T_step,T_eval+T_step))

	def filesize(self):
		#well, working mostly good
		reset_size=235e-9*self.N
		setup_size=49e-9*self.N
		measure_size=0
		print("Filesize per saving: Reset {0:4.2f} GB, Snapshot {1:4.2f} GB".format(reset_size,setup_size))
		


#input parameters
N=np.linspace(64,64,1)**3*4
eq_steps=np.linspace(1000,1000,1)
pr_steps=np.linspace(20,20,1)*1e3
meas_steps=np.linspace(10,20,2)


print("================================")
for i in range(len(N)):
	for j in range(len(eq_steps)):
		for k in range(len(pr_steps)):
			for l in range(len(meas_steps)):
				temp=Measurement(N[i],eq_steps[j],pr_steps[k],meas_steps[l])	
				temp.runtime()
				temp.filesize()		
				print("----------------------------------")





#msub ./bin/login1.nemo.privat.BasicEDMD.x -i init -f in_files/InFile_450_545_4000_eq0k_pr3k_1



