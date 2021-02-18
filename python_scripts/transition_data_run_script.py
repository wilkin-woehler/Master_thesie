import subprocess as sub
import numpy as np
import time
import os
import matplotlib.pyplot as plt
from datetime import date
import h5py

from mpl_toolkits.mplot3d import Axes3D

#removes files for cleanup
def remove_file(file_name):
	try:
		os.remove(file_name)
	except:
		print("no "+file_name+" yet")


class simulation:
	def __init__(self,N_part,phi_init,phi_final,step_array,seed):
		self.N_part = N_part
		self.phi_init = phi_init
		self.phi_final = phi_final
		self.step_array = step_array
		self.seed = seed
		self.tau = 0

		steps_eq,steps_pr,sweep_eq,sweep_co,sweep_pr = self.step_array[0],self.step_array[1],self.step_array[2],self.step_array[3],self.step_array[4]

		self.infoString = "_{0}_{1}_{2}_eq{3}k_pr{4}k_{5}".format(self.N_part,int(self.phi_init*1000),int(self.phi_final*1000),int(steps_eq*1e-3),int(steps_pr*1e-3),self.seed)

		self.in_file_name= "in_files/InFile"+self.infoString
		self.terminal_out_file_name = "output/terminal"+self.infoString+".txt"
		self.ovito_file_name ="ovito_files/ovito"+self.infoString+".dump"
		self.out_file_name = "outfiles/Out"+self.infoString+".h5"
		self.init_file_name = "initfiles/Initfile"+self.infoString+".h5"


	def gen_in_out_files(self):
        #try cleaning, ovito_files, out_files, output, init_files and in_files

		remove_file(self.in_file_name)
		remove_file(self.terminal_out_file_name)
		remove_file(self.ovito_file_name)
		remove_file(self.out_file_name)
		remove_file(self.init_file_name)	

	#generate infile
		in_file = open(self.in_file_name,"w+")
		in_file.write("Seed "+str(self.seed)+"\n")
		in_file.write("RNGInc 0\n")
		in_file.write("Dimension 3\n")
		in_file.write("Temperature 1.00\n")
		in_file.write("Number "+str(self.N_part)+"\n")
		in_file.write("DensityFinal "+str(self.phi_final)+"\n")
		in_file.write("DensityInital "+str(self.phi_init)+"\n")
		in_file.write("EqSteps "+str(steps_eq)+"\n")
		in_file.write("PrSteps "+str(steps_pr)+"\n")
		in_file.write("EqEvalSteps "+str(int(sweep_eq))+"\n")
		in_file.write("CoEvalSteps "+str(int(sweep_co))+"\n")
		in_file.write("PrEvalSteps "+str(int(sweep_pr))+"\n")
		in_file.write("EqSnapSteps "+str(1000000)+"\n")
		in_file.write("CoSnapSteps "+str(1000000)+"\n")
		in_file.write("PrSnapSteps "+str(1000000)+"\n")
		in_file.write("HSDiameter 1\n")
		in_file.write("ParticleMass 1\n")
		in_file.write("MaxNeigh 12\n")
		in_file.write("CosCutoff -0.986\n")
		in_file.write("Solq6q6 0.6\n")
		in_file.write("MinSolBonds 8\n")
		in_file.close()
	
	#ovito_file, out_file, init_file are created in cpp simulation
	#terminal_out_file is created by write_out_sim		

	def run_sim(self):
		self.gen_in_out_files()
		print("running a simulation for: ",self.in_file_name)
		self.simulation = sub.Popen(["./bin/arcturus.BasicEDMD.x", "-i", "init", "-f", self.in_file_name, "-o", "Out.h5"], stdout = sub.PIPE, text=True)

	def write_out_sim(self):
		out,err = self.simulation.communicate()
		out_file=open(self.terminal_out_file_name,"w+")
		out_file.write(out)
		out_file.close()

	def evaluate_sim(self):
		filename = "outfiles/Out"+self.infoString+".h5"
		print('read in file: '+filename)
		try:
			read_file=h5py.File(filename, "r")
		except():
			return


		#read in watchdata from measurements
		data=read_file['observables/watchtable'][()]
		steps=np.zeros(len(data))
		times=np.zeros(len(data))
		solids=np.zeros(len(data))

		for i in range(len(data)):
			steps[i]=data[i][0]
			times[i]=data[i][1]
			solids[i]=data[i][2]
		

		off=1
		fig,axe = plt.subplots(3)
		axe[0].scatter(times[off:-off],solids[off:-off]/N,s=1)

		axe[0].set_ylabel('Solids fration')
		axe[0].set_ylim(0,1)

		axe[1].scatter(times[off+1:-off],steps[off+1:-off]-steps[off:-off-1],s=1)
		axe[1].set_ylabel('steps step')
		axe[1].axhline(0)

		axe[2].scatter(times[off+1:-off],times[off+1:-off]-times[off:-off-1],s=1)
		axe[2].set_ylim(0,max(times[off+1:-off]-times[off:-off-1]))
		axe[2].set_ylabel('time step')

		fig.savefig("evaluations/Transition_for_Out"+self.infoString+".png")

		fig.tight_layout()
		
		plt.close('all')



if __name__ =="__main__":
	########################### General Numbers and Switches ##############################
	#For use of ths script, there are two switches to either resimulate a series, or do a single_evaluation of all simulations.
	# Furthermore the Numbers can be choosen freely when resimulating, but have to be choosen as during simulation when only evaluating. 
	simulate =True # switch for re simulating, only executes the evaluation if set to False
	single_eval=True

	#Frame numbers during simulation (snapshots). 
	N_eq = 30 #30
	N_pr = 250 #250 

	#number of fcc unit cells, and correspondig particle number sampled in the series
	N_cell = np.linspace(16,16,1, dtype =int)
	N_part = 4 * N_cell**3
	print('particle numbers to be sampled')
	print(N_part)

	#densities sampled in the series
	phi_init=0.45
	phi=np.linspace(0.54,0.54,1)
	print('densities to be sampled')
	print(phi)

	#statistics of each simulation
	stat_n=81 #sets the number of different seeds used for one set of N and rho
	seed_offset = 581 #determines the first seed used, i.e simulations with seeds 
                          #from [seed_offset, seed_offset+stat_n-1] are carried out

	#should be choosen to fit the desired number of steps, and resolution
	sweep_eq =  300 #300
	sweep_co =  300
	sweep_pr =  300 #300

	#steps to be executed
	steps_eq = sweep_eq * N_eq
	steps_pr = sweep_pr * N_pr

	#collection of step numbers, for later use
	step_array=[steps_eq,steps_pr,sweep_eq,sweep_co,sweep_pr]

	#preparing the simulation list with all simulations to be used
	simulation_list=[]
	for N in N_part:  
		for phi_final in phi:
			for seed in range(seed_offset,stat_n+seed_offset):
				simulation_list.append(simulation(N,phi_init,phi_final,step_array,seed))

	N_process=3 # number of processes run in parallel during simulation time
	

	#bad looking procedure to work through the simulation list, only executing N_process simulations at ounce
	active=np.zeros(len(simulation_list))
	sim_i=0		
	while (sim_i < len(simulation_list)):
		if (sum(active) < N_process):
			if simulate:
				simulation_list[sim_i].run_sim()
				active[sim_i] = 1
			sim_i+=1
		else:
			time.sleep(20)
			for k in range(sim_i):
				if (active[k]==1) & (simulation_list[k].simulation.poll()==0) :
					print("1give out simulation ", k)
					simulation_list[k].write_out_sim()
					if simulate :
						try:
							simulation_list[k].evaluate_sim()
						except():
							print('1some problem with simulation ',k)
					active[k]=0

	while(sum(active)!=0):
		time.sleep(20) #every few seconds the program looks if a process has finished
		for k in range(len(simulation_list)):
			if (active[k]==1) & (simulation_list[k].simulation.poll()==0) :
				print("2give out simulation ", k)
				simulation_list[k].write_out_sim()
				if simulate :
					try:
						simulation_list[k].tau = simulation_list[k].evaluate_sim()
					except():
						print('2some problem with simulation ',k)
				active[k]=0




