import subprocess as sub
import numpy as np
import time
import os
import matplotlib.pyplot as plt
from datetime import date
import h5py

from mpl_toolkits.mplot3d import Axes3D

#remove files for cleanup
def remove_file(file_name):
	try:
		os.remove(file_name)
	except:
		print("no "+file_name+" yet")

def find_t(val_y,data_x,data_y):
	index = np.arange(len(data_y))[np.abs(data_y-val_y)==np.min(np.abs(data_y-val_y))]

	if len(index)==1 :
		return(index,data_x[index[0]])
	elif len(index)==0 :
		print('no corresponding value found. (Should not happen)')
		return(0,0)
	else:
		print('more than one value is closest to the given value, picking one from:')
		print('indexes: ',index)
		for i in index:
			print(data_x[i],data_y[i],val_y)
		return(index[-1],data_x[index[-1]])


#simple linear regression routine
def lineare_regression(x,y,ey):
    '''

    Lineare Regression.

    Parameters
    ----------
    x : array_like
        x-Werte der Datenpunkte
    y : array_like
        y-Werte der Datenpunkte
    ey : array_like
        Fehler auf die y-Werte der Datenpunkte

    Diese Funktion benoetigt als Argumente drei Listen:
    x-Werte, y-Werte sowie eine mit den Fehlern der y-Werte.
    Sie fittet eine Gerade an die Werte und gibt die
    Steigung a und y-Achsenverschiebung b mit Fehlern
    sowie das chi^2 und die Korrelation von a und b
    als Liste aus in der Reihenfolge
    [a, ea, b, eb, chiq, cov].
    '''

    s   = sum(1./ey**2)
    sx  = sum(x/ey**2)
    sy  = sum(y/ey**2)
    sxx = sum(x**2/ey**2)
    sxy = sum(x*y/ey**2)
    delta = s*sxx-sx*sx
    b   = (sxx*sy-sx*sxy)/delta
    a   = (s*sxy-sx*sy)/delta
    eb  = np.sqrt(sxx/delta)
    ea  = np.sqrt(s/delta)
    cov = -sx/delta
    corr = cov/(ea*eb)
    chiq  = sum(((y-(a*x+b))/ey)**2)

    return(a,ea,b,eb,chiq,corr)

#routine to read out a histogram(name) from h5md file(read_file)
def read_out_hist(name, read_file):
	data=read_file['observables/'+name][()]
	count=[]
	s_count=[]
	bins=[]
	step_check=-1
	for i in range(len(data)):
		if (step_check!=data[i][0]):
			count.append([])
			s_count.append([])
			bins.append([])
			step_check=data[i][0]
		count[-1].append(data[i][3])
		s_count[-1].append(data[i][4])
		bins[-1].append(data[i][2])
	return(count, bins, s_count)

#Enskog Diffusion constant
def D_E(rho,g_d):
	m=1
	k_bT=1
	d=1
	return(3/8/rho/d**2/g_d*(k_bT/np.pi/m)**(1/2))


class simulation:
	def __init__(self,N_part,phi_init,phi_final,step_array,seed):
		self.N_part = N_part
		self.phi_init = phi_init
		self.phi_final = phi_final
		self.step_array = step_array
		self.seed = seed
		self.tau = 0
		self.transition_width = 0

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
		in_file.write("EqSnapSteps "+str(sweep_eq)+"\n")
		in_file.write("CoSnapSteps "+str(sweep_co)+"\n")
		in_file.write("PrSnapSteps "+str(sweep_pr)+"\n")
		in_file.write("HSDiameter 1\n")
		in_file.write("ParticleMass 1\n")
		in_file.write("MaxNeigh 12\n")
		in_file.write("CosCutoff -0.986\n")
		in_file.write("Solq6q6 0.6\n")
		in_file.write("MinSolBonds 8\n")
		in_file.close()
		

	def run_sim(self):
		self.gen_in_out_files()
		print("running a simulation for: ",self.in_file_name)
		#self.simulation = sub.Popen(["./bin2/vega.BasicEDMD.x", "-i", "init", "-f", self.in_file_name, "-o", "Out.h5"], stdout = sub.PIPE, text=True)
		self.simulation = sub.Popen(["./bin2/vega.BasicEDMD.x", "-i", "init", "-f", self.in_file_name, "-o", "Out.h5"])
		#self.simulation = sub.Popen(["./bin2/arcturus.BasicEDMD.x", "-i", "init", "-f", self.in_file_name, "-o", "Out.h5"])

	def write_out_sim(self):
		out,err = self.simulation.communicate()
		out_file=open(self.terminal_out_file_name,"w+")
		out_file.write(out)
		out_file.close()

	def evaluate_sim(self):
		#steps_eq,steps_pr,sweep_eq,sweep_co,sweep_pr = self.step_array[0],self.step_array[1],self.step_array[2],self.step_array[3],self.step_array[4]
		full_eval=False
		filename = "outfiles/Out"+self.infoString+".h5"
		print('read in file: '+filename)
		try:
			read_file=h5py.File(filename, "r")
		except():
			return

		#read in histograms for overview plot
		if full_eval:
			g_r,r,s_g_r = read_out_hist('RDF_hist', read_file)
			msd_count,msd, s_msd_count = read_out_hist('msd_hist', read_file)
			q6q6_count, q6q6, s_q6q6_count= read_out_hist('q6q6_hist', read_file)
			sol_neigh_count, sol_neigh, s_sol_neigh_count = read_out_hist('sol_neigh_hist', read_file)
			vel_count, vel, s_vel_count = read_out_hist('Velocity_hist', read_file)
			vel=np.array(vel)
			vel_count=np.array(vel_count)

		#obtain step numbers and particle numbers
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


		#read in watchdata from measurements
		data=read_file['observables/watchtable'][()]
		steps=np.zeros(len(data))
		times=np.zeros(len(data))
		solids=np.zeros(len(data))
		clus_size=np.zeros(len(data))
		msd_mean=np.zeros(len(data))
		s_msd_mean=np.zeros(len(data))
		vel_2=np.zeros(len(data))

		for i in range(len(data)):
			steps[i]=data[i][0]
			times[i]=data[i][1]
			solids[i]=data[i][2]
			clus_size[i]=data[i][3]
			msd_mean[i]=data[i][4]
			s_msd_mean[i]=data[i][5]
			vel_2[i]=data[i][6]

		
		if full_eval:
			print("number of pr snapshots    : {0}".format(len(q6q6)))
		print("number of all evaluations : {0}".format(len(msd_mean)))
		print("")

		nuc_start_time=times[np.abs(steps-(eq_steps+co_steps)) == min(np.abs(steps-(eq_steps+co_steps)))]
		nuc_end_time=times[np.abs(steps-(co_steps+pr_steps)) == min(np.abs(steps-(co_steps+pr_steps)))]
		#new calculation of nucleation time
		nuc_end_time = times[np.abs(solids-0.5*part_num) == min(np.abs(solids-0.5*part_num)) ][0]

		tau_1 = nuc_end_time-nuc_start_time


		end_val=np.mean(solids[-10:])/part_num
		index_tau, tau_2 = find_t(0.5*end_val, times, solids/part_num)
		self.tau = tau_2

		e=np.exp(1)
		start_index, transition_start = find_t(1/(1+e)*end_val, times, solids/part_num)
		print('start', transition_start)
		stop_index, transition_stop  = find_t(e/(1+e)*end_val, times, solids/part_num)
		print('stop', transition_stop)
		self.transition_width = transition_stop - transition_start

		#print(nuc_end_time)
		#print(solids)

		if full_eval:
			n=5
			ind = np.linspace(0,len(g_r)-1,n,dtype=int)

			fig, axe = plt.subplots(2,3)
			axe[0,0].set_title("g(r)")
		
			#for i in range(len(g_r)):
			for i in range(n): 
			#print(r)
				axe[0,0].errorbar(r[ind[i]],g_r[ind[i]],s_g_r[ind[i]],fmt='.', ms=1,linewidth=0.5, alpha= 0.5)
			#axe[0,0].scatter(r[i],g_r[i])
			axe[0,0].axvline(1)
			axe[0,0].axhline(1)


			axe[0,1].set_title("msd")
		#for i in range(len(msd)): #plotting the msd distributions, very flat, extending quite far
			#axe[0,1].plot(np.array(msd_count[i])*100000000,np.array(msd[i])**2)
			axe[0,1].errorbar(times, msd_mean, s_msd_mean, fmt='.')
			axe[0,1].set_ylim(0,max(msd_mean))


			axe[1,0].set_title("q6q6")
			#for i in range(len(sol_neigh)):
			for i in range(n):
				#axe[1,0].plot(q6q6[i],q6q6_count[i])
				axe[1,0].errorbar(q6q6[ind[i]],q6q6_count[ind[i]],s_q6q6_count[ind[i]],fmt='.')


			axe[1,1].set_title("sol_neigh")
			#for i in range(len(sol_neigh)):
			for i in range(n):
				#axe[1,1].plot(sol_neigh[i],sol_neigh_count[i])
				axe[1,1].errorbar(sol_neigh[ind[i]],sol_neigh_count[ind[i]],s_sol_neigh_count[ind[i]], fmt = '.')


			axe[0,2].set_title('velocity')
			#for i in range(len(vel)):
			for i in range(n):
				axe[0,2].scatter(vel_count[ind[i]]*np.max(times)*10, vel[ind[i]]**2)	
			axe[0,2].scatter(times, vel_2)

	
			axe[1,2].set_title("solids")
			axe[1,2].scatter(times,solids, s=1, label='overall solids')
			axe[1,2].scatter(times,clus_size, s=1, label='largest cluster')
			axe[1,2].axvline(nuc_start_time)
			axe[1,2].axvline(nuc_end_time) #only works if eq_steps and pr_steps are the same
			#axe[1,2].axvline(tau_2+self.transition_width/2,color='r')
			axe[1,2].axvline(transition_start,color='r')
			axe[1,2].axvline(tau_2,color='r',linestyle='dashed')
			#axe[1,2].axvline(tau_2-self.transition_width/2,color='r')
			axe[1,2].axvline(transition_stop,color='r')
			axe[1,2].legend()
			#axe[1,2].set_ylim(0,N)
			axe[1,2].set_xlim(tau_2-50,tau_2+50)
			fig.tight_layout()
			fig.savefig("evaluations/Overview_for_Out"+self.infoString+".png", dpi=500)
		else:

			plt.scatter(times,solids, s=1, label='overall solids')
			plt.scatter(times,clus_size, s=1, label='largest cluster')
			plt.axvline(nuc_start_time)
			plt.axvline(nuc_end_time) #only works if eq_steps and pr_steps are the same
			plt.axvline(tau_2+self.transition_width/2,color='r')
			plt.axvline(transition_start,color='r')
			plt.axvline(tau_2,color='r',linestyle='dashed')
			#plt.axvline(tau_2-self.transition_width/2,color='r')
			plt.axvline(transition_stop,color='r')
			plt.legend()
			plt.title("solids")
			#plt.set_ylim(0,N)
			plt.xlim(tau_2-50,tau_2+50)
			plt.tight_layout()
			plt.savefig("evaluations/Solids_Overview_for_Out"+self.infoString+".png", dpi=500)
			
		plt.close('all')
		


		#return(nuc_end_time-nuc_start_time)

def OverViewPlot():
	tau_raw=np.zeros(shape=(len(N_part),len(phi),stat_n))
	for i in range(len(N_part)):
		for j in range(len(phi)):
			for k in range(stat_n):
				tau_raw[i,j,k]=simulation_list[k+j*stat_n+i*(stat_n*len(phi))].tau
	tau_mean=np.mean(tau_raw,axis=2)
	s_tau_mean=np.std(tau_raw,axis=2)

	trans_raw=np.zeros(shape=(len(N_part),len(phi),stat_n))
	for i in range(len(N_part)):
		for j in range(len(phi)):
			for k in range(stat_n):
				trans_raw[i,j,k]=simulation_list[k+j*stat_n+i*(stat_n*len(phi))].transition_width
	trans_mean=np.mean(trans_raw,axis=2)
	s_trans_mean=np.std(trans_raw,axis=2)


	N_edges=[]
	N_edges.append(1.5*N_part[0]-0.5*N_part[1])
	for i in range(len(N_part)-1):
		N_edges.append((N_part[i+1]+N_part[i])/2)
	N_edges.append(1.5*N_part[-1]-0.5*N_part[-2])

	#print(N_edges)
	phi_edges=[]
	phi_edges.append(1.5*phi[0]-0.5*phi[1])
	for i in range(len(phi)-1):
		phi_edges.append((phi[i+1]+phi[i])/2)
	phi_edges.append(1.5*phi[-1]-0.5*phi[-2])

	#print(phi_edges)
	fig = plt.figure()
	im = plt.pcolormesh(phi_edges,N_edges, np.log((N_part*tau_mean.T).T))
	fig.colorbar(im)
	plt.ylabel('N')
	plt.xlabel('Volume fraction')
	plt.title('log (<tau> * N)')
	plt.xticks(phi,np.array(phi,dtype='str'))
	plt.yticks(N_part)
	plt.savefig("nuc_time_colormesh_N.pdf")

	fig = plt.figure()
	im = plt.pcolormesh(phi_edges,N_edges, np.log(tau_mean))
	fig.colorbar(im)
	plt.ylabel('N')
	plt.xlabel('Volume fraction')
	plt.title('log( <tau> )')
	plt.xticks(phi,np.array(phi,dtype='str'))
	plt.yticks(N_part)
	plt.savefig("nuc_time_colormesh.pdf")


	fig = plt.figure()
	im = plt.pcolormesh(phi_edges, N_edges, trans_mean)
	fig.colorbar(im)
	plt.ylabel('N')
	plt.xlabel('Volume fraction')
	plt.title('transition_width')
	plt.xticks(phi,np.array(phi,dtype='str'))
	plt.yticks(N_part)
	plt.savefig("transition_time_colormesh.pdf")


	plt.close('all')
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')


	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111, projection='3d')

	n = 100

	for i in range(len(N_part)):
		for j in range(len(phi)):
			xs = N_part[i]
			ys = phi[j]

			zs = np.log(tau_mean[i,j]*N_part[i])
			ax.scatter(xs, ys, zs, marker='x')
			zs = np.log(tau_mean[i,j]*N_part[i]+s_tau_mean[i,j]*N_part[i]/np.sqrt(stat_n))
			ax.scatter(xs, ys, zs, marker='o')
			zs = np.log(tau_mean[i,j]*N_part[i]-s_tau_mean[i,j]*N_part[i]/np.sqrt(stat_n))
			ax.scatter(xs, ys, zs, marker='o')
			
			zs2 = np.log(tau_mean[i,j])
			ax2.scatter(xs, ys, zs2, marker='x')
			zs2 = np.log(tau_mean[i,j]+s_tau_mean[i,j]/np.sqrt(stat_n))
			ax2.scatter(xs, ys, zs2, marker='o')
			zs2 = np.log(tau_mean[i,j]-s_tau_mean[i,j]/np.sqrt(stat_n))
			ax2.scatter(xs, ys, zs2, marker='o')
			#print(tau_mean[i,j],s_tau_mean[i,j])
	
	ax.set_xlabel('N')
	ax.set_ylabel('Volume fraction')
	ax.set_zlabel('log( <tau> * N)')
	ax.set_xticks(N_part)
	ax.set_xticklabels(np.array(N_part,dtype='str'))
	ax.set_yticks(phi)

	ax2.set_xlabel('N')
	ax2.set_ylabel('Volume fraction')
	ax2.set_zlabel('log( <tau> )')
	ax2.set_xticks(N_part)
	ax2.set_xticklabels(np.array(N_part,dtype='str'))
	ax2.set_yticks(phi)

	fig.savefig('nuc_time_3D_log_N')
	fig2.savefig('nuc_time_3D_log')
	plt.show()

if __name__ =="__main__":
	########################### General Numbers and Switches ##############################
	#For use of ths script, there are two switches to either resimulate a series, or do a single_evaluation of all simulations.
	# Furthermore the Numbers can be choosen freely when resimulating, but have to be choosen as during simulation when only evaluating. 
	simulate = True# switch for re simulating, only executes the evaluation if set to False
	single_eval=True

	#Frame numbers during simulation (snapshots). Have to be equal for the simple evaluation
	N_eq = 1 #10
	N_pr = 1 #200 

	#number of fcc unit cells, and correspondig particle number sampled in the series
	N_cell = np.linspace(36,36,1, dtype =int)
	N_part = 4 * N_cell**3
	print('particle numbers to be sampled')
	print(N_part)

	#densities sampled in the series
	phi=np.linspace(0.60,0.62,3)

	#statistics of each simulation, and start seed
	stat_n=1
	start_seed=10

	phi_init=0.45
	print('densities to be sampled')
	print(phi)
	#sweep size, at the moment defines measurement and snapshot. Later on only should define snapshot
	#should be choose such to fit the right number of steps
	sweep_eq =  2000 #1000
	sweep_co =  100
	sweep_pr =  1000 #1000

	#steps numbers to be executed, see above for choosing them
	steps_eq = sweep_eq * N_eq
	steps_pr = sweep_pr * N_pr

	#collection of step numbers, for later use
	step_array=[steps_eq,steps_pr,sweep_eq,sweep_co,sweep_pr]

	#preparing the simulation list with all simulation to be used
	simulation_list=[]
	for N in N_part:  
		for phi_final in phi:
			for seed in range(start_seed,stat_n+start_seed):
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
			time.sleep(5)
			for k in range(sim_i):
				if (active[k]==1) & (simulation_list[k].simulation.poll()==0) :
					print("1give out simulation ", k)
					#simulation_list[k].write_out_sim()
					if simulate :
						try:
							simulation_list[k].evaluate_sim()
						except():
							print('1some problem with simulation ',k)
					active[k]=0

	while(sum(active)!=0):
		time.sleep(5) #every few seconds the program looks if a process has finished
		for k in range(len(simulation_list)):
			if (active[k]==1) & (simulation_list[k].simulation.poll()==0) :
				print("2give out simulation ", k)
				#simulation_list[k].write_out_sim()
				if simulate :
					try:
						simulation_list[k].evaluate_sim()
					except():
						print('2some problem with simulation ',k)
				active[k]=0

	if ((single_eval) and (~simulate)) :
		for i in range(len(simulation_list)):
			try:
				simulation_list[i].evaluate_sim()
			except():
				print('3some problem with simulation',i)
		simulation_list=np.array(simulation_list)
		#np.save('auxilary_files/simulation_list',simulation_list)

	#simulation_list= np.load('auxilary_files/simulation_list.npy')

	OverViewPlot()
	

	exit()


	########################## Single Evaluations ##############################
	
	if single_eval:

		for i in range(len(phi)):
			try:
				a,ea,d=evaluate_sim(phi[i], step_array)
				diffusion_slope.append(a)
				s_diffusion_slope.append(ea)
				D_Enskog.append(d)
			except:
				print("unusuable simulation data for phi={0}".format(phi[i]))
				diffusion_slope.append(1)
				s_diffusion_slope.append(1e9)
				D_Enskog.append(1)

		diffusion_slope=np.array(diffusion_slope)/6 #/6 from 6Dt=<x^2>
		s_diffusion_slope=np.array(s_diffusion_slope)
		D_Enskog = np.array(D_Enskog)
		np.savez('auxilary_files/diffusion_data_{0}.npz'.format(N_part),diffusion_slope,s_diffusion_slope,phi,D_Enskog)
	else:
		temp_file=np.load('auxilary_files/diffusion_data_{0}.npz'.format(N_part))
		diffusion_slope=temp_file['arr_0']
		s_diffusion_slope=temp_file['arr_1']
		phi=temp_file['arr_2']
		D_Enskog=temp_file['arr_3']




