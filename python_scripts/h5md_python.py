#author : wilkin woehler
#date : 30.11.2020
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import h5py
import time as time
import os as os

class h5md_dataset:
	def __init__(self,path,N,rho_i,rho_f,eq,pr,seed,plot_parameters):
		#readin of h5mdfile
		self.infostring="_{0}_{1}_{2}_eq{3}k_pr{4}k_{5}".format(N,int(rho_i*1000),int(rho_f*1000),int(eq*1e-3),int(pr*1e-3),seed)
		self.filename = path+"Out"+self.infostring+".h5"

		print("Read in file: "+self.filename)
		self.read_file = h5py.File(self.filename, "r")
		
		#given quantities
		self.N = N
		self.rho_i=rho_i
		self.rho_f=rho_f
		self.eq_steps=eq
		if (self.eq_steps==0): #todo: required for quench step/time could be more eqaisly archived
			self.eq_steps = 100
			#print("Eq_steps could be:", self.read_file['particles/co_all_parts/position/step'][()])
			#print("Or               :", self.read_file['particles/co_all_parts/position/step'][()][0])
			#print("Eq_times could be:", self.read_file['particles/co_all_parts/position/time'][()])
			#print("Or               :", self.read_file['particles/co_all_parts/position/time'][()][0])

		#try:
			
		self.pr_steps=pr
		self.seed=seed
			
		#derived quantities
		self.box_length=(6*self.rho_f/np.pi/self.N)**(-1/3)
		
		#read quantities
		self.reset_times = self.read_file['reset_sim/particles/position/time'][()]
		#print(self.reset_times)
		#self.show_file()

		#hard coded quantites (so far):
		#stop_time must not exceed any simulation time
		#step_number should be choosen of the same order as the actual measurement number in the time frame
		self.step_number = plot_parameters[0]
		self.stop_time = plot_parameters[1]
		self.maxc = plot_parameters[2]
		self.equi_time=np.linspace(0,self.stop_time,self.step_number)

		self.time_edges = np.zeros(self.step_number+1)		
		self.time_edges[0]=(1.5*self.equi_time[0]-0.5*self.equi_time[1])
		for i in range(self.step_number-1):
			self.time_edges[i+1]=((self.equi_time[i+1]+self.equi_time[i])/2)
		self.time_edges[self.step_number]=(1.5*self.equi_time[-1]-0.5*self.equi_time[-2])

	def print_names(self,name,object_name):
		print(name)
		print(object_name)

	def show_file(self):
		self.read_file.visititems(self.print_names)

	def load_watchtable(self):
		self.wdata = self.read_file['observables/watchtable']
		self.len_wdata=len(self.wdata)
		print('watchtable of length: ', self.len_wdata)

		self.wdata_steps=np.zeros(self.len_wdata)
		self.wdata_times=np.zeros(self.len_wdata)

		#slow, because it has to iterate through the table entry by entry
		for i in range(self.len_wdata):
			self.wdata_steps[i] = self.wdata[i][0]
			self.wdata_times[i] = self.wdata[i][1]
		

		try:
			self.quench_step=self.read_file['particles/co_all_parts/position/step'][()][0]
			self.quench_time=self.read_file['particles/co_all_parts/position/time'][()][0]
			print(self.read_file['particles/co_all_parts/position/step'][()],self.read_file['particles/co_all_parts/position/time'][()])
		except:
			print("Fall back on old evaluating quench time/step")
			self.quench_step=np.min(self.wdata_steps[np.where((self.wdata_steps-self.eq_steps)>=0)])
			self.quench_time=np.min(self.wdata_times[np.where(self.wdata_steps==self.quench_step)])

		if (self.quench_step==0) or (self.quench_time==0):
			print("Fall back on old evaluating quench time/step after try block")
			self.quench_step=np.min(self.wdata_steps[np.where((self.wdata_steps-self.eq_steps)>=0)])
			self.quench_time=np.min(self.wdata_times[np.where(self.wdata_steps==self.quench_step)])	
		
		self.mod_time = self.wdata_times[self.wdata_times>=self.quench_time] - self.quench_time

	def load_watchtable2(self):
		pass

	def load_histogram(self,name):
		data=self.read_file['observables/'+name][()]
		count=[]
		pos=[]
		step_check=-1
		for i in range(len(data)):
			if (step_check!=data[i][0]):
				count.append([])
				pos.append([])
				step_check=data[i][0]
			pos[-1].append(data[i][2])
			count[-1].append(data[i][3])
		pos=np.array(pos)
		count=np.array(count)
		return(pos, count)
	
	def load_positions(self):
		pass

	def load_velocities(self):
		pass
	
	def load_ClusterDistribution(self):
		self.ClusterDistribution=self.read_file['observables/Cluster_Distribution/value']
		print("Cluster_distribution of length: ",len(self.ClusterDistribution))
		#chopping into equidistant time steps	

	def load_TensorOfGyration(self):
		self.ToG = self.read_file['observables/TensorOfGyration'][()]
		#chopping into equidistant time steps
		
	def eval_MSD(self):
		self.msd_mean=np.zeros(self.len_wdata)
		for i in range(self.len_wdata):
			self.msd_mean[i] = self.wdata[i][4]

	def equi_distance_data(self, data, write_out = False):
		equi_data=np.zeros(self.step_number)

		mod_data = data[self.wdata_times>=self.quench_time]

		delta_times=self.mod_time[1:]-self.mod_time[:-1]
		delta_data=mod_data[1:]-mod_data[:-1]

		for j in range(self.step_number):
			#search for index in times for data point just before global_time
			mask_sign = (self.mod_time-self.equi_time[j] <= 0)
			index = np.arange(len(self.mod_time))[mask_sign]
			if len(index)<=1:
				#for 0 index, this is required because only a one number array is returned
				#print(index," is only index for ", j)
				index = 0
			else:
				index = index[-1]
			#print(j, index)
			#linear interpolation of data points
			equi_data[j]=mod_data[index]+(self.equi_time[j]-self.mod_time[index])*delta_data[index]/delta_times[index]

		if write_out:
			#requires the directory to exist where it is suppossed to be written to
			filename = "largest_cluster_textfiles/Out"+self.infostring+".txt"

			file=open(filename,'w+')
			for j in range(self.step_number):
				file.write('{0}\t{1}\n'.format(self.equi_time[j],equi_data[j]))
			file.close()
		return(equi_data)


	def find_induction_time(self, dataname="LargestCluster"):
		if (dataname == "Solids"):
			data = self.solids
		elif (dataname == "LargestCluster"):
			data = self.LargestCluster

		
	
		#print()
		cut_times=self.wdata_times[np.where(self.wdata_times>=self.quench_time)]
		cut_data =data[np.where(self.wdata_times>=self.quench_time)]
		
		p1_guess=np.mean(data[-10:])	
		comp_data=np.abs(cut_data-0.5*p1_guess)
		comp_value=np.min(np.abs(cut_data-0.5*p1_guess))
		p0_guess=np.min(cut_times[np.where(comp_data==comp_value)])-self.quench_time
#		if len(p0_guess)
#		p0_guess=p0_guess[np.where(p0_guess>=self.quench_time)][0]

		#p2_guess=50
		#p3_guess=np.mean(data[:10])

		self.tau=p0_guess
		if p1_guess < 100:
			self.tau=-1

	def eval_SolidFraction(self):
		self.solids=np.zeros(self.len_wdata)
		for i in range(self.len_wdata):
			self.solids[i] = self.wdata[i][2]
		self.find_induction_time("Solids")


	def eval_LargestCluster(self):
		self.LargestCluster=np.zeros(self.len_wdata)
		for i in range(self.len_wdata):
			self.LargestCluster[i] = self.wdata[i][3]
		self.find_induction_time("LargestCluster")

	def plot_LargestCluster(self):
		fig = plt.figure()
		#plt.plot(self.wdata_times,self.solids, figure=fig)
		plt.plot(self.wdata_times,self.LargestCluster, figure=fig)
		plt.axvline(self.tau+self.quench_time, figure=fig)
		plt.axvline(self.quench_time, figure=fig, c= 'k')
		plt.title("Quench time is T_q = {0:4.2f}dt and Tau= {1:4.2f}dt".format(self.quench_time,self.tau), figure=fig)
		plt.ylim(min(self.LargestCluster[100:])-100,max(self.LargestCluster[100:])+100)
		fig.savefig("singles/lc"+self.infostring+".png")

		plt.close(fig)

	def eval_ToG(self):
		self.load_TensorOfGyration()

		self.RoG=np.zeros(len(self.ToG)) #Radius of Gyration
		self.b=np.zeros(len(self.ToG)) #asphericity
		self.c=np.zeros(len(self.ToG)) #acylindricity
		self.k2=np.zeros(len(self.ToG)) #relative shape anisotropy
		l=np.zeros(shape=(len(self.ToG),3)) #eigenvalues of ToG
		for i in range(len(self.ToG)):
			xx=self.ToG[i][0]
			yy=self.ToG[i][1]
			zz=self.ToG[i][2]
			xy=self.ToG[i][3]
			xz=self.ToG[i][4]
			yz=self.ToG[i][5]
			M=np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])
			
			eig_val, eig_vec = np.linalg.eig(M)
			try:
				M_diag= np.matmul(np.matmul(np.linalg.inv(eig_vec),M),eig_vec)
				
			except:
				M_diag = M 
				print('Matrix was not diagonalizabale : ', M)
			l=np.sort(M_diag.diagonal()) 
			#print(l)
			self.RoG[i] = np.sqrt(l[0]+l[1]+l[2])
			#print(RoG[i])
			self.b[i]   = l[2]-1/2*(l[0]+l[1])
			self.c[i]   = l[1]-l[0]
			if self.RoG[i] > 1e-10:
				self.k2[i]  = 3/2*(l[0]**2+l[1]**2+l[2]**2)/(l[0]+l[1]+l[2])**2-1/2
			else:
				self.k2[i]  = 0


	def eval_ClusterDistribution(self):
		self.load_watchtable()
		self.load_ClusterDistribution()

		l1=len(self.ClusterDistribution)
#		self.maxc=200 #maximum cluster size for display

		pnt_raw = np.zeros(shape=(self.maxc,l1))
		stable_cluster=[]
		stable_cluster_times=[]
		for i in range(l1):
			#if (l1>10) and ((i+1)%int(l1/10)==0):	
			#	print("{0:3.0f} %".format(i/l1*100))
			temp=self.ClusterDistribution[i][0]
			stable_cluster.append([])
			stable_cluster_times.append([])
			print(temp)
			for j in range(len(temp)):
				#try:
				k=temp[j][0]
				if k < self.maxc:
					pnt_raw[k,i]+=temp[j][1]
				elif (i>10):
				#else:
					for l in range(temp[j][1]):
						stable_cluster[-1].append(k)
						stable_cluster_times[-1].append(self.wdata_times[i])
				#except:
				#	print("Probablny some read error, should not occur anymore")

		pnt_raw[0,:]=0 #must be used to normalize in case multiple meausrements are included

		self.pnt = np.zeros(shape=(self.maxc,self.step_number))

		for i in range(self.maxc):
			if (i%int(self.maxc/10)==0):
				print("By now {0:2.0f} % of p(N,t) is equidistant".format(i/self.maxc*100))
			self.pnt[i,:]=self.equi_distance_data(pnt_raw[i,:])
			self.pnt[self.pnt<0.15] = 0.0


		
		
		self.stable_cluster=stable_cluster
		self.stable_cluster_times=stable_cluster_times

		stable_trajectories=[]
		stable_trajectory_times=[]
		
		sw=True
		len_stable_clusters=len(self.stable_cluster)
		print("number of stable clusters: ", len_stable_clusters)
		while sw:
			stable_trajectories.append([])
			stable_trajectory_times.append([])
			counter = 0;
			for i in range(len_stable_clusters):
				if len(self.stable_cluster[i])>0:
					stable_trajectories[-1].append(self.stable_cluster[i].pop())
					stable_trajectory_times[-1].append(self.stable_cluster_times[i].pop())
				else:
					counter+=1;
			print(counter)
			if counter == len_stable_clusters:
				sw=False
		self.stable_trajectories = stable_trajectories
		self.stable_trajectory_times = stable_trajectory_times
		print('Number of stabl trejectories found: ', len(self.stable_trajectories))

	def eval_SolidNeighbours(self):
		pass
	def eval_Q6Q6(self):
		pass
	def acf(self,data):
		pass

	#def find_induction_time():
	#	pass

	def fit_transition():
		pass

	def plot_measurement(self,path=""): #function for plotting. Should be done individually
		fig, axe = plt.subplots(2,2, figsize = (15,10), sharex=True)

		axe[0,0].set_title("solids from {0} points".format(len(self.solids)))
		axe[0,0].scatter(self.wdata_times-self.quench_time,self.solids, s=1, label='overall solids')
		axe[0,0].scatter(self.wdata_times-self.quench_time,self.LargestCluster, s=1, label='largest cluster')
		axe[0,0].axvline(self.tau+self.quench_time-self.quench_time)
		axe[0,0].legend()


		axe[1,0].set_title("msd")
		axe[1,0].scatter(self.wdata_times-self.quench_time, self.msd_mean,s=1)
		axe[1,0].set_ylim(0,max(self.msd_mean))

	#	axe[1,0].set_title('Radius of Gyration, missing {0}'.format(len(self.RoG)-self.len_wdata))
	#	axe[1,0].set_ylabel('Radius')
	#	axe[1,0].set_xlabel('Time')
	#	axe[1,0].scatter(self.wdata_times[:len(self.RoG)],self.RoG,label= 'Radius of Gyration')
	#	axe[1,0].scatter(self.wdata_times[:len(self.RoG)],self.b,label= 'Asphericity')
	#	axe[1,0].scatter(self.wdata_times[:len(self.RoG)],self.c,label= 'acylindrictiy')
	#	axe[1,0].axhline(self.box_length/2,label='Box length/2', linestyle = 'dashed')
	#	axe[1,0].legend()
	#	#print(len(self.RoG),self.len_wdata)


#		for i in range(10,len(self.stable_cluster)):
#			axe[0,1].scatter(self.stable_cluster_times[i]-self.quench_time,self.stable_cluster[i])
		for i in range(len(self.stable_trajectory_times)):
			axe[0,1].plot(self.stable_trajectory_times[i],np.array(self.stable_trajectories[i])**(1/3))
		#axe[0,1].set_yscale('log')

		#def forward(x):
		#	return x**(1/2)


		#def inverse(x):
		#	return x**2

		#axe[0,1].set_yscale('function')	
		#axe[0,1].yaxis.set_major_locator(FixedLocator(np.arange(0, 100, 10)**(1/3)))
		#axe[0,1].yaxis.set_major_locator(FixedLocator(np.arange(0, 100, 10)))

		temp_pnt=np.log(self.pnt)
		temp_pnt[temp_pnt==temp_pnt[0,-1]]=-1

		im = axe[1,1].pcolormesh(self.time_edges,np.linspace(-0.5,self.maxc+0.5,self.maxc+1),temp_pnt)
		#fig.colorbar(im)
		axe[1,1].set_ylabel('N')
		axe[1,1].set_xlabel('t')
		axe[1,1].set_title('p(n,t) RAW')
		inset=fig.add_axes([0.97,0.07,0.02,0.4], 'off')
		inset.set_xticks([])
		inset.set_yticks([])
		plt.colorbar(im, ax = inset, fraction = 1)
		for i in range(len(self.reset_times)):
			#for j in range(4):
			#	axe[j%2,(int(j/2)%2)].axvline(self.reset_times[i], linestyle='dashed', color='grey',alpha=0.8, linewidth=2)	
			print(i,self.reset_times[i])
			axe[1,1].axvline(self.reset_times[i]-self.quench_time, linestyle='dashed', color='grey',alpha=0.8, linewidth=2)
			axe[0,1].axvline(self.reset_times[i]-self.quench_time, linestyle='dashed', color='grey',alpha=0.8, linewidth=2)
			axe[1,0].axvline(self.reset_times[i]-self.quench_time, linestyle='dashed', color='grey',alpha=0.8, linewidth=2)
			axe[0,0].axvline(self.reset_times[i]-self.quench_time, linestyle='dashed', color='grey',alpha=0.8, linewidth=2)

		fig.tight_layout()
		
		
		plt.savefig(path+"Overview"+self.infostring+".png") 
		plt.close('all')
		#plt.show()

	def close(self):
		self.read_file.close()

def extract_parameters(filename):
	split_string = filename.split("_")
	N = int(split_string[1])
	rho_i = float(split_string[2])*1e-3
	rho_f = float(split_string[3])*1e-3
	eq = int(split_string[4][2:-1])*1e3
	pr = int(split_string[5][2:-1])*1e3
	seed = int(split_string[6][:-3])
	return(N,rho_i,rho_f,eq,pr,seed)


def show_broken_files(path,write_out=True):
	fine=0
	filenames=np.sort(os.listdir(path))

	if write_out:
		log_file=open("missed.txt","w+")

	for i in range(len(filenames)):
		try:
			read_file=h5py.File(path+filenames[i],"r")
			read_file.close()
			fine+=1
		except:
			print(filenames[i])
			if write_out:
#				log_file.write("{0} {1} {2} {3} {4} {5}\n".format(N, rho_i, rho_f, eq, pr, seed))
				log_file.write("{0} {1} {2} {3} {4} {5}\n".format(extract_parameters(filenames[i])))
	
	print("Fine are {0}/{1} = {2:4.2f}".format(fine,len(filenames),fine/len(filenames)))
	if write_out:
		log_file.close()



class h5md_DataSeries:
	#contains and evaluates on a larger scale a series of datasets
	def __init__(self,path, plot_parameters):
		self.path=path

		self.filenames = np.sort(os.listdir(self.path))
		self.num=len(self.filenames)
		self.title_name=path[-15:] #quite a rough naming


		#hard coded quantites (so far, same as for dataset):
		#stop_time must not exceed any simulation time
		#step_number should be choosen of the same order as the actual measurement number in the time frame
		self.step_number = plot_parameters[0]
		self.stop_time = plot_parameters[1]
		self.maxc = plot_parameters[2]
		self.plot_parameters = plot_parameters

		self.equi_time=np.linspace(0,self.stop_time,self.step_number)

		self.time_edges = np.zeros(self.step_number+1)		
		self.time_edges[0]=(1.5*self.equi_time[0]-0.5*self.equi_time[1])
		for i in range(self.step_number-1):
			self.time_edges[i+1]=((self.equi_time[i+1]+self.equi_time[i])/2)
		self.time_edges[self.step_number]=(1.5*self.equi_time[-1]-0.5*self.equi_time[-2])


	def collect_largest_cluster(self, plot=True, number=-1):
		print("Collecting largest cluster data from"+self.path)
		if number == -1 :
			number=self.num
		else:
			number = min(self.num, number)

		self.tau_lc=np.zeros(number)
		self.quench_times=np.zeros(number)

		fig2 = plt.figure()
		for i in range(number):
			N,rho_i,rho_f,eq,pr,seed=extract_parameters(self.filenames[i])
			datenset = h5md_dataset(self.path,N,rho_i,rho_f,eq,pr,seed, self.plot_parameters)
			datenset.load_watchtable()
			datenset.eval_LargestCluster()
			self.tau_lc[i]=datenset.tau
			self.quench_times[i]=datenset.quench_time
			datenset.close()

			if plot:
				datenset.plot_LargestCluster()
				plt.plot(datenset.wdata_times[10:],datenset.LargestCluster[10:], linewidth=0.1, c = "C{0}".format(i%10), alpha = 0.6, figure = fig2)
				plt.axvline(datenset.quench_time, c = 'k', figure = fig2)
				plt.axvline(datenset.tau+datenset.quench_time, c = "C{0}".format(i%10), figure = fig2)
				plt.ylim(0,N+1000)
		if plot: 
			fig2.savefig("lc_plots/Largest_cluster_over_time"+self.title_name.replace('/','_')+".png", dpi = 400)

			fig3 = plt.figure()
			plt.hist(self.tau_lc, bins = 50, figure = fig3)
			fig3.savefig("lc_plots/Tau_distribution"+self.title_name.replace('/','_')+".png")			
			plt.close(fig2)
			plt.close(fig3)

		plt.close('all')		


	def collect_solid_transition(self, plot=True):
		pass

	
	def collect_pnt(self, plot=True, number=-1):
		#collect_cluster_distribution is already taken by higher layer..
		if number == -1 :
			number=self.num
		else:
			number = min(self.num, number)

		self.pnt=np.zeros(shape=(self.maxc, self.step_number))
		for i in range(number):
			N,rho_i,rho_f,eq,pr,seed=extract_parameters(self.filenames[i])
			datenset = h5md_dataset(self.path,N,rho_i,rho_f,eq,pr,seed, self.plot_parameters)
			
			if False:
				print("Load watchtable")
				datenset.load_watchtable()
				print("eval solid frac")
				datenset.eval_SolidFraction()
				print("eval largest clsuter")
				datenset.eval_LargestCluster()
				print("eval cluster distribution")
				datenset.eval_ClusterDistribution()
				print("eval MSD")
				datenset.eval_MSD()
				datenset.eval_ToG()
				print("plot all")
				datenset.plot_measurement("4k_overviews/singles/")
			else:
				datenset.eval_ClusterDistribution()

			self.pnt = self.pnt + datenset.pnt

		if plot:
			temp_pnt = np.log(self.pnt)
			temp_pnt[temp_pnt==temp_pnt[0,0]] = -1
			im=plt.pcolormesh(self.time_edges,np.linspace(-0.5,self.maxc+0.5,self.maxc+1),temp_pnt)
			plt.colorbar(im)
			plt.ylabel('t')
			plt.xlabel('N')
			plt.title('p(n,t) RAW')
			plt.savefig("lc_plots/ClusterDistribution"+self.title_name.replace('/','_')+".png")
			plt.close('all')


class h5md_DataRun:
	#brings together dataseries to compare
	def __init__(self, paths, pos1, pos2, names, plot_parameters, multi_threading=-1, run_name = "", evaluation_number = -1):
		self.num = len(paths) # all lists should have the same length
		self.series=[]
		for i in range(self.num):
			self.series.append(h5md_DataSeries(paths[i], plot_parameters))
		self.pos1=np.array(pos1)
		self.pos2=np.array(pos2)
		self.names=names
	
		self.bins=100
		self.min_tau=1e8
		self.max_tau=0
		self.evaluation_number = evaluation_number
		self.multi_threading = multi_threading
		self.run_name = run_name

		#hard coded quantites (so far, same as for dataset):
		#stop_time must not exceed any simulation time
		#step_number should be choosen of the same order as the actual measurement number in the time frame
		self.step_number = plot_parameters[0]
		self.stop_time = plot_parameters[1]
		self.maxc = plot_parameters[2]

		self.equi_time=np.linspace(0,self.stop_time,self.step_number)

		self.time_edges = np.zeros(self.step_number+1)		
		self.time_edges[0]=(1.5*self.equi_time[0]-0.5*self.equi_time[1])
		for i in range(self.step_number-1):
			self.time_edges[i+1]=((self.equi_time[i+1]+self.equi_time[i])/2)
		self.time_edges[self.step_number]=(1.5*self.equi_time[-1]-0.5*self.equi_time[-2])



	def collect_transition_times_A(self):
	# Iterates through the dataseries' collecting transition times from each
		for i in range(self.num):
			self.series[i].collect_largest_cluster(plot=True, number=self.evaluation_number)
			if (min(self.series[i].tau_lc)<self.min_tau):
				self.min_tau = min(self.series[i].tau_lc)
			if (max(self.series[i].tau_lc)>self.max_tau):
				self.max_tau = max(self.series[i].tau_lc)

	def multi_thread_func_lc(self, dataset, num):
	# Helper function for collect_transition_times_B()
		dataset.collect_largest_cluster(plot = True, number=num)
		return(dataset)

	def collect_transition_times_B(self):
	# Version of collect_transition_times_A(), using multiple threads
		import multiprocessing as mp
		pool = mp.Pool(processes= self.multi_threading)
		jobs = []
		for i in range(self.num):
			jobs.append(pool.apply_async(self.multi_thread_func_lc, args = (self.series[i], self.evaluation_number)))

		pool.close()
		pool.join()

		for i in range(self.num):
			self.series[i] = jobs[i].get()
			if (min(self.series[i].tau_lc)<self.min_tau):
				self.min_tau = min(self.series[i].tau_lc)
			if (max(self.series[i].tau_lc)>self.max_tau):
				self.max_tau = max(self.series[i].tau_lc)


	def collect_transition_times(self):
	# MAKRO for single or multithreaded evaluation
		if (self.multi_threading == -1):
			self.collect_transition_times_A()
		else:
			self.collect_transition_times_B()

	def plot_transition_times(self):
		#note : not usable for single row or column subplots.
		fig, axe = plt.subplots(max(self.pos1)+1,max(self.pos2)+1,sharex = True)
		binning=np.linspace(self.min_tau,self.max_tau,self.bins)
		for i in range(self.num):
			print(self.series[i].tau_lc)
			axe[self.pos1[i],self.pos2[i]].hist(self.series[i].tau_lc,bins=binning)
			axe[self.pos1[i],self.pos2[i]].set_title(self.names[i])
			
		fig.tight_layout()
		fig.savefig("Overall_transition_times"+self.run_name+".png")
		plt.close('all')

	def multi_thread_func_cd(self, dataset, num):
	# Helper function for collect_transition_times_B()
		dataset.collect_pnt(plot = True, number=num)
		return(dataset)

	def collect_cluster_distribution_A(self):
		self.pnt=np.zeros(shape=(self.maxc,self.step_number))
		for i in range(self.num):
			self.series[i].collect_pnt(plot = True, number=self.evaluation_number)
			self.pnt = self.pnt + self.series[i].pnt

	def collect_cluster_distribution_B(self):
	# Version of collect_transition_times_A(), using multiple threads
		import multiprocessing as mp
		pool = mp.Pool(processes= self.multi_threading)
		jobs = []
		for i in range(self.num):
			jobs.append(pool.apply_async(self.multi_thread_func_cd, args = (self.series[i], self.evaluation_number)))

		pool.close()
		pool.join()

		self.pnt=np.zeros(shape=(self.num, self.maxc, self.step_number))
		for i in range(self.num):
			self.series[i] = jobs[i].get()
			temp_pnt = self.series[i].pnt
			print(np.max(temp_pnt), np.min(temp_pnt)*1e10)
			temp_pnt = np.log(temp_pnt)
			temp_pnt[temp_pnt==temp_pnt[-1,0]]=-1
			self.pnt[i,:,:]=temp_pnt

	def collect_cluster_distribution(self):
	# MAKRO for single or multithreaded evaluation
		if (self.multi_threading == -1):
			self.collect_cluster_distribution_A()
		else:
			self.collect_cluster_distribution_B()


	
	def plot_cluster_distribution(self):
		fig, axe = plt.subplots(max(self.pos1)+1,max(self.pos2)+1,sharex = True)
		for i in range(self.num):
			im=axe[self.pos1[i],self.pos2[i]].pcolormesh(self.time_edges,np.linspace(-0.5,self.maxc+0.5,self.maxc+1),self.pnt[i,:,:])
			axe[self.pos1[i],self.pos2[i]].set_title(self.names[i])
			fig.colorbar(im, ax=axe[self.pos1[i],self.pos2[i]])
		fig.tight_layout()
		fig.savefig("Overall_ClusterDistribution"+self.run_name+".png", dpi = 600)

		plt.close('all')




		#im=plt.pcolormesh(self.time_edges,np.linspace(-0.5,self.maxc+0.5,self.maxc+1),self.pnt)
		#plt.colorbar(im)
		#plt.ylabel('t')
		#plt.xlabel('N')
		#plt.title('p(n,t) RAW')
		#plt.savefig("lc_plots/ClusterDistribution"+self.title_name.replace('/','_')+".png")
		#plt.close('all')







if __name__ == "__main__":
	#series=h5md_DataSeries()
	print("Es gibt hier nichts zu sehen.")
	pass
