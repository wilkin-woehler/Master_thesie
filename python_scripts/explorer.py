#author : wilkin woehler
#date : 11.05.2020
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import h5py
import time as time

eval_hist=False
dump=False

#read in h5(md) file from which data is supposed to be read
N=32000
rho_inital=0.45
rho=0.550 #rho_final in production phase
steps_eq=0
steps_pr=1500
seed=1


#N=16384
#rho_inital=0.45
#rho=0.54 #rho_final in production phase
#steps_eq=50000
#steps_pr=50000
#seed=720








def dump_file(positions,step,box_bound,file,typ=1):
	#print(positions)
	#very basic dump-file creator only saving position and time step
	# var: positions, array containing (N,dim) entries for N particles dim
	# var: step
	# var: box_bound
	# var: file
	if type(typ)==type(1):
		typ=np.ones(len(positions))
		#print('something')
	#file=open(file_path,'r+')
	file.read()
	file.write('ITEM: TIMESTEP\n')
	file.write(str(int(step))+'\n')
	file.write('ITEM: NUMBER OF ATOMS\n')
	file.write(str(len(positions))+'\n')
	file.write('ITEM: BOX BOUNDS xx yy zz\n')

	lx=box_bound[0][1]-box_bound[0][0]
	ly=box_bound[1][1]-box_bound[1][0]
	lz=box_bound[2][1]-box_bound[2][0]

	for i in range(len(box_bound)):
		file.write('{0} {1}\n'.format(box_bound[i][0],box_bound[i][1]))

	file.write('ITEM: ATOMS id type radius x y z\n')
	for i in range(len(positions)):
		file.write('{0} {1} {5} {2} {3} {4}\n'.format(i+1,int(typ[i]),positions[i,0],positions[i,1],positions[i,2],0.5))


def read_out_hist(name,read_file):
	data=read_file['observables/'+name][()]
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

def eval_set(N,rho_inital,rho,steps_eq,steps_pr,seed):

	infostring="_{0}_{1}_{2}_eq{3}k_pr{4}k_{5}".format(N,int(rho_inital*1000),int(rho*1000),int(steps_eq*1e-3),int(steps_pr*1e-3),seed)

	filename = "../outfiles/NEMO/16_k_december/eq_steps_500/Out"+infostring+".h5"
	#filename = "../outfiles/Out"+infostring+".h5"
	print('read in file: '+filename)
	read_file=h5py.File(filename, "r")


	#show all folders and subfolders, from which a specific one can be picked
	def print_names(name):
	    print(name)
	#read_file.visit(print_names)

	#show all folders and subfolders, with more information
	def print_names(name,object_name):
	    print(name)
	    print(object_name)
	#read_file.visititems(print_names)


	#data = read_file['particles/pr_all_parts/position/value'][()]
	#step_data = read_file['particles/pr_all_parts/position/step'][()]
	#time_data = read_file['particles/pr_all_parts/position/time'][()]
	#print(step_data[-1])
	#print(time_data[-1])

	#print(time_data)
	#print(len(data))
	#print(data)
	#for i in range(len(data)):
	#	print(data[i])





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
		for j in range(len(data[i])):
			pass
			#print(data[i][j])
		steps[i]=data[i][0]
		times[i]=data[i][1]
		solids[i]=data[i][2]
		clus_size[i]=data[i][3]
		msd_mean[i]=data[i][4]
		s_msd_mean[i]=data[i][5]
		vel_2[i]=data[i][6]





	box_length=(6*rho/np.pi/N)**(-1/3)


#	ToG = read_file['observables/TensorOfGyration'][()]

#	RoG=np.zeros(len(ToG))
#	b=np.zeros(len(ToG))
#	c=np.zeros(len(ToG))
#	k2=np.zeros(len(ToG))
#	l=np.zeros(shape=(len(ToG),3))

#	for i in range(len(ToG)):
#		#print(ToG[i])
#		xx=ToG[i][0]
#		yy=ToG[i][1]
#		zz=ToG[i][2]
#		xy=ToG[i][3]
#		xz=ToG[i][4]
#		yz=ToG[i][5]
#		M=np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])
#		#print(M)
#		#eig_val, eig_vec = np.linalg.eig(M)
#
#		eig_val, eig_vec = np.linalg.eig(M)
#		try:
#			M_diag= np.matmul(np.matmul(np.linalg.inv(eig_vec),M),eig_vec)
#		except:
#			M_diag = M 
#			print('Matrix was not diagonalizabale : ', M)
#		l=np.sort(M_diag.diagonal()) 
#		#print(l)
#		RoG[i] = np.sqrt(l[0]+l[1]+l[2])
#		#print(RoG[i])
#		b[i]   = l[2]-1/2*(l[0]+l[1])
#		c[i]   = l[1]-l[0]
#		if RoG[i] > 1e-10:
#			k2[i]  = 3/2*(l[0]**2+l[1]**2+l[2]**2)/(l[0]+l[1]+l[2])**2-1/2
#		else:
#			k2[i]  = 0

	#Q6R_raw = read_file['observables/Q6_real'][()]
	#Q6C_raw = read_file['observables/Q6_complex'][()]
	#Q6R = np.zeros(shape=(len(Q6R_raw),7))
	#Q6C = np.zeros(shape=(len(Q6C_raw),7))
	#Q6=np.zeros(shape=(len(Q6R_raw),7))
	#for i in range(len(Q6R_raw)):
	#	for j in range(7):
	#		Q6R[i][j] = Q6R_raw[i][j]
	#		Q6C[i][j] = Q6C_raw[i][j]
	#		Q6[i,j]= np.sqrt( Q6C[i,j]**2 + Q6R[i,j]**2 ) 



#exit()


	#file.close() would be better but 


#opening a file making it ready to be written on, in a not very elegant fashion... and the dumping the position time_series

	if dump:
		file_path='test.dump'
		try:
			os.remove(file_path)
		except:
			print('no file yet there but is created')

		new_file=open(file_path,'w+')
		new_file.close()


		l0=box_length
		sigma=1

		lx=l0*sigma
		ly=l0*sigma
		lz=l0*sigma

		#L=[[-lx/2,lx/2],[-ly/2,ly/2],[-lz/2,lz/2]] #box length maybe read it also from h5md file
		L=[[0,lx],[0,ly],[0,lz]]

		dump_time_1=time.time()
		file=open(file_path,'r+')
		#print(len(data1),len(data1[0]),len(data1[0][0]))
		print(data[2][0])
		for i in range(2,len(data)):
			#print(data[i][0])
			dump_file(data[i][0],step_data[i],L,file,typ=np.arange(len(data[i][0]))+1)
		file.close()
		dump_time_2=time.time()

		print('Time to write down all dump files = {0:4.2f}'.format(dump_time_2-dump_time_1))






#exit()


	if eval_hist:

		g_r_t=read_file['observables/RDF_hist'][()]

		g_r=[]
		r=[]
		step_check=-1
		for i in range(len(g_r_t)):
			if (step_check!=g_r_t[i][0]):
				g_r.append([])
				r.append([])
				step_check=g_r_t[i][0]

			r[-1].append(g_r_t[i][2])
			g_r[-1].append(g_r_t[i][3])


		msd_data=read_file['observables/msd_hist'][()]
		msd_count=[]
		msd=[]
		step_check=-1
		for i in range(len(msd_data)):
			if (step_check!=msd_data[i][0]):
				msd_count.append([])
				msd.append([])
				step_check=msd_data[i][0]
			msd[-1].append(msd_data[i][2])
			msd_count[-1].append(msd_data[i][3])

		q6q6_data=read_file['observables/q6q6_hist'][()]
		q6q6_count=[]
		q6q6=[]
		step_check=-1
		for i in range(len(q6q6_data)):
			if (step_check!=q6q6_data[i][0]):
				q6q6_count.append([])
				q6q6.append([])
				step_check=q6q6_data[i][0]
			q6q6[-1].append(q6q6_data[i][2])
			q6q6_count[-1].append(q6q6_data[i][3])

		sol_neigh_data=read_file['observables/sol_neigh_hist'][()]
		sol_neigh_count=[]
		sol_neigh=[]
		step_check=-1
		for i in range(len(sol_neigh_data)):
			if (step_check!=sol_neigh_data[i][0]):
				sol_neigh_count.append([])
				sol_neigh.append([])
				step_check=sol_neigh_data[i][0]
			sol_neigh[-1].append(sol_neigh_data[i][2])
			sol_neigh_count[-1].append(sol_neigh_data[i][3])

		vel_data=read_file['observables/Velocity_hist'][()]
		vel_count=[]
		vel=[]
		step_check=-1
		for i in range(len(vel_data)):
			if (step_check!=vel_data[i][0]):
				vel_count.append([])
				vel.append([])
				step_check=vel_data[i][0]
			vel[-1].append(vel_data[i][2])
			vel_count[-1].append(vel_data[i][3])
		vel=np.array(vel)
		vel_count=np.array(vel_count)



		#Cluster_size, Cluster_count = read_out_hist('Cluster_hist', read_file)


	#print("number of snapshots : {0}".format(len(q6q6)))
	print("number of evaluations : {0}".format(len(msd_mean)))
	#print("steps", steps)

	if eval_hist:

		print("")
		l1=len(read_file['observables/Cluster_Distribution/value'])
		#if l1 != 0:

		maxc=200
		pnt = np.zeros(shape=(maxc,l1))

		for i in range(l1):
			#print(read_file['observables/Cluster_Distribution/value'][i])
			for j in range(len(read_file['observables/Cluster_Distribution/value'][i][0])):
				try:
					k=read_file['observables/Cluster_Distribution/value'][i][0][j][0]
					if k < maxc:
						pnt[k,i]+=read_file['observables/Cluster_Distribution/value'][i][0][j][1]
				except:
					pass;

		time_data = read_file['particles/pr_all_parts/position/time'][()]
		time_edges = []
		reset_times = read_file['reset_sim/particles/position/time'][()]
		

		time_edges.append(1.5*time_data[0]-0.5*time_data[1])
		for i in range(len(time_data)-1):
			time_edges.append((time_data[i+1]+time_data[i])/2)
		time_edges.append(1.5*time_data[-1]-0.5*time_data[-2])

#		print(time_data)
		#exit
		pnt[0,:]=0

		pnt = np.log(pnt)
		pnt[pnt==pnt[0,0]]=-1



		fig, axe = plt.subplots(2,3, figsize = (15,8))
		axe[0,0].pcolormesh(np.linspace(-0.5,maxc+0.5,maxc+1),time_edges,pnt.T)
		#axe[0,0].colorbar(im)
		axe[0,0].set_ylabel('t')
		axe[0,0].set_xlabel('N')
		axe[0,0].set_title('p(n,t) RAW')
		for i in range(len(reset_times)):
			axe[0,0].axhline(reset_times[i], linestyle='dashed', color='grey',alpha=0.8, linewidth=2)
			axe[0,1].axvline(reset_times[i], linestyle='dashed', color='grey',alpha=0.8, linewidth=2)
			axe[0,2].axvline(reset_times[i], linestyle='dashed', color='grey',alpha=0.8, linewidth=2)
			axe[1,2].axvline(reset_times[i], linestyle='dashed', color='grey',alpha=0.8, linewidth=2)

		#		set_title("g(r)")
		#for i in range(len(g_r)):
		#	axe[0,0].plot(r[i],g_r[i])

		axe[0,1].set_title("msd")
		for i in range(len(msd)):
			#axe[0,1].plot(np.array(msd_count[i])*100000000,np.array(msd[i])**2)
			axe[0,1].scatter(times, msd_mean,s=1)
			axe[0,1].set_ylim(0,max(msd_mean))


		axe[1,0].set_title("(q6q6)^sol_neigh_threshold")
		for i in range(len(sol_neigh)):
			axe[1,0].plot(q6q6[i],np.array(q6q6_count[i])**8)
		axe[1,0].axvline(0.6)

		axe[1,1].set_title("sol_neigh")
		for i in range(len(sol_neigh)):
			axe[1,1].plot(sol_neigh[i],sol_neigh_count[i])
		axe[1,1].axvline(8)

		axe[1,2].set_title("solids from {0} points".format(len(solids)))
		#axe[1,2].scatter(times,solids, s=1, label='overall solids')
		axe[1,2].scatter(times,clus_size, s=1, label='largest cluster')
		axe[1,2].legend()
		#axe[1,2].set_ylim(0,N)

		#axe[0,2].set_title('velocity')
		#for i in range(len(vel)):
		#	axe[0,2].scatter(vel_count[i]*np.max(times)*10, vel[i]**2)	
		#axe[0,2].scatter(times, vel_2)


		#axe[0,2].set_title('clusters')
		#for i in range(len(Cluster_size)):
		#	axe[0,2].scatter(Cluster_size[i], Cluster_count[i]*N, s = 1)
			#print(Cluster_count[i][:10])
		#axe[0,2].set_ylim(0,1e-4)
		#axe[0,2].set_xlim(0,30)

		axe[0,2].set_title('Radius of Gyration')
		axe[0,2].set_ylabel('Radius')
		axe[0,2].set_xlabel('Time')
		#axe[0,2].scatter(time_data,RoG,label= 'Radius of Gyration')
		#axe[0,2].scatter(time_data,b,label= 'Asphericity')
		#axe[0,2].scatter(time_data,c,label= 'acylindrictiy')
		axe[0,2].axhline(box_length/2,label='Box length/2', linestyle = 'dashed')
		#axe[0,2].scatter(time_data,k2*box_length/2,label = 'realative shape anisotropy * {0:2.0f}'.format(box_length/2))
		axe[0,2].legend()

		#axe[1,3].set_title('Q6')
		#axe[1,3].set_ylabel('<q6l>')
		#axe[1,3].set_xlabel('Time')
		#for i in range(7):
		#	axe[1,3].scatter(time_data,Q6[:,i],c='C{0}'.format(i),label='Q6{0}'.format(i))
		#	axe[1,3].scatter(time_data,Q6R[:,i],c='k'.format(i),label='Q6{0}'.format(i))
		#	axe[1,3].scatter(time_data,Q6C[:,i],c='orange'.format(i),label='Q6{0}'.format(i))
		#axe[1,3].set_ylim(np.min(Q6),np.max(Q6))
		#axe[1,3].legend()


		fig.tight_layout()
		fig.savefig("Overview_for_Out"+infostring+".png")

	else:
		
		print("")
		fig, axe = plt.subplots(2,2, figsize = (15,10))

		axe[0,1].set_title("msd")
		axe[0,1].scatter(times, msd_mean,s=1)
		axe[0,1].set_ylim(0,max(msd_mean))

		
		axe[0,0].set_title("solids from {0} points".format(len(solids)))
		axe[1,0].plot(times,solids, label='overall solids')
		axe[0,0].plot(times[np.abs(clus_size)<20000],solids[np.abs(clus_size)<20000], label='overall solids')
		axe[0,0].plot(times[np.abs(clus_size)<20000],clus_size[np.abs(clus_size)<20000], label='largest cluster')
		axe[1,0].plot(times,clus_size, label='largest cluster')
		axe[0,0].legend()
		print(clus_size)
		#axe[1,2].set_ylim(0,N)

		#axe[0,2].set_title('velocity')
		#for i in range(len(vel)):
		#	axe[0,2].scatter(vel_count[i]*np.max(times)*10, vel[i]**2)	
		#axe[0,2].scatter(times, vel_2)

		#print(len(time_data),len(RoG))
#		axe[1,0].set_title('Radius of Gyration')
#		axe[1,0].set_ylabel('Radius')
#		axe[1,0].set_xlabel('Time')
		#axe[1,0].scatter(times,RoG,label= 'Radius of Gyration')
		#axe[1,0].scatter(times,b,label= 'Asphericity')
		#axe[1,0].scatter(times,c,label= 'acylindrictiy')
#		axe[1,0].axhline(box_length/2,label='Box length/2', linestyle = 'dashed')
		#axe[1,0].scatter(times,k2*box_length/2,label = 'realative shape anisotropy * {0:2.0f}'.format(box_length/2))
#		axe[1,0].legend()

		#axe[1,1].set_title('Q6')
		#axe[1,1].set_ylabel('<q6l>')
		#axe[1,1].set_xlabel('Time')

		#axe[0,2].set_title('Q6')
		#axe[0,2].set_ylabel('<q6l>')
		#axe[0,2].set_xlabel('Time')

		#axe[1,2].set_title('Q6')
		#axe[1,2].set_ylabel('<q6l>')
		#axe[1,2].set_xlabel('Time')
		#for i in range(7):
				#axe[1,1].scatter(time_data,Q6[:,i],c='C{0}'.format(i),label='Q6_2{0}'.format(i),s=1)
			#axe[0,2].scatter(time_data,Q6R[:,i],c='C{0}'.format(i),label='Q6_R{0}'.format(i),s=1)
			#axe[1,2].scatter(time_data,Q6C[:,i],c='C{0}'.format(i),label='Q6_C{0}'.format(i),s=1)
			#axe[1,1].plot(time_data,Q6[:,i],c='C{0}'.format(i),label='Q6_2{0}'.format(i),linewidth=1)
			#axe[0,2].plot(time_data,Q6R[:,i],c='C{0}'.format(i),label='Q6_R{0}'.format(i),linewidth=1)
			#axe[1,2].plot(time_data,Q6C[:,i],c='C{0}'.format(i),label='Q6_C{0}'.format(i),linewidth=1)
		#axe[1,1].set_ylim(np.min(Q6),np.max(Q6))
		#axe[1,1].legend()

		#axe[0,2].set_ylim(np.min(Q6R),np.max(Q6R))
		#axe[0,2].legend()

#		axe[1,2].set_ylim(np.min(Q6C),np.max(Q6C))
#		axe[1,2].legend()


#		l1=len(read_file['observables/Cluster_Distribution/value'])
		#if l1 != 0:

#		maxc=200
#		pnt = np.zeros(shape=(maxc,l1))
#		print('number of evalutations to be displayed: ', l1)
		#for i in range(l1):
			#if ((i+1)%int(l1/100)==0):	
			#	print("{0:3.0f} %".format(i/l1*100))
			#print(read_file['observables/Cluster_Distribution/value'][i])
#			for j in range(len(read_file['observables/Cluster_Distribution/value'][i][0])):
#				try:
#					k=read_file['observables/Cluster_Distribution/value'][i][0][j][0]
#					if k < maxc:
#						pnt[k,i]+=read_file['observables/Cluster_Distribution/value'][i][0][j][1]
#				except:
#					pass;


#
#		temp_cd=read_file['observables/Cluster_Distribution/value']
#		print('number of evalutations to be displayed: ', l1)
#		for i in range(l1):
#			#if ((i+1)%int(l1/100)==0):	
#			#	print("{0:3.0f} %".format(i/l1*100))
#			#print(read_file['observables/Cluster_Distribution/value'][i])
#			temp_cd1=temp_cd[i][0]
#			for j in range(len(temp_cd1)):
#				try:
#					k=temp_cd1[j][0]
#					if k < maxc:
#						pnt[k,i]+=temp_cd1[j][1]
#				except:
#					pass;
#
#
#		time_data = read_file['particles/pr_all_parts/position/time'][()]
#		time_edges = []
		reset_times = read_file['reset_sim/particles/position/time'][()]
		

#		time_edges.append(1.5*times[0]-0.5*times[1])
#		for i in range(len(times)-1):
#			time_edges.append((times[i+1]+times[i])/2)
#		time_edges.append(1.5*times[-1]-0.5*times[-2])
#
#		print(time_data)
		#exit
#		pnt[0,:]=0
#		print(np.max(pnt),np.min(pnt))
		#exit()
#		pnt = np.log(pnt)
#		pnt[pnt==pnt[0,0]]=-1


#
#		#fig, axe = plt.subplots(2,3, figsize = (15,8))
#		axe[1,1].pcolormesh(np.linspace(-0.5,maxc+0.5,maxc+1),time_edges,pnt.T)
#		#axe[0,0].colorbar(im)
#		axe[1,1].set_ylabel('t')
#		axe[1,1].set_xlabel('N')
#		axe[1,1].set_title('p(n,t) RAW')
		for i in range(len(reset_times)):
			axe[1,1].axhline(reset_times[i], linestyle='dashed', color='grey',alpha=0.8, linewidth=2)
			axe[0,1].axvline(reset_times[i], linestyle='dashed', color='grey',alpha=0.8, linewidth=2)
			axe[1,0].axvline(reset_times[i], linestyle='dashed', color='grey',alpha=0.8, linewidth=2)
			axe[0,0].axvline(reset_times[i], linestyle='dashed', color='grey',alpha=0.8, linewidth=2)







		fig.tight_layout()
		fig.savefig("Overview_for_Out"+infostring+"_without_histograms.png")
		#plt.show()

	#off=3
	#delta_steps=steps[off+1:-off]-steps[off:-off-1]
	#fig,axe = plt.subplots(2)
	#axe[0].scatter(times[off:-off],solids[off:-off]/N,  s = 1)
	#axe[0].set_ylabel('Solids fration')
	#axe[0].set_ylim(0,1)
	#axe[1].scatter(times[off+1:-off][delta_steps>0],delta_steps[delta_steps>0], s = 1)
	#axe[1].set_ylabel('steps step')
	##axe[2].scatter(times[off+1:-off],times[off+1:-off]-times[off:-off-1])
	##axe[2].set_ylim(0,max(times[off+1:-off]-times[off:-off-1]))
	##axe[2].set_ylabel('time step')
	#axe[1].axhline(0)

	#fig.savefig("Transition_for_Out"+infostring+".png")
	read_file.close()
	plt.close('all')
	return

#Out_1048576_450_526_eq1k_pr1k_10.h5
#Out_1048576_450_540_eq0k_pr10k_1.h5
#Out_32000_450_545_eq0k_pr3k_2.h5
#Out_32000_450_538_eq0k_pr1000k_1.h5

numbers=np.linspace(16,16,1,dtype=int)**3*4
#numbers[4]=64**3*4
print("N from ", numbers, "are going to be evaluated ")
rho_inital=0.450
rho=np.linspace(0.540,0.540,1) #rho_final in production phase
steps_eq=500
steps_pr=200000

n_stat=1
start_seed=300
seeds=np.arange(start_seed,start_seed+n_stat)
print("seeds from ", seeds, "are going to be evaluated ")


missed=0
for i in range(len(numbers)):
	for j in range(len(seeds)):
		for k in range(len(rho)):
			try:
			#if True:
				eval_set(numbers[i],rho_inital,rho[k],steps_eq,steps_pr,seeds[j])
			except:
				missed+=1
				
print("Missed {0} sets.".format(missed))		
