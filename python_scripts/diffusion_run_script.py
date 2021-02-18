import subprocess as sub
import numpy as np
import time
import os
import matplotlib.pyplot as plt
from datetime import date
import h5py

def remove_file(file_name):
	try:
		os.remove(file_name)
	except:
		print("no "+file_name+" yet")

#generate infile(s)
def gen_in_out_files(N_part,phi_init,phi_final,step_array):
        #try cleaning, ovito_files, out_files, output, init_files and in_files

	steps_eq,steps_pr,sweep_eq,sweep_co,sweep_pr = step_array[0],step_array[1],step_array[2],step_array[3],step_array[4]

	in_file_name= "in_files/InFile_{0}_{1}_{2}_{3}k".format(N_part,int(phi_init*1000),int(phi_final*1000),int(steps_pr*1e-3))
	terminal_out_file_name = "output/terminal_out_{0}_{1}_{2}_{3}k.txt".format(N_part,int(phi_init*1000),int(phi_final*1000),int(steps_pr*1e-3))
	ovito_file_name ="ovito_files/ovito_{0}_{1}.dump".format(N_part,int(phi_final*1000))
	out_file_name = "outfiles/Out_{0}_{1}.h5".format(N_part,int(phi_final*1000))
	init_file_name = "initfiles/Initfile_{0}_{1}.h5".format(N_part,int(phi_final*1000))

	remove_file(in_file_name)
	remove_file(terminal_out_file_name)
	remove_file(ovito_file_name)
	remove_file(out_file_name)
	remove_file(init_file_name)	


	in_file = open(in_file_name,"w+")
	in_file.write("Seed 4\n")
	in_file.write("RNGInc 0\n")
	in_file.write("Dimension 3\n")
	in_file.write("Temperature 1.00\n")
	in_file.write("Number "+str(N_part)+"\n")
	in_file.write("DensityFinal "+str(phi_final)+"\n")
	in_file.write("DensityInital "+str(phi_init)+"\n")
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
	
	return(in_file_name,terminal_out_file_name)


#run process(es)
def run_sim(phi_init,phi_final,step_array):
	in_file_name,terminal_out_file_name = gen_in_out_files(N_part,phi_init,phi_final,step_array)
	print("running a simulation for: ",in_file_name)
	simulation = sub.Popen(["./bin/arcturus.BasicEDMD.x", "-i", "init", "-f", in_file_name, "-o", "Out.h5"], stdout = sub.PIPE, text=True)

	return(in_file_name,terminal_out_file_name,simulation)

def write_out_sim(in_file_name,terminal_out_file_name,simulation):
	out,err = simulation.communicate()
	out_file=open(terminal_out_file_name,"w+")
	out_file.write(out)
	out_file.close()






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

def D_E(rho,g_d):
	m=1
	k_bT=1
	d=1
	return(3/8/rho/d**2/g_d*(k_bT/np.pi/m)**(1/2))

def evaluate_sim(phi_sim,step_array):
	steps_eq,steps_pr,sweep_eq,sweep_co,sweep_pr = step_array[0],step_array[1],step_array[2],step_array[3],step_array[4]
	filename = "outfiles/Out_{0}_{1}.h5".format(N_part,int(phi_sim*1000))
	print('read in file: '+filename)
	read_file=h5py.File(filename, "r")

	#show all folders and subfolders, from which a specific one can be picked
	def print_names(name):
	    print(name)
	#read_file.visit(print_names)

	#show all folders and subfolders, from which a specific one can be picked
	def print_names(name,object_name):
	    print(name)
	    print(object_name)
	#read_file.visititems(print_names)

	
	g_r,r,s_g_r = read_out_hist('RDF_hist', read_file)
	msd_count,msd, s_msd_count = read_out_hist('msd_hist', read_file)
	q6q6_count, q6q6, s_q6q6_count= read_out_hist('q6q6_hist', read_file)
	sol_neigh_count, sol_neigh, s_sol_neigh_count = read_out_hist('sol_neigh_hist', read_file)
	vel_count, vel, s_vel_count = read_out_hist('Velocity_hist', read_file)
	vel=np.array(vel)
	vel_count=np.array(vel_count)

	#obtain step numbers particle numbers and date
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
	msd_mean=np.zeros(len(data))
	s_msd_mean=np.zeros(len(data))
	vel_2=np.zeros(len(data))

	for i in range(len(data)):
		steps[i]=data[i][0]
		times[i]=data[i][1]
		solids[i]=data[i][2]
		msd_mean[i]=data[i][4]
		s_msd_mean[i]=data[i][5]
		vel_2[i]=data[i][6]


	#V_spheres=N*4/3*np.pi*1/2**3
	#V_box=V_spheres/rho
	#box_len=V_box**(1/3)
	#print("box length is : ", box_len)


	#fdate = date.today().strftime('%d_%m_%Y')


	#do a linear regression to the msd_mean data
	#print(sweep_pr)
	N_subdiffusive=3*sweep_pr # number of steps to not be taken into account for diffusive behaviour
	N_start = (eq_steps + co_steps + N_subdiffusive)
	#print(steps[steps>N_start])
	print(msd_mean[steps>N_start])
	#print(s_msd_mean[steps>N_start])
	a,ea,b,eb,chiq,corr = lineare_regression(times[steps>N_start],msd_mean[steps>N_start],s_msd_mean[steps>N_start])

	print("number of snapshots : {0}".format(len(q6q6)))
	print("number of evaluations : {0}".format(len(msd_mean)))
	#print("steps", steps)
	print("")
	fig, axe = plt.subplots(2,3)
	axe[0,0].set_title("g(r)")
	
	for i in range(len(g_r)):
		#print(r)
		axe[0,0].errorbar(r[i],g_r[i],s_g_r[i],fmt='.', ms=1,linewidth=0.5, alpha= 0.5)
		#axe[0,0].scatter(r[i],g_r[i])
	axe[0,0].axvline(1)
	axe[0,0].axhline(1)

	axe[0,1].set_title("msd")
	#for i in range(len(msd)):
		#axe[0,1].plot(np.array(msd_count[i])*100000000,np.array(msd[i])**2)
	axe[0,1].errorbar(times, msd_mean, s_msd_mean, fmt='.')
	#axe[0,1].scatter(times, msd_mean, s=10, marker = "x", alpha = 0.5, c='k')
	axe[0,1].set_ylim(0,max(msd_mean))
	axe[0,1].plot(times[steps>N_start],a*times[steps>N_start]+b)

	axe[1,0].set_title("q6q6")
	for i in range(len(sol_neigh)):
		#axe[1,0].plot(q6q6[i],q6q6_count[i])
		axe[1,0].errorbar(q6q6[i],q6q6_count[i],s_q6q6_count[i],fmt='.')

	axe[1,1].set_title("sol_neigh")
	for i in range(len(sol_neigh)):
		#axe[1,1].plot(sol_neigh[i],sol_neigh_count[i])
		axe[1,1].errorbar(sol_neigh[i],sol_neigh_count[i],s_sol_neigh_count[i], fmt = '.')

	#axe[0,2].set_title("velocity")
	#for i in range(len(vel)):
		#axe[1,1].plot(sol_neigh[i],sol_neigh_count[i])
	#	axe[0,2].errorbar(vel[i],vel_count[i],s_vel_count[i], fmt = '.')

	axe[0,2].set_title('velocity')
	for i in range(len(vel)):
		axe[0,2].scatter(vel_count[i]*np.max(times)*10, vel[i]**2)	
	axe[0,2].scatter(times, vel_2)

	axe[1,2].set_title("solids")
	axe[1,2].scatter(steps,solids)

	#axe[1,2].set_ylim(0,N)

	fig.tight_layout()
	fig.savefig("evaluations/Overview_for_Out_{0}_{1}.png".format(N_part,int(phi_sim*1000)), dpi=800)
	plt.close('all')
	#evaluate and plot diffusion curve
	#print(g_r)
	#print(D_E(phi_sim,np.max(g_r)))
	#print(phi_sim,np.max(g_r))
	return(a,ea, D_E(phi_sim,np.max(g_r[1:])))



if __name__ =="__main__":
	########################### General Numbers and Switches ##############################
	#For use of ths script, there are two switches to either resimulate a series, or do a single_evaluation of all simulations.
	# Furthermore the Numbers can be choosen freely when resimulating, but have to be choosen as during simulation when only evaluating. 
	simulate =True # switch for re simulating, only executes the evaluation if set to False
	single_eval=True

	#Frame numbers during simulation (snapshots). 
	N_eq = 2 #10
	N_pr = 2 #20

	#number of fcc unit cells, and correspondig particle number
	N_cell = 5
	N_part = 4 * N_cell**3

	#densities sampled in the series
	phi=np.linspace(0.54,0.555,16)
	print(phi)
	#sweep size, at the moment defines measurement and snapshot. Later on only should define snapshot
	#should be choose such to fit the right number of steps
	sweep_eq =  1000 #500
	sweep_co =  100 
	sweep_pr =  1000 #1000

	#steps numbers to be executed, see above for choosing them
	steps_eq = sweep_eq * N_eq
	steps_pr = sweep_pr * N_pr

	#collection of step numbers, for later use
	step_array=[steps_eq,steps_pr,sweep_eq,sweep_co,sweep_pr]

	############################ Running Simulations ##############################
	#
	N_process=4 # number of processes run in parallel during simulation
	i = 0

	in_files=[]
	out_files=[]
	simulations=[]
	active=np.zeros(len(phi))

	while i<len(phi):
		if (sum(active)<N_process):
			if simulate: #implementation of simulation switch
				t1,t2,t3=run_sim(phi[i],phi[i],step_array)
				active[i]=1
				in_files.append(t1)
				out_files.append(t2)
				simulations.append(t3)
				print("activated simulation ", i)
			i+=1
	
		else:
			time.sleep(5) #every few seconds the program looks if a process has finished
			for k in range(i):
				if (active[k]==1) & (simulations[k].poll()==0) :
					print("give out simulation1", k)
					write_out_sim(in_files[k],out_files[k],simulations[k])
					active[k]=0

					

	while(sum(active)!=0):
		time.sleep(5) #every few seconds the program looks if a process has finished
		for k in range(len(phi)):
			if (active[k]==1) & (simulations[k].poll()==0) :
				print("give out simulation2", k)
				write_out_sim(in_files[k],out_files[k],simulations[k])
				active[k]=0


	########################## Single Evaluations ##############################
	
	if single_eval:
		diffusion_slope,s_diffusion_slope,D_Enskog=[],[],[]
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

	##################### Overall Diffusion Evaluation #########################
	phi_eval_cutoff = 0.19
	phi_show_cutoff = 0.51
	a2,ea2,b2,eb2,chiq2,corr2 = lineare_regression(np.log(phi[phi<phi_eval_cutoff]),np.log(diffusion_slope[phi<phi_eval_cutoff]),np.abs(s_diffusion_slope[phi<phi_eval_cutoff]/diffusion_slope[phi<phi_eval_cutoff]))	


	plt.figure()
	#plt.errorbar(phi,diffusion_slope,s_diffusion_slope)
	#plt.errorbar(np.log(phi[phi<phi_show_cutoff]),np.log(diffusion_slope[phi<phi_show_cutoff]),np.abs(s_diffusion_slope[phi<phi_show_cutoff]/diffusion_slope[phi<phi_show_cutoff]), fmt = 'x')
	plt.scatter(np.log(phi[phi<phi_show_cutoff]),np.log(diffusion_slope[phi<phi_show_cutoff]), label='measurement',s=20,alpha=0.9)
	plt.plot(np.log(phi[phi<phi_show_cutoff]),(a2*np.log(phi[phi<phi_show_cutoff])+b2),alpha=0.4, linestyle='dashed')
#	plt.xscale('log')
#	plt.yscale('log')
	plt.xlabel(r'log($\phi$)')
	plt.ylabel('log(MSD / Step)')
	plt.title('At {0} particles with {1:4.0f} Steps/Particle,\n the diffusion exponent is {2:4.2f} +- {3:4.2f}'.format(N_part,steps_pr, a2, ea2))
	plt.axvline(np.log(phi[0]))
	plt.axvline(np.log(phi_eval_cutoff))



#Thermodynamic and dynamical properties of the hard sphere system revisited by molecular dynamics simulation Sławomir Pieprzyk, * a Marcus N. Bannerman, Maciej Chudak c and David M. Heyes d

	lit_phi_1=np.array([0.10,0.15,0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
0.75, 0.80, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01])*np.pi/6
 
	lit_D_1=np.array([1.94786 ,1.24636, 0.89989, 0.69308, 0.55523, 0.45527, 0.37845, 0.31663, 0.265009, 0.221105, 0.183259, 0.150317, 0.121081, 0.095425, 0.072934, 0.053298, 0.049608, 0.046203, 0.042720,
0.039478, 0.036254, 0.033286, 0.030326, 0.027419, 0.024746, 0.022134, 0.019695, 0.017357, 0.015082, 0.013013, 0.011053, 0.009222])



#Self-Diffusion Coefficient of the Hard-Sphere Fluid: System Size Dependence andEmpirical Correlations D. M. Heyes* and M. J. Cass
	lit_phi_2=np.array([0.05,0.1,0.15,0.2,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5])
	lit_D_2=np.array([2.039,0.936,0.5785,0.402,0.287,0.24,0.2021,0.168,0.1396,0.1136,0.0877,0.0790,0.0699,0.0620,0.0549,0.0485,0.0413,0.0351,0.0301,0.0249,0.0201])

	print(len(lit_phi_2))
	print(len(lit_D_2))
	plt.scatter(np.log(lit_phi_1),np.log(lit_D_1),label='lit_val_1',s=5,alpha=0.4)
	plt.scatter(np.log(lit_phi_2),np.log(lit_D_2),label='lit_val_2',s=5,alpha=0.4)
	print(D_Enskog)
	plt.scatter(np.log(phi[phi<phi_show_cutoff]),np.log(D_Enskog[phi<phi_show_cutoff]),label='Enskog',alpha=0.3)
	plt.legend()
	plt.savefig('Diffusion_depending_on_density_{0}.pdf'.format(N_part))

	lit_fit_par=np.polyfit(np.log(lit_phi_1),np.log(lit_D_1),7)
	lit_fit=np.poly1d(lit_fit_par)

	fig, axe = plt.subplots(2,sharex=True)
	axe[0].scatter(np.log(lit_phi_1),np.log(lit_D_1))
	axe[0].plot(np.log(lit_phi_1),lit_fit(np.log(lit_phi_1)))
	axe[1].scatter(np.log(phi[phi<phi_show_cutoff]),np.exp( np.log(diffusion_slope[phi<phi_show_cutoff])- lit_fit(np.log(phi[phi<phi_show_cutoff]))))
	axe[1].set_xlim(-3,-0.5)
	#axe[1].set_ylim(3,4)
	
	fig.savefig('comparison.png', dpi=1000)


