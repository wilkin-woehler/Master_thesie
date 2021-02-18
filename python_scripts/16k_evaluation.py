from h5md_python import *

def single_eval():
	#evaluates single trajectories, one by one

	outfile_path="../outfiles/"
	#outfile_path="../../"
	numbers=np.linspace(64,64,1,dtype=int)**3*4
	#numbers[4]=64**3*4
	print("N from ", numbers, "are going to be evaluated ")
	rho_inital=np.array([0.45])
	rho=np.linspace(0.532,0.532,1) #rho_final in production phase
	steps_eq=1000
	steps_pr=20000
	
	n_stat=1
	start_seed=472
	seeds=np.arange(start_seed,start_seed+n_stat)
	print("seeds from ", seeds, "are going to be evaluated ")
	
	
	missed=0
	for i in range(len(numbers)):
		for j in range(len(seeds)):
			for k in range(len(rho)):
				for l in range(len(rho_inital)):
					#try:
					if True:
						datenset=h5md_dataset(outfile_path,numbers[i],rho_inital[l],rho[k],steps_eq,steps_pr,seeds[j], [2000,600,100])
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
						datenset.close()
					#except:
					#	print("Could not generate dataset")
					#	missed+=1
				
	print("Missed {0} sets.".format(missed))		

def series_eval(path):
	#evaluates a single bunch of trajectories
	series=h5md_DataSeries(path)
	series.collect_largest_cluster(number=5)
	series.collect_cluster_distribution()
	series.collect_solid_transition()

def run_eval(paths):
	#evaluates bunches of trajectories for different start parameters
	pos2=[0,0,0,1,1,1,1,0]
	pos1=[0,1,3,0,1,2,3,2]
	names=["rho_i_300","rho_i_400","rho_i_480","eq_steps_100","eq_steps_200","eq_steps_500","eq_steps_20k","eq_steps_5k"]
	run = h5md_DataRun(paths, pos1, pos2, names, run_name = "_test")
	run.collect_transition_times()
	run.plot_transition_times()

	#run.collect_cluster_distribution()
	#run.plot_cluster_distribution()

def evaluate_run1():
	p1="../outfiles/NEMO/16k_november/rho_i_300/"
	p2="../outfiles/NEMO/16k_november/rho_i_400/"
	p3="../outfiles/NEMO/16k_november/rho_i_480/"
	p4="../outfiles/NEMO/16k_november/eq_steps_100/1/"
	p5="../outfiles/NEMO/16k_november/eq_steps_200/1/"
	p6="../outfiles/NEMO/16k_november/eq_steps_500/"
	p7="../outfiles/NEMO/16k_november/eq_steps_20k/"
	p8="../outfiles/NEMO/16k_november/base_series/"
	#p9="../outfiles/NEMO/16k_november/eq_steps_100/2/"
	#p10="../outfiles/NEMO/16k_november/eq_steps_200/2/"
	#series_eval(p1)

	paths=[p1,p2,p3,p4,p5,p6,p7,p8]

	pos1=[0,1,3,0,1,2,3,2]
	pos2=[0,0,0,1,1,1,1,0]

	names=["rho_i_300","rho_i_400","rho_i_480","eq_steps_100","eq_steps_200","eq_steps_500","eq_steps_20k","eq_steps_5k"]
	run = h5md_DataRun(paths, pos1, pos2, names, multi_threading = 8, run_name = "_run1")
	run.collect_transition_times()
	run.plot_transition_times()


def evaluate_run2():
	base_path="../outfiles/NEMO/16_k_december/"

	p1="rho_i_300/"
	p2="rho_i_400/"
	p3="rho_i_480/"
	p4="eq_steps_100/"
	p5="eq_steps_200/"
	p6="eq_steps_500/"
	p7="eq_steps_20k/"
	p8="eq_steps_5k/"

	paths=[p1,p2,p3,p4,p5,p6,p7,p8]
	
	for i in range(len(paths)):
		paths[i]=base_path+paths[i]
	
	pos1=[0,1,3,0,1,2,3,2]
	pos2=[0,0,0,1,1,1,1,0]

	names=["rho_i_300","rho_i_400","rho_i_480","eq_steps_100","eq_steps_200","eq_steps_500","eq_steps_20k","eq_steps_5k"]
	run = h5md_DataRun(paths, pos1, pos2, names, multi_threading = 8, run_name = "_run2", evaluation_number = -1)
	run.collect_transition_times()
	run.plot_transition_times()

	#run.collect_cluster_distribution()
	#run.plot_cluster_distribution()

def evaluate_run_4k():
	base_path="../outfiles/NEMO/4k_run_out/"

	p1="rho_i_300/"
	p2="rho_i_400/"
	p3="rho_i_480/"
	p4="eq_steps_100/"
	p5="eq_steps_200/"
	p6="eq_steps_500/"
	p7="eq_steps_20k/"
	p8="eq_steps_5k/"

	paths=[p1,p2,p3,p4,p5,p6,p7,p8]
	
	for i in range(len(paths)):
		paths[i]=base_path+paths[i]
	
	pos1=[0,1,3,0,1,2,3,2]
	pos2=[0,0,0,1,1,1,1,0]

	names=["rho_i_300","rho_i_400","rho_i_480","eq_steps_100","eq_steps_200","eq_steps_500","eq_steps_20k","eq_steps_5k"]
	run = h5md_DataRun(paths, pos1, pos2, names, multi_threading = 8, run_name = "_4k_run", evaluation_number = -1)

	run.collect_transition_times()
	run.plot_transition_times()

	#run.collect_cluster_distribution()
	#run.plot_cluster_distribution()




def test():
	p9="../outfiles/NEMO/16k_november/eq_steps_100/2/"
	p10="../outfiles/NEMO/16k_november/eq_steps_200/2/"
	#series_eval(p1)

	paths=[p9,p10]

	#evaluates bunches of trajectories for different start parameters
	pos2=[1,0]
	pos1=[0,1]
	names=["eq_steps_100","eq_steps_200"]
	run = h5md_DataRun(paths, pos1, pos2, names, run_name = "_test", multi_threading= 2, evaluation_number=5)
	#run = h5md_DataRun(paths, pos1, pos2, names, run_name = "_test")
	run.collect_transition_times()
	run.plot_transition_times()

	run.collect_cluster_distribution()
	run.plot_cluster_distribution()


def evaluate_pre_run_1m():
	path="../outfiles/NEMO/1m_pre_runs/"
	
	rho_531_series = h5md_DataSeries(path+"rho_531/")
	#rho_531_series.collect_largest_cluster()
	rho_531_series.collect_pnt()

	rho_532_series = h5md_DataSeries(path+"rho_532/")
	#rho_532_series.collect_largest_cluster()
	rho_532_series.collect_pnt()

	rho_534_series = h5md_DataSeries(path+"rho_534/")
	#rho_534_series.collect_largest_cluster()
	rho_534_series.collect_pnt()

if __name__ == "__main__":
	#evaluate_run1()

	#evaluate_run_4k()
	#evaluate_run2()
	
	#test()
	single_eval()
	#evaluate_pre_run_1m()




