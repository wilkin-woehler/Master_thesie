from h5md_python import *

def evaluate_run():
	base_path="../outfiles/NEMO/1m_pre_results/"

	p1="rho_534/"
	p2="rho_533/"
	p3="rho_532/"
	p4="rho_531/"

	paths=[p1,p2,p3,p4]

	for i in range(len(paths)):
		paths[i]=base_path+paths[i]

	pos1=[0,1,0,1]
	pos2=[0,0,1,1]

	names=["rho_534","rho_533","rho_532","rho_531"]

	run = h5md_DataRun(paths, pos1, pos2, names, [2000,200,100], multi_threading = 4, run_name = "pre_results_14_12_2020", evaluation_number = -1)

	run.collect_transition_times()
	run.plot_transition_times()

	run.collect_cluster_distribution()
	run.plot_cluster_distribution()

def test_eval():
	base_path="../outfiles/NEMO/1m_pre_results/rho_534/"	
	missed=0
	for i in range(20):
		#try:
		if True:
			test_set=h5md_dataset(base_path,4*64**3,0.45,0.534,1000,20000,i+1,[4000,600,100])
			test_set.eval_ClusterDistribution()
			test_set.eval_SolidFraction()
			test_set.eval_LargestCluster()
			test_set.eval_MSD()
			test_set.plot_measurement("")
		#except:
		#	print('Missed: ', i)
		#	missed+=1
	print("Overall missed: ", missed, " sets")
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
	test_eval()
#	evaluate_run()




