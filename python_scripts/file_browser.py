import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import h5py
import os




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
				split_string = filenames[i].split("_")
				N = int(split_string[1])
				rho_i = int(split_string[2])
				rho_f = int(split_string[3])
				eq = int(split_string[4][2:-1])
				pr = int(split_string[5][2:-1])
				seed = int(split_string[6][:-3])
				log_file.write("{0} {1} {2} {3} {4} {5}\n".format(N, rho_i, rho_f, eq, pr, seed))
	
	print("Fine are {0}/{1} = {2:4.2f}".format(fine,len(filenames),fine/len(filenames)))
	if write_out:
		log_file.close()


if __name__ == "__main__":
	show_broken_files("../outfiles/NEMO/16k_november/")

