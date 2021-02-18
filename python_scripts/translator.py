#author : wilkin woehler
#date : 11.05.2020

import h5py
import time as time
import os as os
import numpy as np


#read in h5(md) file from which data is supposed to be read
filename = "../../Out_file1.h5"
#filename = "../../Out_0.75_1.00_sm.h5"
read_file=h5py.File(filename, "r")

#show all folders and subfolders, from which a specific one can be picked
def print_names(name):
    print(name)
read_file.visit(print_names)

#show all folders and subfolders, from which a specific one can be picked
def print_names(name,object_name):
    print(name)
    print(object_name)
read_file.visititems(print_names)


#open a text file to which the data can be written
dump_file1=open("dump_position.txt","w+")
dump_file2=open("dump_velocity.txt","w+")

#take out position data and write it in the dump_file
data1=read_file['particles/test_name/position/value'][()]
data2=read_file['particles/test_name/velocity/value'][()]

for i in range(len(data1[0])):
	dump_file1.write("{0} {1} {2}\n".format(data1[0][i][0],data1[0][i][1],data1[0][i][2]))
	dump_file2.write("{0} {1} {2}\n".format(data2[0][i][0],data2[0][i][1],data2[0][i][2]))
	print("Particle {0}".format(i))
	print("p: {0} {1} {2}".format(data1[0][i][0],data1[0][i][1],data1[0][i][2]))
	print("v: {0} {1} {2}".format(data2[0][i][0],data2[0][i][1],data2[0][i][2]))
dump_file1.close()
dump_file2.close()




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
	#file.close() would be better but 


#opening a file making it ready to be written on, in a not very elegant fashion... and the dumping the position time_series

file_path='test.dump'
try:
	os.remove(file_path)
except:
	print('no file yet there but is created')

new_file=open(file_path,'w+')
new_file.close()


l0=10
sigma=1

lx=l0*sigma
ly=l0*sigma
lz=l0*sigma

#L=[[-lx/2,lx/2],[-ly/2,ly/2],[-lz/2,lz/2]] #box length maybe read it also from h5md file
L=[[0,lx],[0,ly],[0,lz]]

dump_time_1=time.time()
file=open(file_path,'r+')
#print(len(data1),len(data1[0]),len(data1[0][0]))
for i in range(len(data1)):
	dump_file(data1[i,:,:],i,L,file,typ=np.arange(len(data1[0]))+1)
file.close()
dump_time_2=time.time()

print('Time to write down all dump files = {0:4.2f}'.format(dump_time_2-dump_time_1))



