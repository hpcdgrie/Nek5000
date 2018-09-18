# this file will split a bit stl file into serveral small ones for 
# as a stl file is usually too large for a processor memory.
# 
#
from math import *
import os

ntr = 672 # total number of triangles
npf = 672   # numper of triagnles per file
nfiles = int(ceil(ntr/npf)+1.0)

stl_all = open(os.getcwd() +'/stl_all.dat','r')
stl_all.readline() # skip first line

for ifile in range(0,nfiles):
	stl_split = open(os.getcwd() +'/stl_'+str(ifile+1)+'.dat','w')

	for itr in range(0,npf):

		line=stl_all.readline()
		stl_split.write(str(line[15:].split()[0])+','+str(line[15:].split()[1])+','+str(line[15:].split()[2]+'\n'))
		line=stl_all.readline()
		line=stl_all.readline()
		stl_split.write(str(line[13:].split()[0])+','+str(line[13:].split()[1])+','+str(line[13:].split()[2]+'\n'))
		line=stl_all.readline()
		stl_split.write(str(line[13:].split()[0])+','+str(line[13:].split()[1])+','+str(line[13:].split()[2]+'\n'))
		line=stl_all.readline()
		stl_split.write(str(line[13:].split()[0])+','+str(line[13:].split()[1])+','+str(line[13:].split()[2]+'\n'))
		line=stl_all.readline()
		line=stl_all.readline()
	stl_split.close()
