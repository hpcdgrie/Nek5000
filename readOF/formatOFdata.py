# format OpenFOAM data
# easier to be read in Nek5000.
#
# notes:
# special for mus case, because mesh in Nek is rotated to align with axis.
# while OpenFOAM mesh is not. So OpenFOAM coordinates and velocity should
# also be rotated here. 
# But this rotation is not needed for general usage
import os

cxy = [[0 for i in range(2)] for j in range(4)]
cxy[0][0] = 36.88*0.0254
cxy[0][1] = 29.99*0.0254
cxy[1][0] = 44.08*0.0254
cxy[1][1] = 17.78*0.0254
cxy[2][0] = 31.87*0.0254
cxy[2][1] = 10.59*0.0254
cxy[3][0]= 24.67*0.0254
cxy[3][1] = 22.80*0.0254

v32 = [0,0]
v34 = [0,0]

v32[0] = (cxy[1][0] - cxy[2][0])
v32[1] = (cxy[1][1] - cxy[2][1])

v34[0] = (cxy[3][0] - cxy[2][0])
v34[1] = (cxy[3][1] - cxy[2][1])

cosv32 = v32[1]/(v32[0]**2 + v32[1]**2)**0.5
sinv32 = v32[0]/(v32[0]**2 + v32[1]**2)**0.5


casePath = "../min_unit_ke_2percent/"
nProcessors = 64
nSubFiles = 4 # subdivide files in each processor into 16 files
nCellsPerFile = 20000 # max points per file, to reduce memory requirement for Nek
timeFolder = 10000
timeFolder_withC = 10001
maxCells = 310000
C = [[0,0,0] for i in range (maxCells)]
U = [[0,0,0] for i in range (maxCells)]
T = [[0] for i in range (maxCells)]
os.system("mkdir ../OFdata")

for ip in range(0,nProcessors):
	print "converting processor " + str(ip)
	
	processorPath =  casePath+"processor"+str(ip)+"/"
	timePath = processorPath +str(timeFolder)+"/"
	timePath_withC = processorPath +str(timeFolder_withC)+"/"
	writePath = "../OFdata/proc"+str(ip+1)+"/"

	os.system("mkdir "+writePath)

	# first, read C file
	CPath = timePath_withC + "C" 
	Cfile = open(CPath, "r")

	while (True):
		line = Cfile.readline()
		if(len(line)>0):
			if(line.split(" ")[0]=="internalField"):
				break
	line = Cfile.readline()
	numCells = int(line)
	print "cell number: "+str(numCells)
	
	line = Cfile.readline()
	for iCell in range (0,numCells):
		line = Cfile.readline()
		subline = line[1:-2]
		for j in range (0,3):
			C[iCell][j] = float(subline.split()[j])
		# rotating coordinates. this is specific for this case
		xx =  C[iCell][0] - cxy[2][0]
		yy =  C[iCell][1] - cxy[2][1]
		pva =( xx*v32[0] + yy*v32[1]) /((v32[0]**2 + v32[1]**2)**0.5)
		pvb =( xx*v34[0] + yy*v34[1]) /((v34[0]**2 + v34[1]**2)**0.5)
		C[iCell][0] = pva
		C[iCell][1] = pvb 

	Cfile.close()
	iTotCell = 0
	for iSubFile in range(0,nSubFiles):	
		CforNekPath = writePath + "CforNek"+str(iSubFile+1)
		CforNekfile = open(CforNekPath,"w")
		nCells = nCellsPerFile
		if(iSubFile == (nSubFiles-1)):
			nCells = (numCells%nCellsPerFile)
		CforNekfile.write(str(nCells)+"\n")
		for iCell in range (0,nCells):
			CforNekfile.write(str(C[iTotCell][0])+","+str(C[iTotCell][1])+","+str(C[iTotCell][2])+"\n")
			iTotCell= iTotCell+1
		CforNekfile.close()

	# second, read U file
	UPath = timePath + "U" 
	Ufile = open(UPath, "r")

	while (True):
		line = Ufile.readline()
		if(len(line)>0):
			if(line.split(" ")[0]=="internalField"):
				break
	line = Ufile.readline()
	numCells = int(line)
	#print "cell number: "+str(numCells)
	
	line = Ufile.readline()
	for iCell in range (0,numCells):
		line = Ufile.readline()
		subline = line[1:-2]
		for j in range (0,3):
			U[iCell][j] = float(subline.split()[j])
		# rotating coordinates. this is specific for this case
		ux_of = U[iCell][0]
		uy_of = U[iCell][1]
		
		U[iCell][0] = ux_of*cosv32+uy_of*sinv32
		U[iCell][1] = uy_of*cosv32-ux_of*sinv32

	Ufile.close()
	iTotCell = 0
	for iSubFile in range(0,nSubFiles):	
		UforNekPath = writePath + "UforNek"+str(iSubFile+1)
		UforNekfile = open(UforNekPath,"w")
		nCells = nCellsPerFile
		if(iSubFile == (nSubFiles-1)):
			nCells = (numCells%nCellsPerFile)
		UforNekfile.write(str(nCells)+"\n")
		for iCell in range (0,nCells):
			UforNekfile.write(str(U[iTotCell][0])+","+str(U[iTotCell][1])+","+str(U[iTotCell][2])+"\n")
			iTotCell= iTotCell+1
		UforNekfile.close()

	# third, read T file
	TPath = timePath + "T" 
	Tfile = open(TPath, "r")

	while (True):
		line = Tfile.readline()
		if(len(line)>0):
			if(line.split(" ")[0]=="internalField"):
				break
	line = Tfile.readline()
	numCells = int(line)
	#print "cell number: "+str(numCells)
	
	line = Tfile.readline()
	for iCell in range (0,numCells):
		line = Tfile.readline()
		T[iCell] = float(line)

	Tfile.close()
	iTotCell = 0
	for iSubFile in range(0,nSubFiles):	
		TforNekPath = writePath + "TforNek"+str(iSubFile+1)
		TforNekfile = open(TforNekPath,"w")
		nCells = nCellsPerFile
		if(iSubFile == (nSubFiles-1)):
			nCells = (numCells%nCellsPerFile)
		TforNekfile.write(str(nCells)+"\n")
		for iCell in range (0,nCells):
			TforNekfile.write(str(T[iTotCell])+"\n")
			iTotCell= iTotCell+1
		TforNekfile.close()
		
		
#	# first, read U file
#	UPath = timePath + "U" 
#	UforNekPath = timePath + "UforNek" 
#	Ufile = open(UPath, "r")
#	UforNekfile = open(UforNekPath,"w")
#
#	while (True):
#		line = Ufile.readline()
#		#print line
#		if(len(line)>0):
#			if(line.split(" ")[0]=="internalField"):
#				break
#	line = Ufile.readline()
#	numCells = int(line)
#	#print "cell number: "+str(numCells)
#	UforNekfile.write(str(numCells)+"\n")
#	
#	line = Ufile.readline()
#	#print line
#	for iCell in range (0,numCells):
#		line = Ufile.readline()
#		subline = line[1:-2]
#		for j in range (0,3):
#			U[j] = float(subline.split()[j])
#		ux_of = U[0]
#		uy_of = U[1]
#		
#		U[0] = ux_of*cosv32+uy_of*sinv32
#		U[1] = uy_of*cosv32-ux_of*sinv32
#		
#		UforNekfile.write(str(U[0])+","+str(U[1])+","+str(U[2])+"\n")
#	
#	Ufile.close()
#	UforNekfile.close()
#	
#	# second, read T file
#	TPath = timePath + "T" 
#	TforNekPath = timePath + "TforNek" 
#	Tfile = open(TPath, "r")
#	TforNekfile = open(TforNekPath,"w")
#
#	while (True):
#		line = Tfile.readline()
#		#print line
#		if(len(line)>0):
#			if(line.split(" ")[0]=="internalField"):
#				break
#	line = Tfile.readline()
#	numCells = int(line)
#	#print "cell number: "+str(numCells)
#	TforNekfile.write(str(numCells)+"\n")
#	
#	line = Tfile.readline()
#	#print line
#	for iline in range (0,numCells):
#		line = Tfile.readline()
#		TforNekfile.write(line)
#	
#	Tfile.close()
#	TforNekfile.close()
	
