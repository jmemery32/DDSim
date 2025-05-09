import sys
import math

class container:
	foo = 1
	

def readMscale(In):
	nextline = 0
	Out = container()
	for line in In:
		if line.find("Center:") > -1:
			nextline = 1
		elif nextline == 1:
			nextline = 0
			Out.center = [float(x) for x in line.split(" ") if x!='']
	
	numnodes = 0
	Out.nodes = dict()
	In.seek(0)
	for line in In:
		lin = line.split(" ")
		if lin[0].find("Nodes:") >-1:
			numnodes = int(lin[1])
		elif numnodes > 0:
			numnodes = numnodes-1
			key = int(lin.pop(0))
			Out.nodes[key] = [float(x) for x in lin if x!='']
	
	numdata = 0
	numnodes = 0
	In.seek(0)
	for line in In:
		lin = line.split(" ")
		if lin[0].find("Data:") >-1:
			numdata = int(lin[1])
			Out.Datums = dict()
		elif (numdata>0 and numnodes==0):
			numdata = numdata-1
			currentData = lin[0]
			numnodes = int(lin[2])
			Out.Datums[currentData] = dict()
		elif numnodes > 0:
			numnodes = numnodes-1
			key = int(lin.pop(0))
			Out.Datums[currentData][key] = [float(x) for x in lin if x!='']
	
	return Out
			
def readNodes(In,fac,center):
	nodes = dict()
	for line in In:
		lin = line.split(" ")
		key = int(lin.pop(0))
		nodes[key] = [float(x)*fac for x in lin if x!='']
		for i in xrange(len(nodes[key])):
			nodes[key][i] = nodes[key][i] + center[i]
	return nodes

def dist(a,b):
	return math.sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]))

def Interp(input,out):
	struct = container()
	struct.Datums = dict()
	struct.nodes = dict()
	for z in input.Datums.keys():
		struct.Datums[z] = dict()
		for x in out.keys():
			struct.Datums[z][x] = [0 for a in input.Datums[z][input.Datums[z].keys()[0]]]
	for x in out.keys():  #out nodes
		sigma = 0
		for y in input.nodes.keys():  #in nodes
			D = dist(out[x],input.nodes[y])
			weight = 1/D/D
			sigma = sigma + weight
			for z in struct.Datums.keys():  #data fields
				upd = [weight*a for a in input.Datums[z][y]]
				struct.Datums[z][x] = [struct.Datums[z][x][i]+upd[i] for i in xrange(len(upd))]
		for z in struct.Datums.keys():
			struct.Datums[z][x] =  [a/sigma for a in struct.Datums[z][x]]
		struct.nodes[x] = [a for a in out[x]]
	return struct
			
def writeFiles(out,name):
	#write the new nodes file
	nodfile = file(name+".nod",'w')
	for x in out.nodes.keys():
		string = str(x) + " " + str(out.nodes[x][0]) +" "+str(out.nodes[x][1])+" "+str(out.nodes[x][2])+"\n"
		nodfile.write(string)
	nodfile.close()
	# do the datums
	for y in out.Datums.keys():
		ext = ""
		dataname = ""
		if (y.find("STRESS") >-1):
			ext = ".sig"
			dataname = "#STRESS 2\n"
		elif (y.find("STRAIN") >-1):
			ext = ".stn"
			dataname = "#STRAIN 2\n"
		elif (y.find("FORCE") >-1):
			ext = ".frc"
			dataname = "#FORCE 1\n"
		elif (y.find("DISPLACEMENT") >-1):
			ext = ".dsp"
			dataname = "#DISPLACEMENT 1\n"
		elif (y.find("TEMPERATURE") >-1):
			ext = ".tmp"
			dataname = "#TEMPERATURE 0\n"
		elif (y.find("LIFE") >-1):
			ext = ".lif"
			dataname = "#LIFE 0\n"
		datafile = file(name+ext,'w')
		datafile.write(dataname)
		for x in out.Datums[y].keys():
			string = str(x)
			for z in out.Datums[y][x]:
				string = string + " "
				string = string + str(z)
			string = string +"\n"
			datafile.write(string)
		datafile.close()

	

if __name__ == "__main__":
	'''
	drive the system...
	'''

	if '-mscale' in sys.argv:
		index = sys.argv.index('-mscale') + 1
		mscaleFile = sys.argv[index]

	if '-nod' in sys.argv:
		index = sys.argv.index('-nod') + 1
		nodFile = sys.argv[index]

	if '-targetScale' in sys.argv:
		index = sys.argv.index('-targetScale') + 1
		targScale = float(sys.argv[index])
	else:
		targScale = 1;

	if '-outbase' in sys.argv:
		index = sys.argv.index('-outbase') + 1
		outbase = sys.argv[index]
	else:
		outbase = 'output'

	mscaleIn = file(mscaleFile,'r')
	nodIn = file(nodFile,'r')

	indata = readMscale(mscaleIn)
	outNodes = readNodes(nodIn,targScale,indata.center)

	outStruct = Interp(indata,outNodes)
	writeFiles(outStruct,outbase)
	
