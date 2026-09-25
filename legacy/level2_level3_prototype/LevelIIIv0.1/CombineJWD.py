import sys, string



#####################

def HelpPrints():
    print " " 
    print " This script reads two finite element models, the Global and the "
    print " Local, in RDB format and the Local.mpc file.  Then combines them "
    print " by determining the highest node and element number in the local "
    print " model and incrementing all numbering in the global model by as "
    print " much.  "
    print " ---- " 
    print " Usage:"
    print " C:\...>Combine.py [global filename] [local filename] [output filename ](no extension)] "
    print " "
    print " Output:"
    print " Output.con, .nod, etc."
    print ' '

###########################################

class MPCFileJoiner:
	'''
	this class drives the combining of models
	'''

	def __init__(self,Gfile,Lfile,Ofile):
		'''
		load up the info in nodes, and write it out with appropriate increments

		do the same for elements (including incrementing the nodes)

		yada yada yada
		'''

		self.G = Gfile
		self.L = Lfile
		self.O = Ofile

		self.GlobalNodes = self.ReadData(Gfile+'.nod')		
		self.LocalNodes = self.ReadData(Lfile+'.nod')
		self.NodeInc = max(self.LocalNodes.keys())
		self.Nodes = dict([(x+self.NodeInc,self.GlobalNodes[x]) for x in self.GlobalNodes.keys()])
		self.Nodes.update(self.LocalNodes)# = Merge(self.LocalNodes,self.GNodes)
		print 'error checking: ', len(self.LocalNodes) + len(self.GlobalNodes),' = ',len(self.Nodes)
		self.WriteDict(self.Nodes,Ofile+'.nod')

		self.ShapeInc = int(0)
		try:
			self.GlobalShape = self.ReadData(Gfile+'.ShapeMap')
			self.LocalShape = self.ReadData(Lfile+'.ShapeMap')
			self.ShapeInc = max(self.LocalShape.keys())
			self.Shape = dict([(x+self.ShapeInc,self.GlobalShape[x]) for x in self.GlobalShape.keys()])
			self.Shape.update(self.LocalShape)
			self.WriteDict(self.Shape,Ofile+'.ShapeMap')
		except IOError: pass

		self.MatInc = int(0)
		try:
			self.GlobalMat = self.ReadData(Gfile+'.Materials')
			self.LocalMat  = self.ReadData(Lfile+'.Materials')
			self.MatInc = max(self.LocalMat.keys())
			self.Mat = dict([(x+self.MatInc,self.GlobalMat[x]) for x in self.GlobalMat.keys()])
			self.Mat.update(self.LocalMat)
			self.WriteDict(self.Mat,Ofile+'.Materials')
		except IOError: pass
		
		self.GlobalElems = self.ReadData(Gfile+'.con')
		self.LocalElems = self.ReadData(Lfile+'.con')
		self.ElemInc = max(self.LocalElems.keys())
		self.Elems = dict([(x+self.ElemInc,[int(self.GlobalElems[x][0])+self.ShapeInc]+[int(self.GlobalElems[x][1])+self.MatInc]+[self.GlobalElems[x][2]]+[y+self.NodeInc for y in self.GlobalElems[x][3:]]) for x in self.GlobalElems.keys()])
		#self.Elems = dict([(x+self.ElemInc,[self.GlobalElems[x][0]+self.ShapeInc]+[self.GlobalElems[x][1]+self.MatInc]+[self.GlobalElems[x][2]]+[y+self.NodeInc for y in self.GlobalElems[x][3:]]) for x in self.GlobalElems.keys()])
		self.Elems.update(self.LocalElems)
		self.WriteDict(self.Elems,Ofile+'.con')

		try:
			self.GlobalComp = self.ReadData(Gfile+'.edg')
			self.LocalComp = self.ReadData(Lfile+'.edg')
			self.Comps = dict([(x+self.NodeInc,self.GlobalComp[x][:2]+[z+self.NodeInc for z in self.GlobalComp[x][2:]]) for x in self.GlobalComp.keys()])
			self.Comps.update(self.LocalComp)
			self.WriteDict(self.Comps,Ofile+'.edg')
		except IOError: pass

		try:
			self.GlobalSig = self.ReadData(Gfile+'.sig')
			self.LocalSig  = self.ReadData(Lfile+'.sig')
			self.Sigs = dict([(x+self.NodeInc,self.GlobalSig[x]) for x in self.GlobalSig.keys()])
			self.Sigs.update(self.LocalSig)
			self.WriteDict(self.Sigs,Ofile+'.sig')
		except IOError: pass
		


################################################################
	
	def ReadData(self,file):
		'''
		load the contents of a keyed file into dic
		'''
		
		print "reading ", file

		dic = dict()
		fin = open(file,'r')
		buff = fin.readline()
		while buff:
			list = string.split(buff)
			buff = fin.readline()

			dic[list[0]] = list[1:]

		fin.close();
		return dic

##############################################################

	def WriteDict(self,dic,file):
		'''
		Write teh contents of dic to file
		'''

		print "Writing ", file
		outfile = open(file,'w')
		for k, e in dic.iteritems():
			foo = ' '
			for x in e:
				foo += str(x) + ' '
			foo+="\n"
			outfile.write(k+foo)
		outfile.close();

	
################################################################



if __name__=="__main__":
	print " ----------------------------------------------------------- "
	print " | this is a test version, don't expect too much            | "
	print " ----------------------------------------------------------- "

	# do you want help?
	if "-help" in sys.argv or "-h" in sys.argv:
		HelpPrints()
		sys.exit()

	try:
		Globalfilename = sys.argv[1]
		Localfilename  = sys.argv[2]
		Outputfilename = sys.argv[3]
	except IndexError:
		HelpPrints()
		sys.exit()
	
	# try to merge it
	merger = MPCFileJoiner(Globalfilename,Localfilename,Outputfilename)
	#done

		
