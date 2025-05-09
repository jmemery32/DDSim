import sys, string

class Parameters:
    '''
    a class to store the DDSim model parameters from *.par input file.  
    '''

    def __init__(self,filename,path):
        self.monte = 0 # 0 = deterministic; 1 = from distribution;
                       # 2 = particle cracking filter
        self.filename = filename
        self.path = path
        self.material=[]
        self.a_b =[[0.,0.]]
        self.dN=1.0
        self.max_error=0
        self.N_max = 0
        self.min_inc = 0.0
        self.samples = 5
        self.sets = 1
        self.r = 1.1 # default fraction to increase a for building K vs. a
                     # also used in VarAmp as monitor to recompute K
        self.Max_crack_size = 10.0 # default max crack size
        self.__GetParameters()

    def __repr__(self):
        s="  -->  The input parameters are: \n"
        s+=" path: "+self.path+"\n"
        s+=" file: "+self.filename+"\n"
        s+=" Monte = %1i, a = %10.6e, b = %10.6e \n" % \
            (self.monte,self.a_b[0][0],self.a_b[0][1])
        s+=" dN = %6.3e, max_error = %6.3e \n N_max = %6.3e" % \
            (self.dN,self.max_error,self.N_max)
        s+=", min_inc = %3.6e \n" % (self.min_inc)
        s+=" Sets = %i, Samples = %i, r = %3.3e, Max_crack_size = %8.6e" % \
            (self.sets,self.samples,self.r,self.Max_crack_size)
        s+="\n"
        s+=' NASGRO material properties are: \n'
        p=self.__MaterialsRepr()
        s+=p
        return s

    def __MaterialsRepr(self):
        s=' UTS = '+str(self.material[0])
        s+=' YS = '+str(self.material[1])
        s+=' Kie = '+str(self.material[2])
        s+=' Kic = '+str(self.material[3])+"\n"
        s+=' Ak = '+str(self.material[4])
        s+=' Bk = '+str(self.material[5])
        s+=' C = '+str(self.material[6])
        s+=' n = '+str(self.material[7])+"\n"
        s+=" n' = "+str(self.material[18])
        s+=' p = '+str(self.material[8])
        s+=' q = '+str(self.material[9])
        s+=' DK1 = '+str(self.material[10])+"\n"
        s+=' Cthp = '+str(self.material[11])
        s+=' Cthm = '+str(self.material[12])
        s+=' Rcl = '+str(self.material[13])
        s+=' alpha = '+str(self.material[14])+"\n"
        s+=' Smax/SIGo = '+str(self.material[15])
        s+=' R = '+str(self.material[16])
        s+=' a0 = '+str(self.material[17])
        s+=' damp = '+str(self.material[19])+"\n"
        try: s+=' Rso = '+str(self.material[20])+"\n"
        except IndexError: s+=' Rso = 3.0 \n'
        return s

    def __GetParameters(self):
        file = open(self.path+self.filename+'.par','r')
        buff = self.__GetLine(file)
        while buff != None:
            vals = string.splitfields(buff)
            if vals[0] == 'a_b':
                self.a_b =[]
                vals = string.splitfields(self.__GetLine(file))
                for ii in range(len(vals)/2): self.a_b +=[[float(vals[ii*2]), \
                                                        float(vals[ii*2+1])]]
            elif vals[0] == 'material':
                vals = string.splitfields(self.__GetLine(file))
                for val in vals:
                    self.material += [float(val)]
            elif vals[0] == 'dN':
                vals = string.splitfields(self.__GetLine(file))
                self.dN = float(vals[0])
            elif vals[0] == 'max_error':
                vals = string.splitfields(self.__GetLine(file))
                self.max_error = float(vals[0])
            elif vals[0] == 'N_max':
                vals = string.splitfields(self.__GetLine(file))
                self.N_max = float(vals[0])
            elif vals[0] == 'min_inc':
                vals = string.splitfields(self.__GetLine(file))
                self.min_inc = float(vals[0])
            elif vals[0] == 'samples':
                vals = string.splitfields(self.__GetLine(file))
                self.samples = int(vals[0])
            elif vals[0] == 'sets':
                vals = string.splitfields(self.__GetLine(file))
                self.sets = int(vals[0])
            elif vals[0] == 'monte':
                vals = string.splitfields(self.__GetLine(file))
                self.monte = int(vals[0])
            elif vals[0] == 'ratio':
                vals = string.splitfields(self.__GetLine(file))
                self.r = float(vals[0])
                if self.r < 1.0: self.r = 1.1
            elif vals[0] == 'MaxCrack':
                vals = string.splitfields(self.__GetLine(file))
                self.Max_crack_size = float(vals[0])

            buff = self.__GetLine(file)

    def __GetLine(self,file):
        """
        this is a little helper function that scans the input
        file and ignores lines starting with a '#'
        """
        buff = file.readline()
        if len(buff) == 0: return None
        while buff[0] == '#':
            buff = file.readline()
            if len(buff) == 0: return None
        return buff    
