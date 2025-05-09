import sys, string

class Spectrum:
    '''
    A class to store the variable amplitude loading.  Reads "infile" in the
    constructor.  First line of "infile"
    '''

    def __init__(self,infile):
        self.load=[] # pairs of (min,max) percent of maximum applied stress

        inf = open(infile,'r')
        vals = inf.readlines()
        inf.close()
        for val in vals:
            fff=string.splitfields(val)
            fff.sort()
            self.load+=[(float(fff[0]),float(fff[1]))]

        self.length=len(self.load)

        self.firsttime = 1 # boolean to indicate this is the first time through
                           # self.load.  Use in Delta to remove zig-zags.

    def Initialize(self):
        self.count=0

    def Delta(self):
        '''
        return the next delta sigma in the list.  If at the end of the list,
        start over.
        '''

        # strip these out of self.load for easier access below...
        Kmaxc = self.load[self.count][1] # c for current
        Kminc = self.load[self.count][0]

        # first time through self.load, remove any zig-zags where:
        #    Kmax,i ~ Kmin,i+1 or
        #    Kmin,i ~ Kmax,i+1
        zero = 0.0001
        self.count+=1
        if self.count < self.length - 1:
            if self.firsttime:
                Kmaxn = self.load[self.count][1] # n for next... 
                Kminn = self.load[self.count][0]
                if abs(Kmaxc - Kminn) <= zero:
                    self.RemoveZig()
                    return Kmaxn - Kminc, Kminc/Kmaxn
                elif abs(Kmaxn - Kminc) <= zero:
                    self.RemoveZag()
                    return Kmaxc - Kminn, Kminn/Kmaxc

        else:
            self.firsttime = 0 
            self.count = 0

        return Kmaxc - Kminc, Kminc/Kmaxc

        # old... 11/15/05... 
##        # now figure out Delta and return it... 
##        if self.count < self.length:
##            if abs(self.load[self.count][0]) > abs(self.load[self.count][1]):
##                if Kmaxc > 0.0:
##                    DK = self.load[self.count][1]
##                else: DK = self.load[self.count][0] - self.load[self.count][1]
##            else: DK = self.load[self.count][1]-self.load[self.count][0]
##            self.count+=1
##        else:
##            if abs(self.load[self.count][0]) > abs(self.load[self.count][1]):
##                if self.load[self.count][1] > 0.0:
##                    DK = self.load[self.count][1]
##                else: DK = self.load[self.count][0]-self.load[self.count][1]
##            else: DK = self.load[self.count][1]-self.load[self.count][0]
##            self.count=0
##
##        return DK

    def RemoveZag(self):
        '''
        method to condense self.load when:
            Kmin,i ~ Kmax,i+1
        '''

        (Kminc,Kmaxc) = self.load.pop(self.count-1)
        (Kminn,Kmaxn) = self.load.pop(self.count-1)
        self.load.insert(self.count-1,(Kminn,Kmaxc))
        self.length = len(self.load)

    def RemoveZig(self):
        '''
        method to condense self.load when:
            Kmax,i ~ Kmin,i+1
        '''

        (Kminc,Kmaxc) = self.load.pop(self.count-1)
        (Kminn,Kmaxn) = self.load.pop(self.count-1)
        self.load.insert(self.count-1,(Kminc,Kmaxn))
        self.length = len(self.load)

if __name__=="__main__":
    spec = Spectrum('sips3002.val')
    dk=[]; R=[]
    spec.Initialize()
    for i in range(5300):
        dk1,R1=spec.Delta()
        dk+=[dk1]
        R+=[R1]

    print spec.load[0:15]
    print dk[0:15]
    print R[0:15]
    print ' '
    print spec.load[4950:]
    print dk[4950:4957]
    print R[4950:4957]









