import numpy as Numeric
import sys, cPickle

def CalcSampMean(set,values):
    '''
    To calculate the sample mean of r.v., self.SampMean.append.
    '''
    mean = Numeric.cumsum(values)[-1]/len(values)
    return mean

def CalcSampVar(set,values,sampmean):
    '''
    To calculate the sample variance of r.v., self.SampVar.append.
    '''

    sum=0.0

    for i in values:
        sum=sum+(i-sampmean)**2.0

    try: var = sum/(float(len(values))-1.0)

    # i.e. only one sample, can happen with particle filter if only one
    # particle breaks... Gerd's does the statistics there with the whole
    # set of particles anyway, so this is not used... 
    except ZeroDivisionError: var= 0.0 

    return var

def CalcProb(set,x,values):
    '''
    To cacluate the probability P(X<=x). I.e., number of samples in rv that
    are less than or equal to x divided by the total number of samples.
    '''
    less=0.0
    for i in values:
        if i <= x:
            less+=1.0
    return less/(float(len(values)))

def TestStata():
    stat = Stata(2,5)
    list1 = [1.2, 3.0, 3.1, 3.4, 5.4]
    list2 = [1253.4, 1485.1, 6984.2, 4512.2, 7845.3]
    print ' list 1 = ', list1
    print ' mean = 3.22  variance = 2.232'
    print ' list 2 = ', list2
    print ' mean = 4416.04  variance = 9239304.43'

    stat.UpdateSet()
    for i in list1: stat.UpdateSamples(i,0)
    stat.UpdateSet()
    for j in list2: stat.UpdateSamples(j,1)

    print ' '
    print ' As stored...'
    print ' list 1 = ', stat.sample(0)
    print ' list 2 = ', stat.sample(1)
    print ' or     = ', stat

    print ' '
    print ' Statistics...'
    print ' Sample Mean and Variance of list 1 = ', \
          stat.SampleMean(0), stat.SampleVariance(0)
    print ' Sample Mean and Variance of list 2 = ', \
          stat.SampleMean(1), stat.SampleVariance(1)

    print ' P(X <= 3.2) =', stat.CalcProb(0,3.2)
    print ' P(X <= 3000) =', stat.CalcProb(1,3000.)

def TestStatN():
    stat = StatN(2,5)
    list1 = [1.2, 3.0, 3.1, 3.4, 5.4]
    list2 = [1253.4, 1485.1, 6984.2, 4512.2, 7845.3]
    print ' list 1 = ', list1
    print ' mean = 3.22  variance = 2.232'
    print ' list 2 = ', list2
    print ' mean = 4416.04  variance = 9239304.43'

    stat.UpdateSet()
    for i in range(len(list1)): stat.UpdateSamples(list1[i],0,i)
    stat.UpdateSet()
    for j in range(len(list2)): stat.UpdateSamples(list2[j],1,i+j+1)

    print ' '
    print ' As stored...'
    print ' list 1 = ', stat.sample(0)
    print ' list 2 = ', stat.sample(1)
    print ' or     = ', stat

    print ' '
    print ' Statistics...'
    print ' Sample Mean and Variance of list 1 = ', \
          stat.SampleMean(0), stat.SampleVariance(0)
    print ' Sample Mean and Variance of list 2 = ', \
          stat.SampleMean(1), stat.SampleVariance(1)

    print ' P(X <= 3.2) =', stat.CalcProb(0,3.2)
    print ' P(X <= 3000) =', stat.CalcProb(1,3000.)

class Stata:
    '''
    This class stores and operates on sets of samples of a random
    variable.  Written in support of DamModel to be stored with
    DamModel.__DamOro[doid].
    '''

    def __init__(self,sets,samples):
        self.rv=[] # to store all samples of the random variable
        # [{a0:[RID0,RIDn],a1:[RID1...aN-1:[RIDN-1]},{aN:RIDN,aN+1:RIDN+1...}]
        # M sets of N.  NOTE:  keys have LISTS of RIDs because it is possible
        # to have a0=an.  
        self.SampMean={} # to store the sample mean of each set, setid:mean
        self.SampVar={} # to store the sample variance of each set, setid:var
        self.sets=sets # integer number of sets
        self.samples=samples # integer number of samples
        self.RID=0
        self.rid_ai={} # mapping rid_ai[rid] -> [setid,ai]

################ class Stata

    def __repr__(self):
        s = ""
        setid=0
        for diction in self.rv:
            s+="Set %s" % (setid)
            s+="\n"
            iters=diction.keys()
            iters.sort()
            for ai in iters:
                for rid in diction[ai]:
                    s+="(%i, %1.4e) ," % (rid, ai)
            s=s[0:-2]
            s+="\n"
            setid+=1

        return s

################ class Stata

    def UpdateSamples(self,a,set,RID=None): # run UpdateSet first
        if RID:
            self.rid_ai[RID]=[set,a]
            if self.rv[set].has_key(a):
                self.rv[set][a]+=[RID]
            else: 
                self.rv[set][a]=[RID]
        else:
            self.rv[set][a]=[self.RID]
            self.rid_ai[self.RID]=[set,a]
            self.RID+=1

################ class Stata

    def UpdateSet(self):
        self.rv.append({})

################ class Stata

    def sample(self,set):
        return self.rv[set]

################ class Stata

    def SampleMean(self,set):
        if self.SampMean.has_key(set):
            return self.SampMean[set]
        else:
            mean = CalcSampMean(set,self.rv[set].keys())
            self.SampMean[set]=mean
            return mean

################ class Stata

    def SampleVariance(self,set):
        if self.SampVar.has_key(set):
            return self.SampVar[set]
        else:
            sampmean=self.SampleMean(set)
            var=CalcSampVar(set,self.rv[set].keys(),sampmean)
            self.SampVar[set]=var
            return var

################ class Stata

    def CalcProb(self,set,x):
        '''
        To cacluate the probability P(X<=x). I.e., number of samples in rv that
        are less than or equal to x divided by the total number of samples.
        '''
        return CalcProb(set,x,self.rv[set].keys())

################ class Stata

    def GetRID(self,ai):
        '''
        loop through the sets and return the first value for a given ai
        '''
        for set in self.rv:
            if set.has_key(ai):
                return set[ai]

        return -12

################ end class Stata

class StatN:
    '''
    This class stores and operates on sets of samples of a random
    variable.  Written in support of DamModel to be stored with
    DamModel.__DamOro[doid].
    '''

    def __init__(self,sets,samples):
        self.rv=[] # to store all samples of the random variable
        self.SampMean={} # to store the sample mean of each set, setid:mean
        self.SampVar={} # to store the sample variance of each set, setid:var
        self.sets=sets # integer number of sets
        # samples = integer number of samples or str(dynamic) to instruct StatN
        # to set self.samples dynamically as you update samples
        if samples == 'dynamic':
            self.samples=0
            self.flag = 'dynamic'
        else: 
            self.samples=samples
            self.flag = 'static'

################ class StatN

    def __repr__(self):
        s = ""
        setid=0
        for set in self.rv:
            s+="Set %s" % (setid)
            s+="\n"
            iters=set.values()
            iters.sort()
            for ai in iters:
                s+="%1.4e ," %(ai)
            s=s[0:-2]
            s+="\n"
            setid+=1

        return s

################ class StatN

    def UpdateSamples(self,N,set,RID): # run UpdateSet first
        if type(RID)==type([]):
            for rid in RID:
                if self.rv[set].has_key(rid): continue
                else: self.rv[set][rid]=N
        else: # RID is not a list (like from twins)
            self.rv[set][RID]=N
        if self.flag=='dynamic':
            self.samples+=1

################ class StatN

    def UpdateSet(self):
        self.rv.append({})

################ class StatN

    def sample(self,set):
        return self.rv[set]

################ class StatN

    def SampleMean(self,set):
        if self.SampMean.has_key(set):
            return self.SampMean[set]
        else:
            mean = CalcSampMean(set,self.rv[set].values())
            self.SampMean[set]=mean
            return mean

################ class StatN

    def SampleVariance(self,set):
        if self.SampVar.has_key(set):
            return self.SampVar[set]
        else:
            sampmean=self.SampleMean(set)
            var=CalcSampVar(set,self.rv[set].values(),sampmean)
            self.SampVar[set]=var
            return var

################ class StatN

    def CalcProb(self,set,x):
        '''
        To cacluate the probability P(X<=x). I.e., number of samples in rv that
        are less than or equal to x divided by the total number of samples.
        '''
        return CalcProb(set,x,self.rv[set].values())

################ end class StatN

#------------------------------------------------------
# Test Driver
#------------------------------------------------------

if __name__ == "__main__":

    TestStata()
    TestStatN()

