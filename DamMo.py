
import os, sys, math, numpy as Numeric, cPickle, types

# wash
import Vec3D
import ColTensor
import MeshTools
##import dadN
import JohnsVectorTools as JVT
# me
import Statistic, DamClass, VarAmplitude, DamErrors
# John D
import GeomUtils
import Integration 

# give dumpfile the string name for a file and varamp will dump
# a vs N to a file named dumpfile.  Leave as None if you do not
# want it to DumpFile().  
dumpfile = 0
##dumpfile = "h:\\users\jme32\\research\\URETI\\Verification\\119032_will_-a"
Nlimit = 55000.0
Nout = 1

# Some global functions...

def holdit():
    '''
    a little function for debugging.
    '''
    print ' '
    print ' ################################################'
    print " ---> hit enter to close figure and move on <---"
    print ' ################################################'
    print ' '
    c = sys.stdin.read(1)

def DumpFile(xlist,ylist,filename):
    '''
    function to dump items in list out to a text named filename.
    '''

    fil=open(filename,'a')
    for i in range(min(len(xlist),len(ylist))):
        fil.write(str(xlist[i])+' '+str(ylist[i])+"\n")

    fil.write("# *** \n")
    fil.close()

class DamModel:

    #------------------------------------------------------
    # Embedded Class
    #------------------------------------------------------

    class __DamOroContainer:

        def __init__(self,sets,samples):
            self.ai=0.0 # store the initial crack size that gotcha! 
            self.DamEl=[] # list of damage elements for this doid
            self.DamElLength=0 # stored length of DamEl list changes on
               # subsequent DamMo.AddFDam() so i don't store every ai in
               # variable amplitude loading
            self.nextdN = 0.0 # for RK5 scheme. reset in SimDamGrowth

            # initialize to something that doesn't make sense (i.e. willgrow
            # should be -1 - compressive stress field, 0 - subcritical,
            # 1 - will grow, 2 - unstable, 3 - change shape, 4 - outgrown
            # surroundings.  Will be set to appropriate value in SimDamGrowth)
            self.WillGrow = 10 
            self.Life=-1

            # use a Statistic.Stat class to store and operate on sets of
            # samples of random variables.
            self.N = Statistic.StatN(sets,samples)

            # self.ProbFail is filled by DamModel.CalcStats() as:
            # [[Pf11, Pf12, Pf13...], [Pf21, Pf22, Pf23...],...]
            # where the first set of probs corresponds to the probability of
            # failure of the first set of samples and
            # P(N<=ncr1) = Pf11
            # whereas the prob of failure for the second set of samples is:
            # P(N<=ncr1) = Pf21
            self.ProbFail = [] 

            # 'check' = 0.0 - initially 0.0 self.N is not populated.  
            #           Switch to 1.0 when DamModel.UpdateSample is first
            #           called. 
            self.check = 0.0

        def Flush(self):
            '''
            reset some of the parameters for monte carlo variable amplitude
            problems
            '''
            self.ai=0.0 
            self.DamEl=[] 
            self.DamElLength=0 
            self.WillGrow = 10 
            self.Life=-1

######## DamModel

    def AddDamOro(self,doid,ais=None):
        self.DamOro[doid] = self.__DamOroContainer(self.parameters.sets, \
                                                   self.parameters.samples)
        if ais: self.ais=ais # so can change from one doid to the next... 

######## DamModel

    def AddFDam(self,doid,xyz,a,b,material,verbose,verify):
        if len(self.DamOro[doid].DamEl) != 0:
            self.DamOro[doid].DamElLength=len(self.DamOro[doid].DamEl)
        self.DamOro[doid].DamEl+=[DamClass.Fellipse(doid,xyz,self.model,a,b, \
                                                    material,None,verbose, \
                                                    verify)]

######## DamModel

    def BuildDkva(self,doid,did,material,errfile_name,r,iters,Scale=None):
        '''
        Builds the Delta K vs. a curve to be integrated for life prediction.
        Given a did

        (Similar method as SimDamGrowth)
        '''

        if not Scale:
            Scale = 1.0
        Flag=0;
##        while Flag == 0:
        # counting to iters will hopefully get us a will_grow == 2!
        for count in xrange(iters):
            # First answer the question, is this an interesting case?  Also,
            # WillGrow should return change of shape info as necessary.

#######     KEEP THIS CODE for parallel jobs...   #############
#######     |||||||||||||||||||||||||||||||||||   #############
#######     VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV   #############
            try:
                will_grow,dam_elem,have_doubt= \
                    self.__dam_elems[did][-1].WillGrow(self.model,material, \
                                                       0.0,self.CMesh,Scale)
            except KeyboardInterrupt:
                raise 'KeyboardInterrupt'
            except:
                will_grow=2
                # remove the last stuff in DamHistory because we couldn't
                # calc. K to go with it.
                self.__dam_elems[did][-1].PopLast()
                # retrieve the initial a
                a_in = self.__dam_elems[did][0].DamHistory['ab'][0]
                errfile=open(errfile_name,'a')
                errfile.write(str(doid)+' '+str(a_in)+"\n")
                errfile.close()

#######     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   #############
#######     |||||||||||||||||||||||||||||||||||   #############
#######     KEEP THIS CODE for parallel jobs...   #############

#######     KEEP THIS CODE FOR DEBUG PURPOSES...   #############
#######     |||||||||||||||||||||||||||||||||||    #############
#######     VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV    #############
##            will_grow,dam_elem,have_doubt= \
##                    self.__dam_elems[did][-1].WillGrow(self.model,material,\
##                                                       0.0, self.CMesh,Scale)
##
##            if verification:
####                print will_grow
##                af,bf = self.GiveCurrentGeo(did)[1][0],\
##                        self.GiveCurrentGeo(did)[1][1]
##                if will_grow == 3:
##                    print ' --> switch at N=', N, af, bf
##                else: print doid, will_grow, \
##                      self.__dam_elems[did][-1].DamHistory

##            self.__dam_elems[did][-1].PrintInfo(did)
##            print ' (me - DadModel 361) SimDamGrowth counter:',self.Count
##            print ' (me - DadModel 362)', will_grow,dam_elem,have_doubt
##            print ' (me - DadModel 363) doid = ', doid
#######     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    #############
#######     |||||||||||||||||||||||||||||||||||    #############
#######     KEEP THIS CODE FOR DEBUG PURPOSES...   #############

            # will_grow = -1 if Ki <= 0.0, then stop growing crack
            if will_grow == -1: 
                Flag=1
                self.__DamOro[doid].WillGrow = -1
                break

            # will_grow = 2 if Ki > Kic, then stop growing crack we keep 
            # growing in this case because for VarAmp may have less than 
            # sigma_max applied.  Means must check for Ki > Kic in Simpson.

            # will_grow = 4 --> crack has outgrown the extents of the body
            # no ligament left, crack has effectively caused net fracture
            elif will_grow == 4: 
                Flag=1
                self.__DamOro[doid].WillGrow = 4
                break

            # will_grow = 3 means change of shape
            elif will_grow == 3:
                self.__dam_elems[did][-1].PopLast()
                self.__dam_elems[did]+=[dam_elem]

            # for all other cases build the K v. a history.
            else:
                # to approximate the actual growth of the ellipse we want to
                # use the NASGRO equation to extend our ellipse.  However,
                # it is dependent on N.  Come up with some fake N here that
                # is related to "count", the counter in the for loop
                fakeN = count*100000/iters
                self.__dam_elems[did][-1].ExtendEllipse(material,self.model, \
                                                        r,fakeN)
                # to calculate K for the last set of crack sizes.  If crack is
                # too large at here, .CalcKi() dumps, so try it and if it
                # doesn't work, remove the last set of a's and b's from the
                # history.  
                if count == iters-1:
                    try: 
                        self.__dam_elems[did][-1].CalcKI()
                    except:
                        self.__dam_elems[did][-1].PopLast()

######## DamModel

    def GiveCurrentGeo(self,doid):
        as=[]
        for a in self.DamOro[doid].DamEl[-1].a:
            as+=[a[0][-1]]

        return as

######## DamModel

    def GrowDam(self,doid,DamEl,parameters,N,scale,integration): 
        '''
        Calculate individual growth amount for each "corner" of the ellipse
        and advance the crack.
        '''

        assert (len(DamEl.a[1])>0) # i.e there is a SIF History... 

        # Get current crack lengths...
        # Don't use __GiveCurrent because don't want average a and b...
        abab=[DamEl.a[i][0][-1] for i in range(DamEl.CrackFrontPoints)]

        # change to percentage of minimum ellipse dimension
        max_er=parameters.max_error*min(abab)
        # change back to absolute tolerance
        # max_er=max_error
        max_er=min(max_er,parameters.max_error)

        # time step:
        dN=parameters.dN

        if integration=='FWD':
            # use NASGRO eqn.
            # calc.dadN...
            rate=[DamEl.dadN[i].Calc_dadN(DamEl.a[i][1][-1],DamEl.a[i][0][-1],\
                                int(N)) for i in range(DamEl.CrackFrontPoints)]
##            blah=open('kVdadN.txt','a')
##            blah.write(str(DamEl.a[0][1][-1])+' '+str(rate[0])+"\n")
##            blah.close()

        elif integration=='RK4':
            # assume always using NASGRO equations with RK4 scheme...
            rate = Integration.RK4slope_vector(dN,N,abab, \
                                               DamEl.dAdN,[scale])

        else: # if self.IntegrationScheme=='RK5'
            # assume always using NASGRO equations with RK5 scheme...
            if (self.DamOro[doid].nextdN == 0.0): self.DamOro[doid].nextdN = dN

            # the integrations scheme can sometimes become unstable and pass
            # calcki numbers that don't make sense (too big).  Hence, we put in
            # a try and except block.  if RK5 scheme is bunk we do on small
            # step by the fwd euler method in hopes that springs us loose of
            # bad spot.  may need to add a finally or something to this in the
            # future.

            try:
                rate,dN,self.DamOro[doid].nextdN = \
                    Integration.RK5_CKslope_vector(self.DamOro[doid].nextdN,\
                    N,abab,DamEl.dAdN,[scale],max_er,.05)

            except DamErrors.FitPolyError, message:
                if self.verbose:
                    print ' switch to one step fwd Euler for', \
                          message[0],'at doid', doid, 'due to'
                    print '     ', message[1]
                rate = Integration.Eulerslope_vector( \
                    self.DamOro[doid].nextdN/1000.0,N,abab,DamEl.dAdN,[scale])
                self.DamOro[doid].nextdN = self.DamOro[doid].nextdN/1000.0

            except DamErrors.dAdNError, message:
                if self.verbose:
                    print ' switch to one step fwd Euler for', \
                          message[0],'at doid', doid, 'due to',message[1]
                rate = Integration.Eulerslope_vector( \
                    self.DamOro[doid].nextdN/1000.0,N,abab,DamEl.dAdN,[scale])
                self.DamOro[doid].nextdN = self.DamOro[doid].nextdN/1000.0

            except ValueError, message:
                if self.verbose:
                    print ' switch to one step fwd Euler for Fellipse', \
                          'at doid', doid, 'due to'
                    print '     ', message[0] 
                rate = Integration.Eulerslope_vector( \
                    self.DamOro[doid].nextdN/1000.0,N,abab,DamEl.dAdN,[scale])
                self.DamOro[doid].nextdN = self.DamOro[doid].nextdN/1000.0

        # Caclulate the new a, b, a-, b-.
        # a, b and neg_a, neg_b never to take negative value (which shouldn't
        # happen) so abs() is used to insure that here.
        inc=[abs(rate[i])*dN for i in range(len(abab))]
        big = max(inc)

        if integration=='RK5': # don't mess-up the adaptivity
            abab_new=[abab[i]+inc[i] for i in range(len(abab))]
        else: 
            # don't grow less than min_inc...
            if big < parameters.min_inc:
                dN=parameters.min_inc/max(rate)
                inc=[abs(rate[i])*dN for i in range(len(abab))]
                abab_new=[abab[i]+inc[i] for i in range(len(abab))]

            # don't grow more than a,b (average)
            elif max(inc) > \
                 min(DamEl.GiveCurrent()):
                dN=min(DamEl.GiveCurrent())/ \
                       max(rate)
                inc=[abs(rate[i])*dN for i in range(len(abab))]
                abab_new=[abab[i]+inc[i] for i in range(len(abab))]

            # otherwise, just do what seems natural... 
            else:
                abab_new=[abab[i]+inc[i] for i in range(len(abab))]

        # update Damage element
        DamEl.UpdateState(abab_new)

        return dN

######## DamModel

    def IntegrateLife(self,aL,did,doid,aR,N_max):
        '''
        using simpson's rule, calculate the life at did&doid for crack size
        spanning ai to af.

        ai - initial crack length
        af - final crack length
        r - increment factor (ai+1 = r * ai)
        '''

        # a tolerance to compare ai to af and determine if the interval is
        # effectively of zero length.
        tol = aL/1e6

        N=0.0 # start with 0 life.
        for dam_el in self.__dam_elems[did]:
            Nplus=dam_el.IntegrateLife(aL,aR,N_max,N,tol)
            if Nplus == -1:
                return -1
            N+=Nplus
        return N

######## DamModel

    def InterpolateLife(self,ai,doid,N_max):
        '''
        in case of monte carlo simulation with ONLY random ai (ie NOT random
        initial orientation, or material parameters), this is used to calculate
        the life for each ai based on the simulation results from the smallest
        ai (that will grow!).

        ai - initial crack radius
        did - damage i.d.

        return N

        N - the predicted life for the ai passed in

        works like this:

        Given:
        Ai = pi*ai*ai
        A_simulation  = [A0, A1, ..., Aj-1, Aj,..., An]
        dN_simulation = [0, dN1,..., dNj-1, dNj,..., dNn]

        A_simulation is a list of crack area as the crack grew.
        dN_simulation is a list of the delta N per growth cycle.

        compute dNi as:
        find j for Ai such that, Aj-1 < Ai < Aj
        dNi = (Aj - Ai)*(dNj)/(Aj - Aj-1)
        dNi += dNj+1 + ... + dNn

        this seems a little funny because the dNs are already Nj - Nj-1...
        '''

        As,types,as,dNs=self.InterpolateLifeLists(doid)
        

        # if self.DamOro[doid].WillGrow = 0, means the LARGEST crack stopped 
        # growing because DKi < DKth, so return N = -1 (similar to setting
        # self.DamOro[#].Life = -1), or max life for doid,
        # self.DamOro[doid].Life, or perhaps even self.HighLife... ?????
        if self.DamOro[doid].WillGrow == 0 or \
           self.DamOro[doid].WillGrow == -1: return N_max*1.01

        # to determine the appropriate damage type and calc crack area:
        j = -1
        for i in range(len(as)):
            if as[i] > ai:
                j = i
                break

        try: 
            if j == 0: type = types[0]
            else: type = types[j-1]
        except IndexError:
            print j, doid, ai, self.DamOro[doid].WillGrow
            print self.DamOro[doid].DamEl[-1]
            print as
            raise

        if type == 0:
            Ai = math.pi*ai*ai
        elif type == 1:
            Ai = 0.5*math.pi*ai*ai
        else: Ai = 0.25*math.pi*ai*ai

        # find j again by comparing crack areas 
        j = -1
        for i in range(len(as)):
            if As[i] > Ai:
                j = i
                break

        # calculate life
        if j > 0:
            N=Numeric.sum(dNs[j+1:])
            N+=(As[j]-Ai)*dNs[j]/(As[j]-As[j-1])
            if N > N_max:
                N=N_max*1.01

        # if j == 0, means ai passed in is smaller than the initial crack the 
        # simulation was run for (because the ai passed in must have been too 
        # small to end in DKi > DKc).  This means that when the during the
        # simulation for ai passed in, the code exited with will_grow = 0.
        # So return N = 1.01/*N_max for this ai (1% > N_max to be consistent
        # with .SimDamGrowth.
        elif j == 0: N=1.01*N_max
        else: N=0

        return N

######## DamModel

    def InterpolateLifeLists(self,doid):
        As=[]; types=[]; as=[]; dNs=[]
        for DamEl in self.DamOro[doid].DamEl:
            DamEl.AppendInterpolateLifeLists(As,types,as,dNs)
        return As,types,as,dNs

######## DamModel

    def NFile(self,doid,nfile):
        '''
        writes info to nfile
        '''
        stuff=[]
        if self.parameters.monte:
            for set in self.DamOro[doid].N.rv:
                for key in set.keys(): # key is RID
                    if set[key] < self.parameters.N_max:
                        stuff+=[[doid,key,int(set[key])]]
        else: stuff+=[[doid,0,int(self.DamOro[doid].Life)]]

        if len(stuff)>0: self.ToFile(nfile,stuff)
        else:
            stuff+=[[doid,-1,int(self.parameters.N_max)]]
            self.ToFile(nfile,stuff)

######## DamModel

    def PrintDamInfo(self,doid,dam,monte,set):
        '''
        Prints information about DamModel to screen.
        '''

        dashes='---------------------------------------------'
        dashes+='------------------'
        dashplus='---------+-----------------+-----------------+'
        dashplus+='-----------------'

        if doid == 'all':
            for i in self.DamOro.keys():
                xyz,delxyz,sigxyz=self.model.GetNodeInfo(i)
                print ''
                print '---------------- Damage Origin id = %8i' % (i), \
                      '------------------'
                print ' coords         = ', xyz
                print ' initial a      = ', self.DamOro[i].ai
                print ' Predicted Life = ', self.DamOro[i].Life

                if dam == 'all':
                    print ' **** Damage Elements: '
                    for damel in self.DamOro[i].DamEl:
                        print damel

                if monte:
                    print ' '
                    # initial area = initial area of full ellipse! 
                    print '     set |    initial a    |  initial area   |',\
                    '  predicted life'
                    print dashplus
                    if set == 'all':
                        for ij in range(self.parameters.sets):
                            init_a_set=self.ais.rv[ij].keys()
                            init_a_set.sort()
                            for ai in init_a_set:
                                RID=self.ais.rv[ij][ai]
                                for rid in RID:
                                    print '%8i | %15.6e | %15.6e |   %.0f' % \
                                          (ij, \
                                          ai,(ai**2)*math.pi, \
                                          self.DamOro[i].N.rv[ij][rid])
                            s="  Sample Mean = %10i " % \
                               (int(self.DamOro[i].N.SampleMean(ij)))
                            s+="|  Sample Stnd Dev = %10i" % \
                                  (int(math.sqrt(\
                                      self.DamOro[i].N.SampleVariance(ij))))
                            print s
                            print dashplus
                    else:
                        init_a_set=self.ais.rv[set].keys()
                        init_a_set.sort()
                        for ai in init_a_set:
                            RID=self.ais.rv[set][ai]
                            print '%8i | %15.6e | %15.6e |   %.0f' % (ij, \
                                  ai,(ai**2)*math.pi, \
                                  self.DamOro[i].N.rv[set][RID])
                        print "         Sample Mean = ", \
                              int(self.DamOro[i].N.SampleMean(set))
                        print dashplus
                print dashes,"\n"

        elif doid != 'none':
            print ''
            print '---------------- Damage Origin id = %8i' % (doid), \
                  '------------------'
            print ' coords         = ', self.__DamOro[doid].coords
            print ' dam_list       = ', self.__DamOro[doid].dam_list
            print ' initial a      = ', self.__DamOro[doid].ai
            print ' Predicted Life = ', self.__DamOro[doid].Life
            print dashes,"\n"

######## DamModel

    def ReturnDK(self,aa,a,K):
        '''
        CURRENTLY UNUSED CODE!! (11/4/06)
        return the appropriate list of K's that match aa, linear interpolation.

        return None if an element of aa exceeds the largest a, we have outgrone
        the "critical" a.  
        '''

        # if crack has outgrown the SIF history, return None.  Only check this
        # for a & b (i.e. not -a, -b) because the crack can transition and
        # still be stable.  
        if aa[0] > a[0][-1]: return None
        if aa[1] > a[1][-1]: return None

        # get the indices corresponding to the first a greater than aa
        jj=self.Greater(a,aa)
        lenjj=len(jj)

        # if jj=[None,None,None,None], return None...
        if max(jj) == None: return None

        # negative ii's are accounted for below... for i in range(len(M)):...
        ii = [jj[k]-1 for k in range(lenjj)] # JVT.Minus(jj,[1,1,1,1])

        # get the bounding K's
        Kjj=[K[j][jj[j]] for j in range(lenjj)]
        Kii=[K[i][ii[i]] for i in range(lenjj)]

        # get the bounding a's
        ajj=[a[j][jj[j]] for j in range(lenjj)]
        aii=[a[i][ii[i]] for i in range(lenjj)]

        # compute the slope for linear interpolation and DK
        M=JVT.Divide(JVT.Minus(Kjj,Kii),JVT.Minus(ajj,aii))
        DK=JVT.Plus(Kii,JVT.Star(M,JVT.Minus(aa,aii)))

        # ** aa[i] SMALLER than a[i][0], which can happen due to transition,
        # extrapolate... this is not very eloquent, but shouldn't happen often
        for im in range(lenjj):
            if ii[im] < 0.0:
                # it is possible that a[im][jj[im]+1] = ajj[im] if a NASGRO
                # growth rate is zero.  If so, assign dk = K[im][jj[im]+1]
                if a[im][jj[im]+1]-ajj[im] <= 0.0: dk = K[im][jj[im]+1]
                else:
                    m=(K[im][jj[im]+1]-Kjj[im])/ \
                      (a[im][jj[im]+1]-ajj[im])
                    dk = Kjj[im]-m*(ajj[im]-aa[im])

                # is possible that when extrapolate, value falls below zero
                # in this case, just use the smallest K in the list
                if dk < 0.0: DK[im] = K[im][0]
                else: DK[im]=dk

##        # if there are some negative M's there must be some ii's that are
##        # negative, which means aa[i] is smaller than a[i][0].  If this occurs
##        # assign DK[i]=K[i][0]
##        for i in range(len(M)):
##            if M[i] < 0.0:
##                DK[i]=K[i][0]

        # ** aa[i] LARGER than largest a[i]
        # This should only be possible for i=2,3 because of ifs at beginning of
        # this method.  What happens to cause this is if the crack transitions,
        # say from a semi to a quarter, the SIF history for the -a location is
        # shorter than for +a and +b.  Assign DK = 0 
        try:
            if aa[3] > a[3][-1]: DK[3]=0.0
        except IndexError:
            try:
                if aa[2] > a[2][-1]: DK[2]=0.0
            except IndexError: pass

        return DK

######## DamModel

    def SIFFile(self,doid,sfile):
        stuff=[]
        ai=self.DamOro[doid].DamEl[0].a[0][0][0]
        RID=self.ais.GetRID(ai)
        for damel in self.DamOro[doid].DamEl:
            for i in range(len(damel.a[0][1])):
                this_stuff=[]
                if damel.names[0]=='Fellipse':
                    this_stuff+=[doid,RID,damel.dN[i], \
                                 damel.a[0][0][i],damel.a[0][1][i], \
                                 damel.a[1][0][i],damel.a[1][1][i], \
                                 damel.a[2][0][i],damel.a[2][1][i], \
                                 damel.a[3][0][i],damel.a[3][1][i]]
                elif damel.names[0]=='Hellipse':
                    this_stuff+=[doid,RID,damel.dN[i], \
                                 damel.a[0][0][i],damel.a[0][1][i], \
                                 damel.a[1][0][i],damel.a[1][1][i], \
                                 damel.a[2][0][i],damel.a[2][1][i]]
                else: 
                    this_stuff+=[doid,RID,damel.dN[i], \
                                 damel.a[0][0][i],damel.a[0][1][i], \
                                 damel.a[1][0][i],damel.a[1][1][i]]
                stuff+=[this_stuff]
        self.ToFile(sfile,stuff)


######## DamModel

    def SimDamGrowth(self,doid,parameters,errfile_name,scale,R):
        '''
        Simulates the growth of damage until it is not interesting or it
        reaches some critical, life limiting condition.  Results in the life
        prediction resulting from damage at the given damage origin.

        *** Integration in "time" ***,
        '''

        # collect the corner and midside nodes of the surface mesh for use with
        # WillGrow in finding intersection of ellipses with surface facets.

        Flag=0; N=0.
        CurrentElement=self.DamOro[doid].DamEl[-1]
        self.DamOro[doid].ai=CurrentElement.a[0][0][0]
        while Flag == 0:
            # check crack shape...
            dam_elem,have_doubt= \
                    CurrentElement.GeometryCheck(self.CMesh)

            if dam_elem == 0:
                Ki = CurrentElement.CalculateKi(scale)
                will_grow = CurrentElement.WillGrow(N,R,Ki)
                # will_grow = -1 if Ki < 0.0, then stop growing crack
                if will_grow == -1: ## and monte == 0: #(see ** below.)
                    Flag=1
                    self.DamOro[doid].WillGrow = -1
                    dN=0.0

                # will_grow = 0 if Ki is less than Kth then stop growing crack
                elif will_grow == 0: ## and monte == 0: #(see ** below.)
                    Flag=1
                    self.DamOro[doid].WillGrow = 0
                    dN=0.0

                # Parameter will_grow is set to 1 in WillGrow if the crack,
                # indeed, will grow stably.  (ie. Kth<Ki<Kic)
                elif will_grow == 1:
                    dN=self.GrowDam(doid,CurrentElement,parameters,N,scale, \
                                    self.IntegrationScheme)
                    N+=dN

                    # if N has reached a practical maximum... 
                    if N >= parameters.N_max:
                        Flag=1
                        self.DamOro[doid].WillGrow = 1

                # will_grow = 2 means unstable growth (ie. Ki>=Kic)
                elif will_grow == 2:
                    Flag=1
                    self.DamOro[doid].WillGrow = 2
                    dN=0.0

                # will_grow = 4 crack grew outside body, assume net fracture
                elif will_grow == 4:
                    Flag=1
                    self.DamOro[doid].WillGrow = 4
                    dN=0.0

                # Update dN
                CurrentElement.dN+=[dN]

            elif dam_elem == 4:
                # Update dN
                CurrentElement.dN+=[0.0]
                Flag=1
                self.DamOro[doid].WillGrow = 4

            else:
                # Update dN
                CurrentElement.dN+=[0.0]
                self.DamOro[doid].DamEl+=[dam_elem]
                CurrentElement=dam_elem

        #end while

        # if will_grow = 0 or 2, compute life at DamOro and 
        # set self.DamOro[doid].Life and update HighLife and LowLife
        if N > 0.0:

            if N > parameters.N_max:
                N = parameters.N_max*1.01
                self.DamOro[doid].Life=N

            if N > self.HighLife:
                self.HighLife = N
                self.DamOro[doid].Life=N

            if N < self.LowLife and will_grow != 0:
                self.LowLife = N
                self.DamOro[doid].Life=N

        # if N is less than one then the crack was never grown.  this could
        # happened for two reasons.  1) Ki was never larger than Kth or 2)
        # it was unstable at initial size.  
        else:
            # if 1) set to -1.0 and deal with it later when drawing contour
            # probably by setting it equal to the highest life prediction!?
            # yah, that's what i do. 
            if will_grow == 0 or \
               will_grow == -1: self.DamOro[doid].Life=parameters.N_max*1.01

            # if 2) predicted life should be set to zero
            else:
                self.DamOro[doid].Life = 0.0
                self.LowLife = 0.0

##        DumpFile(self.DamOro[doid].DamEl[-1].a[0][0], \
##                 self.DamOro[doid].DamEl[-1].a[0][1],'Kva.txt')
##        DumpFile(self.DamOro[doid].DamEl[-1].a[1][0], \
##                 self.DamOro[doid].DamEl[-1].a[1][1],'Kvb.txt')

        return N

######## DamModel

    def ToFile(self,f,stuff):
        '''
        Saves stuff to text files to be loaded into SQL:

        arguements:
        f - file object to write to
        stuff - a list of lists of the stuff to be written to the file
        '''

        for line in stuff:
            for item in line:
                f.write(str(item)+' ')
            f.write("\n")

######## DamModel

    def ToMAPFile(self,MAPFile,set):
        '''
        Prints file for contouring in MAP

        MAPFile - file object to write to
        set  - the set from which to write mean value for
        '''

        DamKeys=self.DamOro.keys()
        DamKeys.sort()

        MAPFile.write('LIFE 0'+"\n")
        for i in DamKeys: # i = doid
            if self.parameters.monte:
                if self.DamOro[i].N.SampleMean(set) == -1.0:
                    MAPFile.write(str(i)+' '+str(HighLife)+"\n")
                else:
                    MAPFile.write(str(i)+' '+ \
                                  str(self.DamOro[i].N.SampleMean(set))+"\n")
            else:
                if self.DamOro[i].Life == -1.0:
                    MAPFile.write(str(i)+' '+str(HighLife)+"\n")
                else:
                    MAPFile.write(str(i)+' '+str(self.DamOro[i].Life)+"\n")
        MAPFile.close()

######## DamModel

    def NToPickle(self,filename,doid):
        '''
        Pickle the life predicitons.  First build a dictionary of the N's,
        then pickle.
        '''
        if doid=='all':
            dumpfile=open(filename,'w+b')
            damage_list={}
            for i in self.DamOro.keys(): # doid 
                damage_list[i] = self.DamOro[i].N
            cPickle.dump(damage_list,dumpfile,2)
            dumpfile.close()

        else:
            try:
                readfile=open(filename,'r+b')
                damage_list=cPickle.load(readfile)
                readfile.close()
                damage_list[doid]=self.DamOro[doid].N
                dumpfile=open(filename,'w+b')
                cPickle.dump(damage_list,dumpfile,2)
                dumpfile.close()
            except IOError:
                damage_list={doid:self.DamOro[doid].N}
                dumpfile=open(filename,'w+b')
                cPickle.dump(damage_list,dumpfile,2)
                dumpfile.close()

######## DamModel

    def __init__(self,model,node_list,verbose,verify,extension,saveall, \
                 parameters,DebugGeomUtils=False,SVIEW=False, \
                 integration='RK5'):
        self.DamOro = {} # keys = doid, values = __DamOroContainer instance
        self.node_list = node_list # ordered list of all node ids in vol. mesh
        self.LowLife = 1e+100 # fictitously high numba
        self.HighLife = -1.0 # Miller time!
        self.model = model # MeshTools object
        self.verbose = verbose # parameter for printing
        self.verify = verify
        self.IntegrationScheme = integration # 'FWD', 'RK4', 'RK5', or 'Simp'
        self.extension=extension # corresponds to MSTI_RANK or None
        self.saveall=saveall # boolean to save files as we progress
        self.parameters=parameters
        self.ais=None # pass in later as necessary... 

        # build CMesh object
        self.SurfaceNodesCoords = [] # coordinates of surface facets (no mids)
        self.SurfaceNodeIds = []
        self.SurfaceElements = {} # key is seid, values are corner nids
        self.__MakeCMeshObj()
        if DebugGeomUtils: self.__WriteGeomUtilsDebugFile(DebugGeomUtils)
        if SVIEW: self.__WriteSVIEWFile(SVIEW)

######## DamModel

    def __MakeCMeshObj(self):
        # loop over all nodes to populate self.SurfaceElements
        for nid in self.node_list:
            if self.model.IsSurfaceNode(nid) == 1:
                self.SurfaceNodeIds += [nid]
                list_seids = self.model.GetAdjacentSurfElems(nid)
                for seid in list_seids:
                    if self.SurfaceElements.has_key(seid):
                        continue
                    else:
                        SurfaceNidList = self.model.GetSurfElemInfo(seid)
                        if len(SurfaceNidList) == 6: # triangular surface facet
                            self.SurfaceElements[seid] = ((SurfaceNidList[0], \
                                                           SurfaceNidList[1], \
                                                           SurfaceNidList[2]),\
                                                          (SurfaceNidList[3], \
                                                           SurfaceNidList[4], \
                                                           SurfaceNidList[5]))
                        elif len(SurfaceNidList) == 8: # quad surface facet
                            self.SurfaceElements[seid] = ((SurfaceNidList[0], \
                                                           SurfaceNidList[1], \
                                                           SurfaceNidList[2], \
                                                           SurfaceNidList[3]),\
                                                          (SurfaceNidList[4], \
                                                           SurfaceNidList[5], \
                                                           SurfaceNidList[6], \
                                                           SurfaceNidList[7]))
                        elif len(SurfaceNidList) == 3: # linear tri surf facet
                            self.SurfaceElements[seid] = ((SurfaceNidList[0], \
                                                           SurfaceNidList[1], \
                                                           SurfaceNidList[2]),\
                                                          (None))
                        else: # linear quadrilateral surface facet
                            self.SurfaceElements[seid] = ((SurfaceNidList[0], \
                                                           SurfaceNidList[1], \
                                                           SurfaceNidList[2], \
                                                           SurfaceNidList[3]),\
                                                          (None))
        # build self.SurfaceNodesCoords
        for key in self.SurfaceElements:
            self.SurfaceNodesCoords.append([])
            SurfaceNidList = self.SurfaceElements[key][0]
            for snid in SurfaceNidList: # snid = surface node id
                try: 
                    self.SurfaceNodesCoords[-1] += \
                                            [self.model.GetNodeInfo(snid)[0]]
                except MeshTools.InvalidNodeId:
                    print snid
                    raise

        # object for GeomUtils used in DamClass.Fellipse.__BuildPhi
        self.CMesh = GeomUtils.BuildSurfMeshCObject(self.SurfaceNodesCoords)

######## DamModel

    def __WriteGeomUtilsDebugFile(self,DebugGeomUtils):
        # code to write to file so can copy & paste into f00*.py to
        # debug GeomUtils (comment out most of time).
        blah=open(DebugGeomUtils,'w')
        blah.write('msh=[\\'+"\n")
        for list in self.SurfaceNodesCoords:
            if len(list)==3:
                txt='[Vec3D.Vec3D( '
                txt+=str(list[0].x())+','+str(list[0].y())+','+ \
                      str(list[0].z())
                txt+='),Vec3D.Vec3D('
                txt+=str(list[1].x())+','+str(list[1].y())+','+ \
                      str(list[1].z())
                txt+='),Vec3D.Vec3D('
                txt+=str(list[2].x())+','+str(list[2].y())+','+ \
                      str(list[2].z())
                txt+=')],\\'+"\n"
                blah.write(txt)
            else:
                txt='[Vec3D.Vec3D( '
                txt+=str(list[0].x())+','+str(list[0].y())+','+ \
                      str(list[0].z())
                txt+='),Vec3D.Vec3D('
                txt+=str(list[1].x())+','+str(list[1].y())+','+ \
                      str(list[1].z())
                txt+='),Vec3D.Vec3D('
                txt+=str(list[2].x())+','+str(list[2].y())+','+ \
                      str(list[2].z())
                txt+='),Vec3D.Vec3D('
                txt+=str(list[3].x())+','+str(list[3].y())+','+ \
                      str(list[3].z())
                txt+=')],\\'+"\n"
                blah.write(txt)
        blah.write(']')
        blah.close()
        # exit without further processing...
        sys.exit()

######## DamModel

    def __WriteSVIEWFile(self,SviewFile):
        '''
        code to write surface mesh to svw file so i can load with an old
        version of the MAP that displays ellipses.  This is pretty much just
        for me.
        '''
        file = open(SviewFile,'w')
        for facet in self.SurfaceNodesCoords:
            if len(facet)==3:
                file.write('p'+' 3 '+str(facet[0].x())+' '+ \
                                     str(facet[0].y())+' '+ \
                                     str(facet[0].z())+' '+ \
                                     str(facet[1].x())+' '+ \
                                     str(facet[1].y())+' '+ \
                                     str(facet[1].z())+' '+ \
                                     str(facet[2].x())+' '+ \
                                     str(facet[2].y())+' '+ \
                                     str(facet[2].z())+' '+ \
                                     '0.58 0.58 1.0'+"\n")
            else:
                file.write('p'+' 4 '+str(facet[0].x())+' '+ \
                                     str(facet[0].y())+' '+ \
                                     str(facet[0].z())+' '+ \
                                     str(facet[1].x())+' '+ \
                                     str(facet[1].y())+' '+ \
                                     str(facet[1].z())+' '+ \
                                     str(facet[2].x())+' '+ \
                                     str(facet[2].y())+' '+ \
                                     str(facet[2].z())+' '+ \
                                     str(facet[3].x())+' '+ \
                                     str(facet[3].y())+' '+ \
                                     str(facet[3].z())+' '+ \
                                     '0.58 0.58 1.0'+"\n")
        file.close()

######## DamModel

    def UpdateSample(self,doid,a,N,set):
        '''
        method to update the __DamOroInfo.statistics the sample information in
        the instance of the Statistic.Stat class.

        a - initial crack size
        N - associated life prediction
        set - set id
        '''

        if self.DamOro[doid].check == 0.0:
            self.DamOro[doid].check = 1.0

        RID=self.ais.rv[set][a]
        self.DamOro[doid].N.UpdateSamples(N,set,RID)

######## DamModel

    def UpdateSet(self,doid):
        '''
        Originally, there were more than one Statistic class associated with a
        damage origin.  This method allowed one call from DDSim.py to update
        them all.  It is not currently necessary, i could make this call
        directly in DDSim.py (cracks.DamOro[doid].N.UpdateSet()), but i want
        to keep this method, incase things change.  
        '''

        self.DamOro[doid].N.UpdateSet()

######## DamModel

    def CalcStats(self,doid,set,ncr):
        '''
        Calculate mean, variance and prob. of failure for sets of samples
        generated by the Monte Carlo simulation and stored in __DamOro.

        ncr - a list of critical life values
        '''

        self.DamOro[doid].N.SampleMean(set)
        self.DamOro[doid].N.SampleVariance(set)

        # assuming CalcStats is called sequentially for set = 0, 1, 2, etc.
        # must update self.__DamOro[doid].ProbFail for new info... 
        self.DamOro[doid].ProbFail.append([])

        for i in ncr:
            self.DamOro[doid].ProbFail[set] += \
                        [self.DamOro[doid].N.CalcProb(set,i)]

######## DamModel

    def FrontPoints(self,frontfile):
        '''
        i only want to run this for one doid runs, so pick the first damage
        element in DamOro.
        '''
        keys=self.DamOro.keys()
        keys.sort()

        for Damage in self.DamOro[keys[0]].DamEl:
            for i in range(len(Damage.a[0][0])):
                qpnts=Damage.ComputeFrontPoints(i)
                for qpnt in qpnts:
                    frontfile.write(str(qpnt.x())+" "+str(qpnt.y())+" "+ \
                                    str(qpnt.z())+"\n")
                frontfile.write('################# \n')


######## DamModel

    def ComputeEulerAngles(self,rot):
        '''
        compute the three euler angles equivalent to the rotation matrix, rot.

        From C++ code borrowed from Wash (via JD)
        '''

        # first of all, let's permute rot so that it mathces JD's definition of
        # the rotation matrix (recall, local z is orthog. to ellipse for him)
        # AND transpose it so columns are eig. vectors (to match wash's code
        # i'm copying for the euler angles).  
        rot=[[rot[1][0],rot[2][0],rot[0][0]], \
             [rot[1][1],rot[2][1],rot[0][1]], \
             [rot[1][2],rot[2][2],rot[0][2]]]

        # make sure there are no entries in rot > 1.0
        for i in range(len(rot)):
            for j in range(len(rot[i])):
                if rot[i][j] > 1.0: rot[i][j] = 1.0
                if rot[i][j] < -1.0: rot[i][j] = -1.0

        # possible thetas:  t0 = theta_z; t1 = theta_x; t2 = theta_y
        t0=[]; t1=[]; t2=[]

        # compute two possible values of theta_x (t1)
        if rot[2][1] == 1.0:
            t1 += [math.pi/2.0]
            t1 += [t1[0]]
        else:
            t1 += [math.asin(rot[2][1])]
            if t1[0] > 0.0: t1 += [math.pi - t1[0]]
            else: t1 += [-math.pi - t1[0]]

        if t1[0] != t1[1]:
            # compute 4 possible values for theta_z (t0)
            tmp = rot[1][1]/math.cos(t1[0])
            if tmp >  1.0: tmp =  1.0
            if tmp < -1.0: tmp = -1.0
            t0 += [math.acos(tmp)]
            t0 += [-t0[0]]
            tmp = rot[1][1]/math.cos(t1[1])
            if tmp >  1.0: tmp =  1.0
            if tmp < -1.0: tmp = -1.0
            t0 += [math.acos(tmp)]
            t0 += [-t0[2]]

            # compute 4 possible values of theta_y (t2)
            tmp = rot[2][2]/math.cos(t1[0])
            if tmp >  1.0: tmp =  1.0
            if tmp < -1.0: tmp = -1.0
            t2 += [math.acos(tmp)]
            t2 += [-t2[0]]

            tmp = rot[2][2]/math.cos(t1[1])
            if tmp >  1.0: tmp =  1.0
            if tmp < -1.0: tmp = -1.0
            t2 += [math.acos(tmp)]
            t2 += [-t2[2]]
        else:
            t0 += [math.acos(rot[0][0])]
            t0 += [-t0[0]]
            t0 += [t0[0]]
            t0 += [-t0[0]]

            t2 += [0.0,0.0,0.0,0.0]

        # make tolerance and check
        tol = 0.0001
        groups = [[[1,1,1,1],[1,1,1,1]],\
                  [[1,1,1,1],[1,1,1,1]],\
                  [[1,1,1,1],[1,1,1,1]],\
                  [[1,1,1,1],[1,1,1,1]]]
        for i in range(len(groups)):
            for j in range(len(groups[i])):
                for k in range(len(groups[i][j])):

                    if not groups[i][j][k]: continue
                    val = math.cos(t0[i])*math.cos(t2[k]) - \
                          math.sin(t0[i])*math.sin(t1[j])*math.sin(t2[k])
                    if abs(val - rot[0][0]) > tol: groups[i][j][k] = 0

                    if not groups[i][j][k]: continue
                    val = math.sin(t0[i])*math.cos(t2[k]) + \
                          math.cos(t0[i])*math.sin(t1[j])*math.sin(t2[k])
                    if abs(val - rot[1][0]) > tol: groups[i][j][k] = 0

                    if not groups[i][j][k]: continue
                    val = -math.cos(t1[j])*math.sin(t2[k])
                    if abs(val - rot[2][0]) > tol: groups[i][j][k] = 0

                    if not groups[i][j][k]: continue
                    val = -math.sin(t0[i])*math.cos(t1[j])
                    if abs(val - rot[0][1]) > tol: groups[i][j][k] = 0

                    if not groups[i][j][k]: continue
                    val = math.cos(t0[i])*math.cos(t1[j])
                    if abs(val - rot[1][1]) > tol: groups[i][j][k] = 0

                    if not groups[i][j][k]: continue
                    val = math.sin(t1[j])
                    if abs(val - rot[2][1]) > tol: groups[i][j][k] = 0

                    if not groups[i][j][k]: continue
                    val = math.cos(t0[i])*math.sin(t2[k]) + \
                          math.sin(t0[i])*math.sin(t1[j])*math.cos(t2[k])
                    if abs(val - rot[0][2]) > tol: groups[i][j][k] = 0

                    if not groups[i][j][k]: continue
                    val = math.sin(t0[i])*math.sin(t2[k]) - \
                          math.cos(t0[i])*math.sin(t1[j])*math.cos(t2[k])
                    if abs(val - rot[1][2]) > tol: groups[i][j][k] = 0

                    if not groups[i][j][k]: continue
                    val = math.cos(t1[j])*math.cos(t2[k])
                    if abs(val - rot[2][2]) > tol: groups[i][j][k] = 0

        # at this point there should be two sets of possible angles, return the
        # first one (in radians)
        for i in range(len(groups)):
            for j in range(len(groups[i])):
                for k in range(len(groups[i][j])):
                    if groups[i][j][k]:
                        return [t1[j]*180.0/math.pi,\
                                t2[k]*180.0/math.pi,\
                                t0[i]*180.0/math.pi]

        return [0.0,0.0,0.0]

######## DamModel

    def Greater(self,a,aa):
        '''
        return a Numeric.array containing the indices of a larger than the
        corresponding entry in aa.
        '''

##        z=[Numeric.greater(a[k],aa[k]) for k in range(len(aa))]
##        j=[]
##        for zk in z:
##            zkk=zk.tolist()
##            try:
##                j+=[zkk.index(1)]
##            except ValueError:
##                j+=[0]

        j=[]
        for ii in range(len(a)):
            index=0
            try:
                while a[ii][index]<=aa[ii]: index += 1
                j+=[index]
            except IndexError: j+=[0]

        return j

######## DamModel

    def VarAmp(self,doid,aR,Nmax,spec,nore,scale,r):
        '''
        Function to do cycle-by-cycle integration for spectrum loading.  

        Use Kic for failure criteria since the NASGRO equation for Kc includes
        thickness, t, that will be hard to estimate here.  

        spec = VarAmplitude.Spectrum instance.
        nore = true, don't use willenborg, false (default) use it!
        scale = to scale the sig file
        r    = when to recompute K:  if a+ >= r*acurrent: recompute! 
        '''

        # reInitialize self.DamOro[doid].WillGrow = 10.
        # In the cycle-by-cycle while loop set to :
        # -1 - compressive stress field:  
        # 0 - subcritical: if through load program with no growth
        # 1 - will grow: while crack will still  grow
        # 2 - unstable: when DK > DKc
        # 3 - change shape: NEVER within VarAmp()
        # 4 - outgrown: net fracture
        self.DamOro[doid].WillGrow = 1

        # initialize some stuff
        spec.Initialize()
        DamEl=self.DamOro[doid].DamEl[-1]
        self.DamOro[doid].ai=DamEl.a[0][0][0]
        # Note:  the plastic zone size is kept track of within the damage
        # element. As yet there is no way to update this so if the damage
        # changes shape within the while loop below, the history at the crack
        # tip is lost.  5/27/06

        # do the cycle-by-cycle integration (vectorized)
        twice=0; count=0; # used looping through the sprectrum w/o growth
        N=0.0; Nold=0.0; inc=[]; Ks=[]

        while self.DamOro[doid].WillGrow == 1:
            # first check the current damage elements geometry... 
            if inc:
                # this is a little sloppy, RANGE is defined below... 
                test=[aCurrent[i]*r>aa[i] for i in RANGE]
                if min(test): pass
                else: damel,doubt=DamEl.GeometryCheck(self.CMesh)
            else: damel,doubt=DamEl.GeometryCheck(self.CMesh)
    
            # either compute K,da/dN, and new a's, exit or change damage shapes
            if damel == 0:
                lena=DamEl.CrackFrontPoints
                RANGE=range(lena)
                aa=[DamEl.a[i][0][-1] for i in RANGE]
                percentDK,R=spec.Delta()
                if inc:
                    if min(test):
                        Ks=[KsCurrent[i]*percentDK for i in RANGE]
                        wg=DamEl.WillGrow(N,R,Ks)

                    else:
                        aCurrent=aa
                        Ks=DamEl.CalculateKi(scale*percentDK)
                        wg=DamEl.WillGrow(N,R,Ks)
                        if wg != 4:
                            KsCurrent=[Ks[i]/percentDK for i in RANGE]

                else:
                    aCurrent=aa
                    Ks=DamEl.CalculateKi(scale*percentDK)
                    wg=DamEl.WillGrow(N,R,Ks)
                    if wg != 4:
                        KsCurrent=[Ks[i]/percentDK for i in RANGE]

                DamEl.dN+=[N-Nold]
                Nold=N
                N+=1.0

                # unstable crack growth
                if wg == 2:
##                    N+=1.0
                    self.DamOro[doid].WillGrow=2

                # net fracture
                elif wg == 4:
##                    N+=1.0
                    self.DamOro[doid].WillGrow=4

                # if it is a compressive cycle, increase N and continue
                # still don't have amplification for under load... 5/27/06
                elif wg == -1:
                    if twice == 0:
                        twice = 1
                        count +=1
##                        N+=1.0
                    elif count < spec.length:
                        count+=1
##                        N+=1.0
                    else:
##                        N+=1.0
                        self.DamOro[doid].WillGrow=-1

                else:
                    # compute the increment of growth
                    if nore:
                        inc=[DamEl.dadN[i].Calc_dadN(Ks[i],aa[i],int(N-1.0),R) \
                             for i in RANGE]
                    else: 
                        inc=[DamEl.dadN[i].Compute_dadN(Ks[i],aa[i],int(N-1.0),R) \
                             for i in RANGE]

                    aa = JVT.Plus(aa,inc)
                    # Update the damage element
                    for i in range(lena): DamEl.a[i][0]+=[aa[i]]

                    # Check if aa exceeds aR
                    if max(aa) > aR:
                        # only happens when called for Nminus, set back
                        # to old one
                        self.DamOro[doid].WillGrow = 1
##                        N+=1.0
                        break # 

                    # if N exceeds max N return N = -1 for processing at
                    # higher level
                    elif N >= Nmax:
                        N = -1
                        self.DamOro[doid].WillGrow = 0

                    # to prevent repetitive times through the load program
                    # with inc = 0.0...
                    elif max(inc)==0.0 and twice == 0:
                        twice = 1
                        count +=1
##                        N+=1.0
                    elif max(inc)==0.0:
                        if count < spec.length:
                            count+=1
##                            N+=1.0
                        else:
                            N = -1
                            self.DamOro[doid].WillGrow = 0
                    else:
                        twice = 0
                        count = 0
##                        N+=1.0

            elif damel == 4:
                self.DamOro[doid].WillGrow==4

            # if DK = None, stop growing, aa has exceeded a
            else:
                self.DamOro[doid].DamEl+=[damel]
                DamEl=damel
        # end while

        # if wg = 4 for the first time through for a new shape damage
        # element... 
        try:
            a,b,ka,kb=self.DamOro[doid].DamEl[-1].a[0][0][-1],\
                      self.DamOro[doid].DamEl[-1].a[1][0][-1],\
                      self.DamOro[doid].DamEl[-1].a[0][1][-1],\
                      self.DamOro[doid].DamEl[-1].a[1][1][-1]
        except:
            a,b,ka,kb=self.DamOro[doid].DamEl[-2].a[0][0][-1],\
                      self.DamOro[doid].DamEl[-2].a[1][0][-1],\
                      self.DamOro[doid].DamEl[-2].a[0][1][-1],\
                      self.DamOro[doid].DamEl[-2].a[1][1][-1]

##        DumpFile(self.DamOro[doid].DamEl[-1].a[0][0], \
##                 self.DamOro[doid].DamEl[-1].a[0][1],'Kva.txt')
##        DumpFile(self.DamOro[doid].DamEl[-1].a[1][0], \
##                 self.DamOro[doid].DamEl[-1].a[1][1],'Kvb.txt')

        return N,(a,b,ka,kb)

######## DamModel

    def WriteFinalAs(self,AfFile,doid):
        '''
        Method to write Final crack size information to a file.
        '''

        if doid=='all':
            AfFile.write('Final_a 0'+"\n")

            DamKeys=self.DamOro.keys()
            DamKeys.sort()

            for i in DamKeys:
                af,bf=self.DamOro[i].DamEl[-1].GiveCurrent()
                AfFile.write(str(i)+' '+str(af)+' '+str(bf)+"\n")

        else:
            af,bf=self.DamOro[doid].DamEl[-1].GiveCurrent()
            AfFile.write(str(doid)+' '+str(af)+' '+str(bf)+"\n")

######## DamModel

    def WriteInitialAs(self,AFile,Nmax,doid):
        '''
        Method to write initial crack size information to a file.
        '''

        if doid == 'all':
            AFile.write('Initial_a 0'+"\n")

            DamKeys=self.DamOro.keys()
            DamKeys.sort()

            for i in DamKeys:
                AFile.write(str(i)+' '+str(self.DamOro[i].ai)+"\n")
        else:
            AFile.write(str(doid)+' '+str(self.DamOro[doid].ai)+"\n")

######## DamModel

    def WriteRotations(self,RotFile,doid):
        '''
        Method to write the damage origin's orientation to a file to be read by
        MAP.
        '''

        if doid=='all':
            RotFile.write('Orient 1'+"\n")

            DamKeys=self.DamOro.keys()
            DamKeys.sort()

            for i in DamKeys:
                rot=self.DamOro[i].DamEl[0].Rotation
                euler=self.ComputeEulerAngles(rot)
                RotFile.write(str(i)+' '+str(euler[0])+' '+str(euler[1])+ \
                              ' '+str(euler[2])+"\n")
        else:
            rot=self.DamOro[doid].DamEl[0].Rotation
            euler=self.ComputeEulerAngles(rot)
            RotFile.write(str(doid)+' '+str(euler[0])+' '+str(euler[1])+ \
                          ' '+str(euler[2])+"\n")









