
import os, sys, math, pickle, types
import numpy as np

# wash
from . import Vec3D
from . import ColTensor
from . import MeshTools
##import dadN
from . import JohnsVectorTools as JVT
# me
from . import Statistic
from . import DamClass
from . import VarAmplitude
from . import DamErrors
# John D
from . import GeomUtils
from . import Integration
from . import exodus_io

# give dumpfile the string name for a file and varamp will dump
# a vs N to a file named dumpfile.  Leave as None if you do not
# want it to DumpFile().
dumpfile = 0
##dumpfile = "h:\\users\jme32\\research\\URETI\\Verification\\119032_will_-a"
Nlimit = 55000.0
Nout = 1

# Added 2026 (see docs/PORTING_NOTES.md and GrowDam's RK5 branch): a floor
# on how small a step GrowDam will ever take, as a fraction of the nominal
# parameters.dN, regardless of the max_growth_fraction-based cap or the
# exception fallback -- so neither can shrink the step all the way to zero.
_MIN_STEP_FRACTION = 0.001

# Some global functions...

def holdit():
    '''
    a little function for debugging.
    '''
    print(' ')
    print(' ################################################')
    print(" ---> hit enter to close figure and move on <---")
    print(' ################################################')
    print(' ')
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

class _StrippedDamEl:
    '''
    Picklable stand-in for a DamClass.Fellipse/Hellipse/Qellipse instance,
    used only after a doid's processing is fully complete and it's about to
    cross a process boundary (see ddsim.parallel). Every real DamEl holds a
    direct reference to the entire shared MeshTools model (self.model),
    which is what makes a live DamEl catastrophically expensive to pickle --
    but every *downstream* output method only ever reads two small facts off
    it: WriteRotations reads DamEl[0].Rotation, WriteFinalAs reads
    DamEl[-1].GiveCurrent(). This carries just those two, so
    DamOro[doid].DamEl can be safely replaced with a single-element list
    (index 0 and -1 alike) -- see _DamOroContainer.StripDamElForTransport()
    and docs/PORTING_NOTES.md. PrintDamInfo's optional full DamEl dump
    (dam=='all') goes empty for any doid that passed through this --
    accepted, screen-diagnostics only, not a file output.
    '''
    __slots__ = ('Rotation', '_af', '_bf')

    def __init__(self, rotation, af, bf):
        self.Rotation = rotation
        self._af = af
        self._bf = bf

    def GiveCurrent(self):
        return self._af, self._bf


class _DamOroContainer:
    '''
    Module-level (not nested in DamModel) so it's picklable: a class nested
    inside another class with a double-underscore name gets its *storage
    key* mangled (DamModel.__dict__['_DamModel__DamOroContainer']), but its
    __qualname__ stays the unmangled 'DamModel.__DamOroContainer' -- pickle
    resolves classes by walking __qualname__ via getattr from the module, so
    it looks for (and doesn't find) DamModel.__DamOroContainer and raises
    PicklingError. A single leading underscore at module level isn't
    mangled at any nesting depth, so this is picklable as long as DamEl
    doesn't hold anything that isn't (see _StrippedDamEl above).
    '''

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

    def StripDamElForTransport(self):
        '''
        Mutates self.DamEl down to a picklable single-element stand-in, for
        returning this doid's results across a ddsim.parallel worker process
        boundary. Call only after this doid's SimDamGrowth + all
        UpdateSample calls are done: UpdateSample resolves an initial crack
        size against the live, shared DamModel.ais at call time, but nothing
        after that point ever needs the live DamEl again.
        '''
        rotation = self.DamEl[0].Rotation
        af, bf = self.DamEl[-1].GiveCurrent()
        self.DamEl = [_StrippedDamEl(rotation, af, bf)]


class DamModel:

    #------------------------------------------------------
    # Embedded Class
    #------------------------------------------------------

######## DamModel

    def AddDamOro(self,doid,ais=None):
        self.DamOro[doid] = _DamOroContainer(self.parameters.sets, \
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
        for count in range(iters):
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
                raise
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
        alist=[]
        for a in self.DamOro[doid].DamEl[-1].a:
            alist+=[a[0][-1]]

        return alist

######## DamModel

    def GrowDam(self,doid,DamEl,parameters,N,scale,integration,R=None):
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

            # Bug fix, 2026 (see docs/PORTING_NOTES.md): this used to
            # always attempt RK5 first and only shrink the step
            # *reactively*, after an exception, by retrying with a step
            # divided by a fixed factor some number of times. That retry
            # count is sensitive to numerical happenstance unrelated to
            # the real (smoothly-varying) stress field, and was found --
            # visually, then quantitatively (adjacent hot-spot nodes
            # differing by thousands of cycles) -- to inject severe,
            # non-physical node-to-node noise.
            #
            # A first attempt fixed this by capping the step based on
            # proximity to Kic (Kmax/Kic) instead. That didn't hold up
            # either: RK5's own adaptive error control, even with the
            # Cash-Karp coefficients corrected, let the crack *double in
            # size in a single cycle* while still classified as "far from
            # instability" by that metric -- and right at the boundary,
            # the capping fraction came out to ~1.0, essentially no cap
            # at all exactly where it was needed most.
            #
            # Direct, unconditional fix instead: bound da/dN itself, as a
            # fraction of the current crack size, regardless of proximity
            # to Kic -- da/dN <= max_growth_fraction * a / dN, i.e. never
            # let a single step grow the crack by more than
            # max_growth_fraction (user-configurable via the .par file,
            # default 10%). This targets the actual failure mode directly
            # (a step too large relative to the current, possibly very
            # steep, local slope) instead of trying to predict it
            # indirectly via a proxy that turned out to be unreliable.
            # RK5 only ever shrinks a step *within* a single call (via its
            # own recursive refinement) and never grows it beyond what it
            # was given to attempt -- growth only affects the *next*
            # call's starting guess -- so capping nextdN here, every call,
            # against the current local slope is sufficient to bound the
            # actual step taken.
            try:
                rate_now = DamEl.dAdN(N,abab,[scale])
                positive = [rate_now[i] for i in range(len(abab)) if rate_now[i] > 0.0]
                if positive:
                    safe_dN = min(parameters.max_growth_fraction*abab[i]/rate_now[i] \
                                  for i in range(len(abab)) if rate_now[i] > 0.0)
                else:
                    safe_dN = self.DamOro[doid].nextdN
            except (DamErrors.FitPolyError, DamErrors.dAdNError, ValueError):
                rate_now = None
                safe_dN = self.DamOro[doid].nextdN

            if rate_now is not None and safe_dN < self.DamOro[doid].nextdN:
                dN = max(safe_dN, parameters.dN*_MIN_STEP_FRACTION)
                self.DamOro[doid].nextdN = dN
                rate = rate_now

            else:
                try:
                    rate,dN,self.DamOro[doid].nextdN = \
                        Integration.RK5_CKslope_vector(self.DamOro[doid].nextdN,\
                        N,abab,DamEl.dAdN,[scale],max_er,.05)

                except (DamErrors.FitPolyError, DamErrors.dAdNError, \
                        ValueError) as message:
                    if self.verbose:
                        print((' switch to one step fwd Euler for', \
                              'at doid', doid, 'due to', message))
                    dN = max(self.DamOro[doid].nextdN/10.0, parameters.dN*_MIN_STEP_FRACTION)
                    self.DamOro[doid].nextdN = dN
                    rate = Integration.Eulerslope_vector( \
                        dN,N,abab,DamEl.dAdN,[scale])

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
        DamEl.state_N.append(N+dN)  # see state_N's docstring (DamClass.py)

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

        As,types,alist,dNs=self.InterpolateLifeLists(doid)
        

        # if self.DamOro[doid].WillGrow = 0, means the LARGEST crack stopped 
        # growing because DKi < DKth, so return N = -1 (similar to setting
        # self.DamOro[#].Life = -1), or max life for doid,
        # self.DamOro[doid].Life, or perhaps even self.HighLife... ?????
        if self.DamOro[doid].WillGrow == 0 or \
           self.DamOro[doid].WillGrow == -1: return N_max*1.01

        # to determine the appropriate damage type and calc crack area:
        j = -1
        for i in range(len(alist)):
            if alist[i] > ai:
                j = i
                break

        try: 
            if j == 0: type = types[0]
            else: type = types[j-1]
        except IndexError:
            print((j, doid, ai, self.DamOro[doid].WillGrow))
            print((self.DamOro[doid].DamEl[-1]))
            print(alist)
            raise

        if type == 0:
            Ai = math.pi*ai*ai
        elif type == 1:
            Ai = 0.5*math.pi*ai*ai
        else: Ai = 0.25*math.pi*ai*ai

        # find j again by comparing crack areas 
        j = -1
        for i in range(len(alist)):
            if As[i] > Ai:
                j = i
                break

        # calculate life
        if j > 0:
            N=np.sum(dNs[j+1:])
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
        As=[]; types=[]; alist=[]; dNs=[]
        for DamEl in self.DamOro[doid].DamEl:
            DamEl.AppendInterpolateLifeLists(As,types,alist,dNs)
        return As,types,alist,dNs

######## DamModel

    def NFile(self,doid,nfile):
        '''
        writes info to nfile
        '''
        stuff=[]
        if self.parameters.monte:
            for set in self.DamOro[doid].N.rv:
                for key in list(set.keys()): # key is RID
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
            for i in list(self.DamOro.keys()):
                xyz,delxyz,sigxyz=self.model.GetNodeInfo(i)
                print('')
                print(('---------------- Damage Origin id = %8i' % (i), \
                      '------------------'))
                print((' coords         = ', xyz))
                print((' initial a      = ', self.DamOro[i].ai))
                print((' Predicted Life = ', self.DamOro[i].Life))

                if dam == 'all':
                    print(' **** Damage Elements: ')
                    for damel in self.DamOro[i].DamEl:
                        print(damel)

                if monte:
                    print(' ')
                    # initial area = initial area of full ellipse! 
                    print(('     set |    initial a    |  initial area   |',\
                    '  predicted life'))
                    print(dashplus)
                    if set == 'all':
                        for ij in range(self.parameters.sets):
                            init_a_set=list(self.ais.rv[ij].keys())
                            init_a_set.sort()
                            for ai in init_a_set:
                                RID=self.ais.rv[ij][ai]
                                for rid in RID:
                                    print(('%8i | %15.6e | %15.6e |   %.0f' % \
                                          (ij, \
                                          ai,(ai**2)*math.pi, \
                                          self.DamOro[i].N.rv[ij][rid])))
                            s="  Sample Mean = %10i " % \
                               (int(self.DamOro[i].N.SampleMean(ij)))
                            s+="|  Sample Stnd Dev = %10i" % \
                                  (int(math.sqrt(\
                                      self.DamOro[i].N.SampleVariance(ij))))
                            print(s)
                            print(dashplus)
                    else:
                        init_a_set=list(self.ais.rv[set].keys())
                        init_a_set.sort()
                        for ai in init_a_set:
                            RID=self.ais.rv[set][ai]
                            print(('%8i | %15.6e | %15.6e |   %.0f' % (ij, \
                                  ai,(ai**2)*math.pi, \
                                  self.DamOro[i].N.rv[set][RID])))
                        print(("         Sample Mean = ", \
                              int(self.DamOro[i].N.SampleMean(set))))
                        print(dashplus)
                print((dashes,"\n"))

        elif doid != 'none':
            print('')
            print(('---------------- Damage Origin id = %8i' % (doid), \
                  '------------------'))
            print((' coords         = ', self.__DamOro[doid].coords))
            print((' dam_list       = ', self.__DamOro[doid].dam_list))
            print((' initial a      = ', self.__DamOro[doid].ai))
            print((' Predicted Life = ', self.__DamOro[doid].Life))
            print((dashes,"\n"))

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
                                    self.IntegrationScheme,R)
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
                # dam_elem is a brand-new object (its __init__ set
                # state_N=[0.0]) -- it's the SAME physical state as
                # CurrentElement's current one, just re-expressed in the new
                # element's own geometry, so its absolute N is N, not 0.0.
                dam_elem.state_N = [N]
                self.DamOro[doid].DamEl+=[dam_elem]
                CurrentElement=dam_elem

        #end while

        # if will_grow = 0 or 2, compute life at DamOro and
        # set self.DamOro[doid].Life and update HighLife and LowLife
        if N > 0.0:

            if N > parameters.N_max:
                N = parameters.N_max*1.01

            # .Life is unconditional (this doid's own computed N, always) --
            # previously it was only set inside the three ifs below, so a
            # doid whose N was neither capped at N_max nor a new model-wide
            # running max/min never got .Life set at all (stayed at the
            # container's -1 "never computed" default, later silently
            # reported as self.HighLife via LifeValues's fallback, or
            # literally -1 via NFile's non-monte branch). Confirmed present,
            # byte-for-byte, in the original 2007 source (commit bd2d7d0) --
            # not introduced by the port. Only affects deterministic
            # (non-Monte-Carlo) runs: Monte Carlo life comes from a separate
            # object (DamOro[doid].N, via SampleMean) and never touches
            # .Life. Also a correctness prerequisite for multiprocessing:
            # pre-fix, .Life depended on what order *other* doids were
            # processed in (whether this one beat the running record at the
            # time), which doesn't even make sense once doids are split
            # across worker processes with their own independent
            # HighLife/LowLife. See docs/PORTING_NOTES.md.
            self.DamOro[doid].Life=N

            if N > self.HighLife:
                self.HighLife = N

            if N < self.LowLife and will_grow != 0:
                self.LowLife = N

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

    def RefreshLifeBounds(self):
        '''
        Recompute self.HighLife/self.LowLife from self.DamOro. Needed after
        a ddsim.parallel multiprocessing merge, where every doid was
        actually processed by some worker's own DamModel instance with its
        own process-local running trackers that started fresh (-1.0/1e100)
        -- the parent's own HighLife/LowLife never saw any of that, and
        wouldn't reflect the real merged result without this. A no-op-ish
        safety net when called after a normal serial run (every doid's
        SimDamGrowth call already updated these directly). Since the
        SimDamGrowth .Life fix (see docs/PORTING_NOTES.md), LifeValues's
        `self.HighLife if v==-1.0 else v` sentinel fallback is dead code in
        practice -- every processed doid's .Life is a real value -- but this
        keeps HighLife/LowLife themselves internally consistent for any
        other code that reads them directly.
        '''
        lives = [e.Life for e in self.DamOro.values() if e.Life != -1]
        if lives:
            self.HighLife = max(self.HighLife, max(lives))
            self.LowLife = min(self.LowLife, min(lives))

######## DamModel

    def LifeValues(self,set):
        '''
        Per-doid predicted life, shared by ToMAPFile and ToExodusFile: mean
        life across Monte Carlo samples (set), or the deterministic single-run
        Life, with self.HighLife as the sentinel for "no valid life computed".

        (Bug fix, 2026: ToMAPFile used to reference a bare 'HighLife' name here
        instead of 'self.HighLife' -- a latent NameError on the sentinel path,
        never exercised by the original test suite. Fixed by centralizing the
        lookup here.)
        '''
        values={}
        for i in list(self.DamOro.keys()): # i = doid
            if self.parameters.monte:
                v=self.DamOro[i].N.SampleMean(set)
            else:
                v=self.DamOro[i].Life
            values[i]=self.HighLife if v==-1.0 else v
        return values

    def ToMAPFile(self,MAPFile,set):
        '''
        Prints file for contouring in MAP

        MAPFile - file object to write to
        set  - the set from which to write mean value for
        '''

        values=self.LifeValues(set)
        DamKeys=list(values.keys())
        DamKeys.sort()

        MAPFile.write('LIFE 0'+"\n")
        for i in DamKeys: # i = doid
            MAPFile.write(str(i)+' '+str(values[i])+"\n")
        MAPFile.close()

######## DamModel

    def ToExodusFile(self,path,set=None):
        '''
        Write predicted life as an Exodus II nodal variable ("life"), viewable
        directly in ParaView (or any other Exodus-aware tool) as a contour
        plot -- reuses the exact same per-node values as ToMAPFile.

        path - output Exodus file path.  Always a NEW file; the mesh this
               model was built from (RDB or Exodus) is never modified.
        set  - the Monte Carlo set to use (ignored in deterministic mode)

        Nodes never run at all (not in self.DamOro) are written as NaN, which
        ParaView shows as masked/blank -- distinct from self.HighLife, which
        means a node WAS run and never reached failure (a real, meaningful
        long-life value, same sentinel ToMAPFile already uses).
        '''
        values=self.LifeValues(set)
        exodus_io.write_exodus(path,self.model.to_mesh_data(),{'life':values}, \
                               default=float('nan'))

######## DamModel

    def NToPickle(self,filename,doid):
        '''
        Pickle the life predicitons.  First build a dictionary of the N's,
        then pickle.
        '''
        if doid=='all':
            dumpfile=open(filename,'w+b')
            damage_list={}
            for i in list(self.DamOro.keys()): # doid 
                damage_list[i] = self.DamOro[i].N
            pickle.dump(damage_list,dumpfile,2)
            dumpfile.close()

        else:
            try:
                readfile=open(filename,'r+b')
                damage_list=pickle.load(readfile)
                readfile.close()
                damage_list[doid]=self.DamOro[doid].N
                dumpfile=open(filename,'w+b')
                pickle.dump(damage_list,dumpfile,2)
                dumpfile.close()
            except IOError:
                damage_list={doid:self.DamOro[doid].N}
                dumpfile=open(filename,'w+b')
                pickle.dump(damage_list,dumpfile,2)
                dumpfile.close()

######## DamModel

    def __init__(self,model,node_list,verbose,verify,extension,saveall, \
                 parameters,DebugGeomUtils=False,SVIEW=False, \
                 integration='RK5'):
        self.DamOro = {} # keys = doid, values = _DamOroContainer instance
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
                    if seid in self.SurfaceElements:
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
                    print(snid)
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
        keys=list(self.DamOro.keys())
        keys.sort()

        for Damage in self.DamOro[keys[0]].DamEl:
            for i in range(len(Damage.a[0][0])):
                qpnts=Damage.ComputeFrontPoints(i)
                for qpnt in qpnts:
                    frontfile.write(str(qpnt.x())+" "+str(qpnt.y())+" "+ \
                                    str(qpnt.z())+"\n")
                frontfile.write('################# \n')

######## DamModel

    def WriteCrackPathVTK(self,path,doid=None):
        '''
        Write one doid's full crack-growth history as a legacy VTK PolyData
        file (ASCII) -- one polyline per recorded growth step, tracing the
        crack front's shape at that step in real (global) coordinates, with
        CELL_DATA giving each step's cumulative cycle count (N) and crack
        type (0=Fellipse/embedded, 1=Hellipse/surface, 2=Qellipse/corner).
        Load directly in ParaView alongside the mesh/stress Exodus file
        (ToExodusFile) and color by "N" to see the predicted crack path grow
        over the node's life.

        Added 2026 (see docs/PORTING_NOTES.md). Like FrontPoints, only
        meaningful for one doid at a time -- if doid isn't given, the first
        (sorted) key in DamOro is used, same convention as FrontPoints.

        Only works for a doid whose full step-by-step history is still
        live, i.e. a serial (no -j) run: -j workers return a stripped-down
        DamEl holding only the final state (see _StrippedDamEl /
        StripDamElForTransport), by which point the history this needs is
        already gone.
        '''
        if doid is None:
            doid = sorted(self.DamOro.keys())[0]
        dam_els = self.DamOro[doid].DamEl
        if dam_els and isinstance(dam_els[0], _StrippedDamEl):
            raise RuntimeError(
                "doid %d's crack-growth history isn't available to plot -- "
                "it went through a -j worker process, which only keeps the "
                "final state (see _StrippedDamEl). Re-run this doid without "
                "-j to get a full crack path." % doid)

        type_index = {'Fellipse': 0, 'Hellipse': 1, 'Qellipse': 2}
        points = []   # flattened (x,y,z), shared across every step/line
        lines = []    # one list of point indices per recorded step
        line_N = []   # cumulative cycle count, one entry per line
        line_type = []
        line_step = []

        step = 0
        for stage in dam_els:
            n_states = len(stage.a[0][0])
            for i in range(n_states):
                cum_N = stage.state_N[i]  # already absolute (see
                                           # state_N's own docstring)
                qpnts = stage.ComputeFrontPoints(i)
                start = len(points)
                points += [(q.x(), q.y(), q.z()) for q in qpnts]
                idx = list(range(start, start + len(qpnts)))
                if stage.names[0] == 'Fellipse':
                    # Fellipse (fully embedded) is a closed loop; its 20
                    # sample points don't quite reach back around (phi runs
                    # [-1.0, 0.9], not [-1.0, 1.0)) -- repeat the first index
                    # to close it visually. Hellipse/Qellipse are genuinely
                    # open arcs that terminate at the free surface.
                    idx.append(start)
                lines.append(idx)
                line_N.append(cum_N)
                line_type.append(type_index.get(stage.names[0], -1))
                line_step.append(step)
                step += 1

        with open(path, 'w') as f:
            f.write("# vtk DataFile Version 3.0\n")
            f.write("DDSim predicted crack path, doid %d\n" % doid)
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write("POINTS %d float\n" % len(points))
            for x,y,z in points:
                f.write("%.8e %.8e %.8e\n" % (x,y,z))

            total_ints = sum(len(l)+1 for l in lines)
            f.write("LINES %d %d\n" % (len(lines), total_ints))
            for l in lines:
                f.write(str(len(l)) + " " + " ".join(str(i) for i in l) + "\n")

            f.write("CELL_DATA %d\n" % len(lines))
            f.write("SCALARS N float 1\n")
            f.write("LOOKUP_TABLE default\n")
            for n in line_N:
                f.write("%.8e\n" % n)
            f.write("SCALARS crack_type int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for t in line_type:
                f.write("%d\n" % t)
            f.write("SCALARS step int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for s in line_step:
                f.write("%d\n" % s)


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
        return a np.array containing the indices of a larger than the
        corresponding entry in aa.
        '''

##        z=[np.greater(a[k],aa[k]) for k in range(len(aa))]
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
                RANGE=list(range(lena))
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
                    DamEl.state_N.append(N)  # see state_N's docstring (DamClass.py)

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
                # damel is a brand-new object (its __init__ set
                # state_N=[0.0]) -- it's the SAME physical state as DamEl's
                # current one, just re-expressed in the new element's own
                # geometry, so its absolute N is N, not 0.0.
                damel.state_N = [N]
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

            DamKeys=list(self.DamOro.keys())
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

            DamKeys=list(self.DamOro.keys())
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

            DamKeys=list(self.DamOro.keys())
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









