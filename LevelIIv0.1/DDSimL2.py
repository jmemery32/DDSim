import sys, VarAmplitude, Parameters, fx, Integration
import dadN

class DamEl:
    '''
    class to store K vs. a history and return K's
    '''

    def __init__(self,dadN,SIFFile,ai,SIFFilescale):
        self.dadN=dadN
##
##        # initial crack size, SIFFile is crack lenght - ai
##        self.ai = ai

        # store a and K pairs when ever CalculateKi is called.  Here are the
        # empty lists
        self.a = [ai]
        self.K = []

        # store dn's during simulations, here is the empty list. 
        self.dN = []

        # to scale SIF up
        self.SIFFilescale=SIFFilescale

        # initiate the method to calc K... 
        buff = SIFFile.readline()
        a=[]; K=[]
        while len(buff) > 0:
            line = buff.split()
            a+=[float(line[0])+self.a[0]]
            K+=[float(line[1])*SIFFilescale]
            buff = SIFFile.readline()

        self.amax = max(a)

        self.KCalculator = fx.OneVariable(a,K)

########### DamEl ###############

    def CalculateKi(self,scale=None,a=None):
        '''
        scale is the %DK for variable amplitude, 1.0 for constant amp.  
        '''

        if not scale: scale = 1.0

        if not a:
            a=self.a[-1]
            K = scale*self.KCalculator.Evaluate(a)
            self.K += [K]

        else: K = scale*self.KCalculator.Evaluate(a)

        return K

########### DamEl ###############

    def WillGrow(self,N,R,Ki): 
        '''
        Given the appropriate damage parameter,
        in this case stress intensity factor, and sets will_grow = :

        will_grow = -1 --> stress field is compressive (ie. Ki < 0.0)
        will_grow = 0  --> Ki insufficient to drive crack (Ki < K threshold)
        will_grow = 1  --> Stable, fatigue growth
        will_grow = 2  --> Unstable growth
        will_grow = 4  --> Assume net fracture (crack grew outside body)

        '''

        if Ki==4: return 4

        Kic=self.dadN.Kc
        Reff=self.dadN.CalcReff(int(N),R)
        Kth=self.dadN.Calc_DKth(self.a[-1],Reff)
        Kmax=Ki/(1.0-R)

        # max K is negative = compression, will not grow.
        if Kmax <= 0.: return -1

        # will not grow, even under cyclic loading if...
        elif Ki <= Kth: return 0

        # as in unstable growth if...
        elif Kmax >= Kic: return 2

        else: return 1

########### DamEl ###############

    def dAdN(self,N,size,args):
        '''
        Used by GrowDam 
        N - current number of cycles of loading
        size - a list of ellipses dimensions: [a, b, -a, -b]
        args - tuple of supporting arguments, in this case just scale! 
        '''
        # unpack args
        scale=args[0]
        abab= max(0,size[i])
        KIs = self.CalculateKi(abab)
        growth = self.dadN.Calc_dadN(KIs,size[0],int(N))

        return [growth]

#################### End of DamEl class ################

def VariableAmplitude(DamEl,Nmax,spec,nore,r):
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
    # -1 - compressive stress field:  *** WON'T ENTER VarAmp()
    # 0 - subcritical: if through load program with no growth
    # 1 - will grow: while crack will still  grow
    # 2 - unstable: when DK > DKc
    # 3 - change shape: NEVER within VarAmp()
    # 4 - outgrown: net fracture
    WillGrow = 1

    # initialize some stuff
    spec.count=0

    # do the cycle-by-cycle integration (vectorized)
    twice=0; count=0; # used looping through the sprectrum w/o growth
    N=0.0; Nold=0.0; inc=0; Ks=0
    aCurrent = DamEl.a[-1]

    while WillGrow == 1:
    # first check the current damage elements geometry... 
        if inc:
            # this is a little sloppy, RANGE is defined below... 
            test = (aCurrent*r > aa)
##            if test: pass
##            else: damel,doubt = DamEl.GeometryCheck(self.CMesh)
##        else: damel,doubt = DamEl.GeometryCheck(self.CMesh)

        aa = DamEl.a[-1]
        percentDK, R = spec.Delta()
        if inc:
            if test:
                Ks = KsCurrent*percentDK
                wg = DamEl.WillGrow(N,R,Ks)

            else:
                aCurrent=aa
                Ks = DamEl.CalculateKi(percentDK)
                wg = DamEl.WillGrow(N,R,Ks)
                if wg != 4:
                    KsCurrent = Ks/percentDK

        else:
            aCurrent = aa
            Ks = DamEl.CalculateKi(percentDK)
            wg = DamEl.WillGrow(N,R,Ks)
            if wg != 4:
                KsCurrent = Ks/percentDK

        DamEl.dN+=[N-Nold]
        Nold=N

        # unstable crack growth
        if wg == 2:
            N+=1.0
            WillGrow=2

        # net fracture
        elif wg == 4:
            N+=1.0
            WillGrow=4

        # if it is a compressive cycle, increase N and continue
        # still don't have amplification for under load... 5/27/06
        elif wg == -1:
            if twice == 0:
                twice = 1
                count +=1
                N+=1.0
            elif count < spec.length:
                count+=1
                N+=1.0
            else:
                N+=1.0
                WillGrow=-1

        else:
            # compute the increment of growth
            if nore: inc = DamEl.dadN.Calc_dadN(Ks,aa,int(N),R)
            else: inc = DamEl.dadN.Compute_dadN(Ks,aa,int(N),R)

            aa += inc

            # Update the damage element
            DamEl.a += [aa]

            # Check if aa exceeds aR
            if aa >= DamEl.amax: break 

            # if N exceeds max N return N = -1 for processing at
            # higher level
            elif N >= Nmax:
                N = -1
                WillGrow = 0

            # to prevent repetitive times through the load program
            # with inc = 0.0...
            elif inc==0.0 and twice == 0:
                twice = 1
                count +=1
                N+=1.0
            elif inc==0.0:
                if count < spec.length:
                    count+=1
                    N+=1.0
                else:
                    N = -1
                    WillGrow = 0
            else:
                twice = 0
                count = 0
                N+=1.0
    # end while

    return N, aa, WillGrow

#################### 

def ConstantAmplitude(CurrentElement,parameters,scale,R,IntType):
    '''
    Simulates the growth of damage until it is not interesting or it
    reaches some critical, life limiting condition.  Results in the life
    prediction resulting from damage at the given damage origin.

    *** Integration in "time" ***,
    '''

    # collect the corner and midside nodes of the surface mesh for use with
    # WillGrow in finding intersection of ellipses with surface facets.

    Flag=0; N=0.

    while Flag == 0:

        Ki = CurrentElement.CalculateKi(scale)
        will_grow = CurrentElement.WillGrow(N,R,Ki)
        # will_grow = -1 if Ki < 0.0, then stop growing crack
        if will_grow == -1: ## and monte == 0: #(see ** below.)
            Flag=1
            WillGrow = -1
            dN=0.0

        # will_grow = 0 if Ki is less than Kth then stop growing crack
        elif will_grow == 0:
            Flag=1
            WillGrow = 0
            dN=0.0

        # Parameter will_grow is set to 1 in WillGrow if the crack,
        # indeed, will grow stably.  (ie. Kth<Ki<Kic)
        elif will_grow == 1:
            dN = GrowDam(CurrentElement,parameters,N,scale,IntType)
            N += dN

            # if N has reached a practical maximum... 
            if N >= parameters.N_max:
                Flag=1
                WillGrow = 1

        # will_grow = 2 means unstable growth (ie. Ki>=Kic)
        elif will_grow == 2:
            Flag=1
            WillGrow = 2
            dN=0.0

        # will_grow = 4 crack grew outside body, assume net fracture
        elif will_grow == 4:
            Flag=1
            WillGrow = 4
            dN=0.0

        # Update dN
        CurrentElement.dN+=[dN]
    #end while

    return N, DamEl.a[-1], WillGrow

#################### 

def GrowDam(DamEl,parameters,N,scale,IntegrationScheme):
    '''
    scheme to calculate the grow increment for constant amplitude loading.
    '''

    # for RK5 integration
    nextdN = 0.0

    assert (len(DamEl.a)>0) # i.e there is a SIF History... 

    # Get current crack length...
    abab = DamEl.a[-1]

    # change to percentage of minimum dimension
    max_er=parameters.max_error*abab
    max_er=min(max_er,parameters.max_error)

    # time step:
    dN=parameters.dN

    if integration=='FWD':
        # use NASGRO eqn.
        # calc.dadN...
        rate=DamEl.dadN.Calc_dadN(DamEl.a[-1],DamEl.K[-1],int(N)) 

    elif integration=='RK4':
        # assume always using NASGRO equations with RK4 scheme...
        rate = Integration.RK4slope_vector(dN,N,abab, \
                                           DamEl.dAdN,[scale])

    else: # if self.IntegrationScheme=='RK5'
        # assume always using NASGRO equations with RK5 scheme...
        if (nextdN == 0.0): nextdN = dN

        # the integrations scheme can sometimes become unstable and pass
        # calcki numbers that don't make sense (too big).  Hence, we put in
        # a try and except block.  if RK5 scheme is bunk we do on small
        # step by the fwd euler method in hopes that springs us loose of
        # bad spot.  may need to add a finally or something to this in the
        # future.

        rate,dN,nextdN = \
            Integration.RK5_CKslope_vector(nextdN,\
            N,abab,DamEl.dAdN,[scale],max_er,.05)


    # Caclulate the new a, b, a-, b-.
    # a, b and neg_a, neg_b never to take negative value (which shouldn't
    # happen) so abs() is used to insure that here.
    inc=abs(rate)*dN

    if integration=='RK5': # don't mess-up the adaptivity
        abab_new=abab+inc
    else: 
        # don't grow less than min_inc...
        if inc < parameters.min_inc:
            dN=parameters.min_inc/rate
            inc=abs(rate)*dN
            abab_new=abab+inc

        # otherwise, just do what seems natural... 
        else:
            abab_new = abab+inc

    # update Damage element
    DamEl.a+=[abab_new]

    return dN

#################### 

def DumpFile(xlist,ylist,filename,freq):
    '''
    function to dump items in list out to a text named filename.
    '''

    if len(xlist) == len(ylist)+1: xlist.pop(-1)
    elif len(ylist) == len(xlist)+1: ylist.pop(-1)

    s=''
    for i in range(len(xlist)):
        if i%freq == 0: s+=(str(xlist[i])+' '+str(ylist[i])+"\n")
        elif i == len(xlist)-1: s+=(str(xlist[i])+' '+str(ylist[i])+"\n")

    fil=open(filename,'a')
    fil.write(s) 
    fil.write("# *** \n")
    fil.close()

#################### 

def TakeCommandLineArgs():
    '''
    takes the command line flags.
    '''

    surface = 0 # 0 = interior, 1 surface point for NASGRO eqn
    if "-surface" in sys.argv: surface = 1
        
    filename = None
    if "-base" in sys.argv:
        index = sys.argv.index('-base')+1
        filename = sys.argv[index]

    if "-siffile" in sys.argv:
        index = sys.argv.index('-siffile')+1
        siffilename = sys.argv[index]

    ai = 0.015 # default initial crack size
    if "-ai" in sys.argv:
        index = sys.argv.index('-ai')+1
        ai = float(sys.argv[index])

    VarAmp = False
    if "-VarAmp" in sys.argv: VarAmp = True

    # when to recompute K in the variable amplitude loading scenario
    GrowthRatio = 0.
    if "-GrowthRatio" in sys.argv:
        index = sys.argv.index('-GrowthRatio')+1
        GrowthRatio = float(sys.argv[index])

    # scale the applied scale for the SIFs.  K's
    # should for the maximum load in the spectrum (or cylce).
    SIFFilescale = 1.0
    if "-SIFFilescale" in sys.argv:
        index = sys.argv.index('-SIFFilescale')+1
        SIFFilescale = float(sys.argv[index])

    nore = 0 # to turn off crack retardation effects
    if "-nore" in sys.argv: nore = 1

    # specify the type of integration for a costant amplitude simulation
    # default --> adaptive, 5 pnt. runga-kutta scheme
    # -Fwd --> forward euler
    # -RK4 --> 4 point RK scheme
    # -Simp --> simpson's rule on generated SIF history
    Int_type = 'RK5'
    if "-Fwd" in sys.argv: Int_type = 'FWD'
    elif "-RK4" in sys.argv: Int_type = 'RK4'
    elif "-Simp" in sys.argv: Int_type = 'Simp'

    Dump = None
    if "-Dump" in sys.argv:
        index = sys.argv.index('-Dump')+1
        Dump = sys.argv[index]

    freq = 1
    if "-freq" in sys.argv:
        index = sys.argv.index('-freq')+1
        freq = int(sys.argv[index])

    return surface,filename,ai,VarAmp,SIFFilescale,nore,Int_type, \
           GrowthRatio,Dump,freq,siffilename

def PrintHelps():
    print " "
    print " -surface to direct NASGRO equation to use plain strain or plane"
    print "         value of alpha"
    print " -base [string] to indicate the filename that stores the parameters "
    print "       and variable amplitude spectrum "
    print " -siffile [string.string] to indicate the filename that stores the sif history "
    print " -ai [float] to inidicate initial crack size (default 0.015)"
    print " -VarAmp to indicate variable amplitude loading"
    print " -GrowthRatio [float] maximum size to grow crack between computing K"
    print "              default is 0."
    print " -SIFFilescale [float] to scale SIF up" 
    print " -nore to suppress the crack retardation model"
    print " "
    print " * Constant amplitude integration schemes (adaptive RK5 is default):"
    print " -Fwd to force forward integration"
    print " -RK4 to force RK4 scheme."
    print " "
    print " -Dump [string] file name to write a vs. K data to."
    print " -freq [int] when N is even divisible by -freq dump string"
    print "             (default = 1)"
    sys.exit()

#################### 

if __name__=="__main__":

    if "-h" in sys.argv: PrintHelps()

    # command line arguments
    surface,filename,ai,VarAmp,SIFFilescale,nore,IntType,GrowthRatio,Dump, \
    freq,siffilename = TakeCommandLineArgs()

    # initialize some stuff
    Spec = VarAmplitude.Spectrum(filename+'.val')
    parameters = Parameters.Parameters(filename,'')
    R = parameters.material[16]
    dadN = dadN.dadN(parameters.material,surface)

    SIFFile = open(siffilename)
    DamEl = DamEl(dadN,SIFFile,ai,SIFFilescale)
    SIFFile.close()

    # do the life prediction
    if VarAmp: N,af,WillGrow = VariableAmplitude(DamEl,parameters.N_max,Spec, \
                                                 nore,GrowthRatio)
    else: N,af,WillGrow = ConstantAmplitude(DamEl,parameters,scale,R,IntType)

    print "\n Life prediction is: "
    print " N = %i, af = %f, WillGrow = %i, Kmax = %f" % (N,af,WillGrow,DamEl.K[-1])

    if Dump: DumpFile(DamEl.a,DamEl.K,Dump,freq)



















