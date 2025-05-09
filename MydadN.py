import math
import Numeric

import dadN as PYDdadN

import Newton # for Willenborg()

### Exceptions
##NewtonSolve_ExceedMaxIts = "Newton.Solve() no convergence"

class dadN:

    ''' This class to calculate da/dN by the NASGRO equation.
    ref. NASGRO Reference Manual, version 4.02, Sept. '02
    '''

    def __init__(self,materials,surf):
        ''' Inputs are:
        materials - a list of material properties
        '''

        #  materials**    = [UTS, YS, Kie, Kic, Ak, Bk, C, n, p, q, DKo, \
        #                    Cthp, Cthm, Rcl, alpha, Smax/SIGo, R, a0, np]
        #
        #  **notes on that:
        #  1. These parameters can be found in NASGRO and the comments here are
        #     pretty much from the NASGRO manual.
        #  2. It is conservative to use Cthp, Cthm as zero.
        #  3. Kie is used for Kc (in NASGRO dadN equation) for part-through
        #     cracks.
        #  4. For through cracks the NASGRO manual has an equation for Kc which
        #     is a function of Kic (eqn 2.14).
        #  5. For Fellipse class of cracks Kc = Kic.

        # materials, as per NASGRO
        self.UTS   = float(materials[0]) # mat'l ultimate strength
        self.YS    = float(materials[1]) # mat'l yield strength
        self.Kie   = float(materials[2]) # part-through fracture toughness
        self.Kic   = float(materials[3]) # plane strain fracture toughness
        self.Ak    = float(materials[4]) 
        self.Bk    = float(materials[5])
        self.C     = float(materials[6])
        self.n     = float(materials[7])
        self.p     = float(materials[8])
        self.q     = float(materials[9])
        self.DK1   = float(materials[10]) # modified 3/4/06 jme32 per NASGRO4.2
        self.Cthp  = float(materials[11])
        self.Cthm  = float(materials[12])
        self.Rcl   = float(materials[13]) # what the hell is this?(jme 1/17/05)
        self.alfa  = float(materials[14]) # plane strn/strss constraint factor
        self.ratio = float(materials[15]) # Smax/SIGo
        self.R     = float(materials[16]) # Smax/Smin
        self.a0    = float(materials[17]) # small crack param (NASGRO eqn 2.12)
        self.np    = float(materials[18])
        self.damp  = float(materials[19]) # to reduce N in small crack
                                          # modification to NASGRO eqn., see \
                                          # CalcReff...
        try: self.Rso   = float(materials[20]) # stress ratio shutoff (see \
                                               # Willenborg.Willenborg())
        except IndexError: self.Rso = 3.0 # recommended value in NASGRO manual

        self.surf  = surf # 1 if surface crack, 0 if NOT

        # Newmans surface correction for plastic zone.
        #       See Willenborg.Willenborg()
        if surf: 
            self.alphag = 1.15
        else: self.alphag = 2.55 

##        self.monte = monte # monte = 1.0 if monte carlo simulation, 0.0 o.w.

        # calculate Ai's, they don't change with different a or DKi or R...
        self.A0 = (0.825-0.34*self.alfa+0.05*self.alfa*self.alfa)* \
                  ((math.cos((math.pi/2.0)*self.ratio))**(1/self.alfa))
        self.A1 = (0.415-0.071*self.alfa)*self.ratio
        self.A3 = 2.0*self.A0+self.A1-1.0
        self.A2 = 1.0-self.A0-self.A1-self.A3

        if self.surf==1:
            self.Kc=self.Kie
        else:
            self.Kc=self.Kic

######### dadN

    def Calc_dadN(self,DKi,a,N,R=None,DKth=None):
        ''' this method caclulates dadN per the NASGRO equation.
        Arguments:

        DKi -  delta K_I
        a -  crack length
        N -  current number of loading cycles.  To account for small/new cracks
        that have no or little plastic wake induced closure effects.
        R - R ratio (if not passed in, constant amplitude and R is located in
            marterials, so dadN gets from self.R) 
        DKth - Can pass in or calculate here.  
        '''

        # .CalcReff() is my effective R to account for small/new crack growth.
        # This is NOT NEWMAN'S CRACK CLOSURE effective stuff.  It
        # assumes that with lower N, the crack will grow, regardless of if it
        # is below threshold, because there is not plasticity induce closure.

        if R != None and R < 100.:
            Reff = self.CalcReff(N,R)
            # compute Kmax here using original R, NOT Reffective.  We want
            # Kmax to be the same.  Must then recompute DKi using Reff.
            Kmax = DKi/(1.0-R)
##            print Reff, Kmax
##            DKi = Kmax*(1.0-Reff)
        else:
            Reff = self.CalcReff(N,9999)
            # compute Kmax here using original R, NOT Reffective.  We want
            # Kmax to be the same.  Must then recompute DKi using Reff.
            Kmax = DKi/(1.0-self.R)
##            DKi = Kmax*(1.0-Reff)

        if not DKth:
            DKth=self.Calc_DKth(a,Reff)

        # add n' branch to the R = 0.7 curve
        DKn = 0. # initialize it so DKi is always greater

        # note to self:  if want to plot a pure R = 0.7 curve can force by
        # adding:  [and N < 100.0] here and below in next logic statement... 
        if abs(0.7 - Reff) <= 0.05: # and N < 100.0:
            R=math.exp(0.5*(math.log(DKth)+math.log(self.Kc)))
            DKn = self.__CalcDKn(DKth+0.001,R,Reff,DKth,Kmax)
            fn = self.__NASGRO(DKn,Reff,DKth,Kmax)
            Cn = math.exp(math.log(fn)-self.np*(math.log(DKn)))

        # calculate dadN per NASGRO eqn. (2.1) with modification to R = 0.65 -
        # 0.7 to make it Paris type (linear in log-log space) with slope n'
        if DKi <= DKn and abs(0.7 - Reff) <= 0.05: # and N < 100.0:
            dadN = Cn*((DKi)**self.np)

        else:
            # if DKi < DKth the applied delta SIF is less than the
            # threshold value so the growth rate should be zero.
            # Seems trivial, like why call this routine if
            # DKi < DKth, but can happen if one point on crack front
            # wants to grow, while the others do not.
            if DKi < DKth:
                return 0.0 
            # calculate Kmax
            dadN = self.__NASGRO(DKi,Reff,DKth,Kmax)

        return dadN

######### dadN

    def __NASGRO(self,DKi,Reff,DKth,Kmax):
        '''
        compute dadN per NASGRO eqn. (2.1) straight up given all parameters.
        '''

        # get f
        f = self.Setf(Reff)

        DKeff=((1.0-f)/(1.0-Reff))*DKi

        left=self.C*((DKeff)**self.n)
        ned=(1.0-DKth/DKi)**self.p
        donkey=(1.0-Kmax/self.Kc)**self.q

        return left*ned/donkey

######### dadN

    def CalcKeff(self,DKi,N,Reff=None):
        '''
        To calculate Delta(K) effective.  This is calculated above for da/dN
        but i need this method so i can vectorize it for plots.

        *** Basically, this is totally a gravy function and isn't called from
        wihtin DDSim. ***
        '''

        if not Reff: Reff = self.CalcReff(N) # R effective

        f = self.Setf(Reff)

        return ((1.0-f)/(1.0-Reff))*DKi

######### dadN

    def Calc_DKth(self,a,R):
        '''
        method to caclulate the threshold delta(K) value.  R is passed in, in
        lieu of self.R because it should be an effective R calculated by
        self.CalcReff().

        a - crack size
        R - R effective
        '''

        # calculate f first. 
        f = self.Setf(R)

        # calculate DK1
##        DK1=self.DKo*((1.0-self.A0)**(1.0+self.Cthp)) # NASGRO eqn. (2.13)

        # we need DK1*.  Since i'm starting with a very small initial flaw
        # i.e. O(intrinsic flaw size of mat'l) a could be smaller than a0.
        # In which case, use a/2*a instead of a/(a+a0) (i.e. say a0 = a).
        if a > self.a0: 
            DK1_star=self.DK1*((a/(a+self.a0))**0.5) # NASGRO eqn. (2.12)
        else:
            DK1_star=self.DK1*(math.sqrt(0.5)) # alternate NASGRO eqn. (2.12)

        # calculate the threshold SIF range, DKth
        if R >= 0.0: Cthm = self.Cthp
        else: Cthm = self.Cthm
        ned=((1.0-R)/(1.0-f*R))**(1.0+R*Cthm)
        donkey=(1-self.A0)**(self.Cthp-R*Cthm)

        DKth=DK1_star*ned/donkey # NASGRO eqn. (2.11)
        return DKth

######### dadN

    def Paris(self,Dki):
        '''
        compute da/dN from simple paris law.
        '''

        dadN=self.C*(Dki**self.n)

        return dadN

######### dadN

    def PrintMaterialModel(self):
        print ' UTS = ', self.UTS
        print ' YS = ', self.YS
        print ' Kie = ', self.Kie
        print ' Kic = ', self.Kic
        print ' Ak = ', self.Ak
        print ' Bk = ', self.Bk
        print ' C = ', self.C
        print ' n = ', self.n
        print " n' = ", self.np
        print ' p = ', self.p
        print ' q = ', self.q
        print ' DK1 = ', self.DK1
        print ' Cthp = ', self.Cthp
        print ' Cthm = ', self.Cthm
        print ' Rcl = ', self.Rcl
        print ' alpha = ', self.alfa
        print ' Smax/SIGo = ', self.ratio
        print ' R = ', self.R
        print ' a0 = ', self.a0
        print ' Is Surface = ', self.surf

######### dadN

    def CalcReff(self,N,R=None):
        '''
        Reff is R effective for small/new cracks in
        which there is no or little plastic wake induced closure effects.  The
        idea is that for high N, Reff = R.  Otherwise, Reff = hard (see below).
        Therefore, for this to be valid R <= 0.7 to begin with.
        
        N - number of loading cycles the crack has currently been load for
        make damp << 1.0 for no Reffective (except for N = 0)
        '''

        hard = 0.7 # hard coded value of Reff when N is small

        if R > 100.:
            R = self.R

        A = (1.0 - math.exp(-N/self.damp))
        B = (R - hard)

        if R <= 0.7:
            Reff = A*B+hard
        else:
            Reff = R

        return Reff

######### dadN

    def Setf(self,Reff):
        '''
        return f  
        '''

        # calculate the crack opening function
        if Reff >= 0.0:
            f=max(Reff,(self.A0 + \
                             self.A1*Reff + \
                             self.A2*Reff*Reff + \
                             self.A3*Reff*Reff*Reff))
        else:
            f=self.A0+self.A1*Reff

        return f

######### dadN

    def SetKol(self):
        ''' dummy function so I can interchange with dadN.pyd'''
        self.Kol=0.0

######### dadN

    def __CalcDKn(self,L,R,Reff,DKth,Kmax):
        '''
        bisection routine to calculate DKn, the DK at which the n' line
        intersects the NASGRO equation.  
        '''
        def Fprime(DK,DKth,Kmax,Reff):
            '''
            private function to evaluate d/dDK(da/dN) at DK for use in the
            bisection.
            '''

            # get f
            f = self.Setf(Reff)

            # calculate Kmax
            DKeff=((1.0-f)/(1.0-Reff))*DK
            left=self.C*((DKeff)**self.n)
            ned=(1.0-DKth/DK)
            donkey=(1.0-Kmax/self.Kc)
            return (left*(ned**self.p)/(donkey**self.q))* \
                   (self.n/DK + \
                    self.p*DKth/DK/DK/ned + \
                    self.q/donkey/self.Kc)

        tol = L/1000.0

        while abs(R-L) > tol:
            m = L + (R - L)/2.0
            FL = (Fprime(L,DKth,Kmax,Reff) / self.__NASGRO(L,Reff,DKth,Kmax))*L
            Fm = (Fprime(m,DKth,Kmax,Reff) / self.__NASGRO(m,Reff,DKth,Kmax))*m
            if Numeric.sign(FL-self.np) == Numeric.sign(Fm-self.np):
                L = m
            else:
                R = m

        return L+(R-L)/2.0

####### end dadN

class Willenborg(dadN):
    '''
    Class to use for variable amplitude loading situations where retardation
    after overloads is important to include.

    ** There is no constructor because i want this class to inherit the
    constructor from it's parent class, dadN.
    '''

####### Willenborg

    def dDKthbydReff(self,a,Reff_in):
        '''
        method to calculate the d(DeltaK_th)/dReff.
        '''

        # set f
        f = self.Setf(Reff_in)

        # we need DK1*.  Since i'm starting with a very small initial flaw
        # i.e. O(intrinsic flaw size of mat'l) a could be smaller than a0.
        # In which case, use a/2*a instead of a/(a+a0) (i.e. say a0 = a).
##        DK1=self.DKo*((1.0-self.A0)**(1.0+self.Cthp)) # NASGRO eqn. (2.13)
        if a > self.a0: 
            DK1_star=self.DK1*((a/(a+self.a0))**0.5) # NASGRO eqn. (2.12)
        else:
            DK1_star=self.DK1*(math.sqrt(0.5)) # alternate NASGRO eqn. (2.12)

##        print "f,Reff_in",f,Reff_in

        if Reff_in >= 0.0:
            Cthm = self.Cthp
##            if Reff_in >= self.A0+self.A1*Reff_in+self.A2*Reff_in*Reff_in+ \
##                          self.A3*Reff_in*Reff_in*Reff_in:
            # if f = Reff_in
            if abs(Reff_in/f - 1.0) < 0.0001:
                Gp = -1.0/(1.0+Reff_in)
            # else if f = A1*R+A2*R*R+A3*R*R*R
            else:
                Gp = ((f*Reff_in-1.0)+(1.0-Reff_in)*(self.A0+ \
                      2.0*self.A1*Reff_in+3.0*self.A2*Reff_in*Reff_in+\
                      4.0*self.A3*Reff_in*Reff_in*Reff_in)) / \
                      ((1.0-Reff_in)*(1.0-f*Reff_in))
        else:
            Cthm = self.Cthm
            Gp = ((f*Reff_in-1.0)+(1.0-Reff_in)*(self.A0+ \
                      2.0*self.A1*Reff_in)) / \
                      ((1.0-Reff_in)*(1.0-f*Reff_in))

        
        B = (1.0-Reff_in)/(1.0-f*Reff_in)
        C = 1.0 - self.A0
        D = 1.0+Reff_in*Cthm 
        Dp = Cthm
        G = math.log((1.0-Reff_in)/(1.0-f*Reff_in))

        U = (B)**(D)
        Up = U*(Dp*G+D*Gp)
        V = (C)**(Reff_in*Cthm - self.Cthp)
        Vp = V*Cthm*math.log(C)

        deriv = DK1_star*(Up*V+U*Vp)

        return deriv

####### Willenborg

    def SetKol(self):
        '''
        quick method to reset Kol for each application of the load program.
        '''

        self.Kol = 0.0 # the overload SIF to start with.
        self.zo = 0.0 # the overload plastic zone size.
        self.aol = 0.0 # crack length at overload.

####### Willenborg

    def Willenborg(self,DKi,a,N,R=None):
        '''
        Method that returns da/dN for variable amplitude loading and positive
        R ratios.  Computes retardation using:

        gallagher, AFFDL-TM-FBR-74-28
        willenborg et al., AFFDL-TM-FBR-71-1
        NASGRO 4.0 reference manual (for funny phi thing)

        DKi - current delta K
        a - current crack length
        N - current number of cummulative cycles (used for new/small crack
            stuff)
        R - Kmin/Kmax
        alpha - the funny alpha thing in the plastic zone size estimate given 
                in the NASGRO manual and attributed to Newman in private comm.
        '''

        # if R is not passed in, use the one in the .par file that is for
        # constant amplitude.  **Note:  Willenborg's model is for retardation
        # due to overload cycles in variable amplitude loading, hence, it only
        # makes sense that R should be passed in!  
        if not R:
            R = self.R

        # compute incoming Kmax and Kmin
        Kmax = DKi/(1.0-R)
        Kmin = Kmax - DKi

        # if Kmax > current self.Kol (overload) store at Kol and return da/dN
        # with no retardation
        if Kmax >= self.Kol:
            self.Kol = Kmax
            self.zo = (math.pi/8)*((Kmax/(self.alphag*self.YS))**2.0)
            self.aol = a
            return self.Calc_dadN(DKi,a,N,R)

        # and if Kmax is less than zero, return zero.  
        elif Kmax < 1e-10: return 0.0

        # if not,  compute some parameter for more logic
        da=a - self.aol
        if da < 0.0: raise "dadN.Willenborg() delta_a is negative."
        elif da > self.zo:
            self.Kol = Kmax
            self.zo = (math.pi/8)*((Kmax/(self.alphag*self.YS))**2.0)
            self.aol = a
            return self.Calc_dadN(DKi,a,N,R)
        else: Kstar = self.Kol*(math.sqrt(1.0-(da/self.zo)))

##        self.PrintMaterialModel()
##        print "Kmin  , Kmax  , Kstar ,   Kol ,    da ,   zo , a"
##        print "%.4f, %.4f, %.4f, %.4f, %1.4e, %1.4e, %1.4e" \
##              %(Kmin, Kmax, Kstar,self.Kol,da,self.zo,a)
##        print ' '

        # if DKi exceeds Kstar continue with no retardation
        if Kmax >= Kstar:
            self.Kol = Kmax
            self.zo = (math.pi/8)*((Kmax/(self.alphag*self.YS))**2.0)
            self.aol = a
            return self.Calc_dadN(DKi,a,N,R)
        elif self.Kol/Kmax > self.Rso: return 0.0

        # for all other possible scenarios, calc Kr and Reff.  
        else:
            Krw = Kstar - Kmax # willenborg's Kreduction
            # extra stuff that needs to be passed along to Newton.Solve() used
            # in self.RandK_Func() & self.RandK_Jacob
            cdata = (a,Krw,DKi,Kmax,Kmin)
            # Initial guesses for Newton solve:
            # *** Note --> turns out 0., 0. is a more stable initial guess. 
            Kri = (1.0/(self.Rso-1.0))*Krw
            Reffi = (Kmin - Kri)/(Kmax - Kri)
            # Kr is phi*Krw
            Kr,Reff = Newton.Solve(2,(0.,0.),1e-10,self.RK_Func,\
                                     self.RK_Jacob,cdata,100)
##            print " ---> Reffi, Kri:", Reffi, Kri
##            print " ---> R, Krw:", R, Krw
##            print " ---> Reff, Kr:", Reff, Kr
##            print " ---> da/dN:", self.Calc_dadN(DKi,a,N,Reff,None)
            # passing None for DKth, might want to check if is compatible with
            # newton solve above.  
            return self.Calc_dadN(DKi,a,N,Reff,None)

####### Willenborg

    def RK_Func(self,cdata,del_X,X,RandK_value):
        '''
        Evaluate the functions:

        F1(Kr,Reff) = F(Reff) - Kr
        F2(Kr,Reff) = G(Kr) - Reff

        where: F(Reff) = phi*Krw
               G(Kr) = (Kmin - Kr)/(Kmax - Kr)
        '''

        # Sum X (initial guess) with total change in X (del_X)
        Kr_in = X[0] + del_X[0]
        Reff_in = X[1] + del_X[1]

        # constants per load cycle:
        a = cdata[0] # current crack size 
        Krw = cdata[1] # Willenborg's Kreduction 
        DKi = cdata[2] # current DeltaK
        Kmax = cdata[3] # current Max K
        Kmin = cdata[4] # current Min K

        # find solution to F1: 
        DKth = self.Calc_DKth(a,Reff_in)
        phi = (1.0 - DKth/DKi)/(self.Rso - 1.0)
        RandK_value[0] = phi*Krw - Kr_in

        # find solution to F2:
        RandK_value[1] = (Kmin - Kr_in)/(Kmax - Kr_in) - Reff_in

####### Willenborg

    def RK_Jacob(self,cdata,del_X,X,Jacob,GBL_count):
        '''
        Evaluate the Jacobian of derivatives of:

        F1(Kr,Reff) = F(Reff) - Kr
        F2(Kr,Reff) = G(Kr) - Reff

        J = [ dF1/dKr dF1/dReff   = [ -1     dF/dReff
              dF2/dKr dF2/DReff ]     dG/dKr    -1   ]
        '''

        # Sum X (initial guess) with total change in X (del_X)
        Kr_in = X[0] + del_X[0]
        Reff_in = X[1] + del_X[1]

        # constants per load cycle:
        a = cdata[0] # current crack size 
        Krw = cdata[1] # Willenborg's Kreduction 
        DKi = cdata[2] # current DeltaK
        Kmax = cdata[3] # current Max K
        Kmin = cdata[4] # current Min K

        # set diag of Jacob:
        Jacob[0][0] = -1.0
        Jacob[1][1] = -1.0

        # Compute dF/dReff:
        derivDKth = self.dDKthbydReff(a,Reff_in)
        A = -1.0*(Krw)/(DKi*(self.Rso - 1.0))
        Jacob[0][1] = A*derivDKth

        # Compute dG/dKr:
        Jacob[1][0] = (-DKi)/((Kmax-Kr_in)**2)

##        # do finite difference to check dF/dReff:
##        # dF/dReff
##        inc = 0.00001
##        f1_val = Numeric.zeros((2),Numeric.Float64)
##        self.RK_Func(cdata,del_X,X,f1_val)
##        f2_val = Numeric.zeros((2),Numeric.Float64)
##        self.RK_Func(cdata,del_X,(X[0],X[1]+inc),f2_val)
##        Fdif = f2_val-f1_val
##        fbyr = Fdif[0]/inc
##
##        # dG/dKr
##        f1_val = Numeric.zeros((2),Numeric.Float64)
##        self.RK_Func(cdata,del_X,X,f1_val)
##        f2_val = Numeric.zeros((2),Numeric.Float64)
##        self.RK_Func(cdata,del_X,(X[0]+inc,X[1]),f2_val)
##        Fdif = f2_val-f1_val
##        gbyk = Fdif[1]/inc
##
##        Jacob[0][1] =fbyr
##        Jacob[1][0] =gbyk
##
##        print "Jacob(0,1) = %.8f" %(Jacob[0][1])
##        print "Finite diff = %.8f" %(fbyr)
##        print "Jacob(1,0) = %.8f" %(Jacob[1][0])
##        print "Finite diff = %.8f" %(gbyk)
##        x

####### Willenborg

    def __Reff(self,Rso,Kred,Kmin,Kmax,R,a):
        '''
        Method to calculate Reff
        '''

        DKi = Kmax - Kmin
        DKth = self.Calc_DKth(a,R)
        phi = (1.0 - DKth/DKi)/(Rso - 1.0)
        nKred = phi * Kred
        return (Kmin - nKred)/(Kmax - nKred),DKth

def FiniteDiff(material):
    '''
    just a little function to do a little verification of the jacobian.
    '''
    fatigue1=Willenborg(material,0)
##    fatigue1.PrintMaterialModel()
##    print "Kic", fatigue1.Kic,"\n"
    a = 0.000271260083761
##    a = 0.1
    DKi = 1.5
    N = 10
    # check dDKth/dR by comparison of finite difference with analytical
    # to help debug Willenborg()...
    inc=0.1

    for i in xrange(6): 
        R = 0.0
        DKth2=fatigue1.Calc_DKth(a,R+inc)
        DKth1=fatigue1.Calc_DKth(a,R)
        deriv = fatigue1.dDKthbydReff(a,R)
        findiff = (DKth2-DKth1)/inc
        err = 100*((deriv-findiff)/deriv)
        print "inc =", inc
        print "R: %1.1f analytical = %f, finite diff = %f, %% error = %f" % \
              (R,deriv,findiff,err)

        R = -0.11718140 
        DKth2=fatigue1.Calc_DKth(a,R+inc)
        DKth1=fatigue1.Calc_DKth(a,R)
        deriv = (fatigue1.dDKthbydReff(a,R))
        findiff = ((DKth2-DKth1)/inc)
        err = 100*((deriv-findiff)/deriv)
        print "R: %1.1f analytical = %f, finite diff = %f, %% error = %f" % \
              (R,deriv,findiff,err)


        R = 0.4
        DKth2=fatigue1.Calc_DKth(a,R+inc)
        DKth1=fatigue1.Calc_DKth(a,R)
        deriv = fatigue1.dDKthbydReff(a,R)
        findiff = (DKth2-DKth1)/inc
        err = 100*((deriv-findiff)/deriv)
        print "R: %1.1f analytical = %f, finite diff = %f, %% error = %f" % \
              (R,deriv,findiff,err)


        print "\n"
        inc=inc/1000.0

####### end Willenborg

if __name__== '__main__':

##    material = [84,75,37,27,1.0,1,0.209e-7,2.947,0.5,1,2,2.0,0.1,0.7,1.9, \
##                 0.3, 0, 0.0015] # Al 7075 T-6
##    material=[84,75,57,50,1.0,1,0.209e-4,1.9,0.5,1,2,2.0,0.1,0.7,1.9, \
##                 0.3, 0, 0.00015]
    #  materials**    = [UTS, YS, Kie, Kic, Ak, Bk, C, n, p, q, DKo, \
    #                    Cthp, Cthm, Rcl, alpha, Smax/SIGo, R, a0, np, damp]
##    material=[84.0, 75.0, 17.0, 7.8, 1.0, 1.0, 0.209e-4, 1.9, 0.5, 1.0, 2.0, \
##              2.0, 0.1, 0.7, 1.9, 0.3, 0.7, 0.0015, 10.,100.0]
    # SIPS3002
##    material=[85.0, 76 .0, 38.0, 28.0, 1, 1, 0.233e-7, 2.885, 0.5, 1.0, 3.0, \
##                   2.0, 0.1, 0.7, 1.9, 0.3, 0., 0.00015, 3.5, 1000]
    # for sips coupon01: 7075-T651 (NASGRO 4.2)
    # [UTS, YS, Kie, Kic, Ak,  Bk, C,     n,   p,  q,  DK1, Cthp, Cthm, Rcl, \
    #                                alpha, Smax/SIGo, R, a0,     np, damp Rso]
    material=[85.0,75.0,38.0,28.0,0.75,2.0,3.0e-8,2.80,0.5,1.0,0.7,1.3,0.1, \
             0.7,2.0,0.3,0.,0.00015,3.5,1000,3.0]

    fat1=Willenborg(material,0)
    fat2=PYDdadN.dadN(material,0)
    fat1.SetKol()
    fat2.SetKol()

    a = 0.000271260083761
    DK = [10.5,13.0,10.5]
    R = [0.15,0.0,0.15]
    N = 10

    for i in range(len(DK)):

        print " no retardation:"
        print fat1.Calc_dadN(DK[i],a,int(N),R[i]), \
              fat2.Calc_dadN(DK[i],a,int(N),R[i])
        print "**"
        print " Willenborg retardation:"
        print fat1.Willenborg(DK[i],a,int(N),R[i]), \
              fat2.Compute_dadN(DK[i],a,int(N),R[i])
        print "** ** ** ** \n"

        N+=1
        a+=fat2.Compute_dadN(DK[i],a,int(N),R[i])













