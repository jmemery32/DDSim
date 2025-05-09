import os, sys, math, numpy as Numeric
import numpy.oldnumeric.linear_algebra as LinearAlgebra
import cPickle

# wash
import Vec3D
import ColTensor
import MeshTools
import JohnsVectorTools as JVT
import dadN as dadN
# John D
import GeomUtils
import Integration
# me
import Statistic, DamClass, DamHistory, DamErrors

# some global functions that come in handy (usually) when debugging... 

def holdit():
    '''
    a little function for debugging.
    '''
    print ' '
    print ' ################################################'
    print "         ---> hit enter to move on <---"
    print ' ################################################'
    print ' '
    c = sys.stdin.read(1)

def DumpFile(xlist,ylist,filename):
    '''
    function to dump items in list out to a text named filename.
    '''

    fil=open(filename,'a')
    for i in range(len(xlist)):
        fil.write(str(xlist[i])+' '+str(ylist[i])+"\n")

    fil.write("# *** \n")
    fil.close()

##############################################################################
##############################################################################

class Damage:
    '''
    A parent class to Fellipse, Hellipse, Qellipse.  Does some basic things...
    '''

    def __repr__(self):
        '''
        Print Damage info to screen.
        '''
        s="\n"
        s+=' ---------------- Damage Element: '+self.names[0]+" -------\n"
        s+=' doid          = %i' % (self.doid)+"\n"
        s+=' Current center: '+repr(self.Trans)+"\n"
        s+=" deltaN list   : \n "+repr(self.dN)[1:-1]+"\n"
        for i in range(self.CrackFrontPoints):
            s+=' '+self.names[i+1]+'      = '+repr(self.a[i])+"\n"
        s=s[0:-1]
        s+=' ---------------- End Damage Element: '+self.names[0]+" -------\n"
        s+="\n"
        return s 

######## Damage

    def AppendInterpolateLifeLists(self,As,types,as,dNs):
        length=len(self.a[0][1])
        for i in range(length):
            if self.dN[i] >= 0.0:
                ab=[self.a[k][0][i] for k in range(self.CrackFrontPoints)]
                as+=[ab[0]]
                As+=[self.CalcArea(ab)]
                dNs+=[self.dN[i]]
                if self.names[0]=='Fellipse': types+=[0]
                elif self.names[0]=='Hellipse': types+=[1]
                else: types+=[2]

######## Damage

    def ParentComputeFrontPoints(self,a,b,center,rotation,yzfunc):
        '''
        to compute Vec3D objects that describe the crack front in global coords
        for the state at which a,b,center and rotation describe the crack.
        '''

        phis=[-1.0+float(i)*0.1 for i in range(20)]
        qpnts=[]

        for phi in phis:
            y,z = yzfunc(a,b,1.0,phi)
            u0 = rotation[1][0]*y+rotation[2][0]*z
            u1 = rotation[1][1]*y+rotation[2][1]*z
            u2 = rotation[1][2]*y+rotation[2][2]*z
            qpnt=center+Vec3D.Vec3D(u0,u1,u2)
            qpnts += [center+Vec3D.Vec3D(u0,u1,u2)]

        return qpnts

######## Damage

    def Lists(self,ai,af,x,K,tol):
        '''
        Return appropriate lists of a and K within the interval [ai,af] to be
        integrated in .Simpson().

        tol = a tolerance to compare ai to af and determine if the interval is
        effectively of zero length.  
        '''

        # if x and K are len() == 0 or 1 return empty lists and neglect
        if len(x) < 2: return [],[]

        # create newx by adding ai+x[start:end]+af
        z = Numeric.greater(x,ai)
        z = z.tolist()
        # sometimes when the damage transitions, ai can be larger than any a's
        # in x.  If this happens return empty lists and, in Simpson, return 0
        # when an empty list is passed in.

        # also... 
        # if Xlist and Klist are return None, we must be trying to integrate a
        # a damage element that does not contain data on the interval [ai,af].
        # This is possible since we integrate Fellipse, Hellipse, Qellipse, in
        # that order and can occur when calculating life during a monte sim 
        # when calc'ng Nminus... this used to be taken care of differently, but
        # now the ValueError exception covers this too.  
        try: 
            start = z.index(1)
        except ValueError:
            return [],[]

        if af == 'end' or af==x[-1]:
            end = len(x)-2
            af = x[-1]
        else:
            z = Numeric.greater(x,af)
            z = z.tolist()
            end = z.index(1)-1

        # If both ai and af are smaller than x[0], the interval of integration
        # does not have any values.  This could happen during monte simulation
        # if the crack won't grow for the smallest ai or if the crack
        # necessarily get's a little bigger from transition.  returning the 
        # following, upon integration, should amount to zero contribution to
        # life or Nminus... 
        if abs(ai - x[0]) < tol and abs(x[0] - af) < tol:
            return [ai,af],[K[0],K[0]]

        # check to see if ai and af are the same.  If they are make Ki and Kf
        # the same
        if abs(ai - af) < tol:
            Ki=K[start-1]
            Kf=Ki

        # Linearly interpolate for Ki and Kf
        else: 
            Ki = K[start-1] + (ai-x[start-1])* \
                 ((K[start]-K[start-1])/(x[start]-x[start-1]))
            try: 
                Kf = K[end] + (af-x[end])* \
                     ((K[end+1]-K[end])/(x[end+1]-x[end]))
            except IndexError, message:
                print x, K, ai, af
                print len(x), len(K)
                raise DamErrors.ListsIndexError, (ai,af,start,end,message)
            except ZeroDivisionError:
                Kf = 0.0

        # if start == 0 it must be that due to transition x[0] is a little
        # larger than the initial a, or this is the second, third, damage
        # element in a sequence
        if start == 0:
            newx=x[start:end+1]
            newx+=[af]
            newK=K[start:end+1]
            newK+=[Kf]
        else: 
            newx=[ai]
            newx+=x[start:end+1]
            newx+=[af]
            newK=[Ki]
            newK+=K[start:end+1]
            newK+=[Kf]

##        print 'ai,af', ai,af
##        print 'newx',newx
##        print 'newK', newK
##        holdit()

        return newx,newK

######## Damage

    def NewNorm(self,ya,za,yb,zb,y00,z00,TRotation):
        '''
        compute the new normal in the plane of the ellipse.

        (ya,za) - coordinate pairs wrt coord. system of ellipse demarking the
                  the location of the intersectio of the ellipse with the body
                  in the general y' direction.
        (yb,zb) - ditto, z' direction.
        (y00,z00) - location of new vertex of ellipse
        TRotation - the ellipse's rotation matrix (transposed)
        '''

        v1=Vec3D.Vec3D(0., ya-y00, za-z00)
        v1=v1.Normalize()
        v2=Vec3D.Vec3D(0., yb-y00, zb-z00)
        v2=v2.Normalize()
        v_ave=(v1+v2).Normalize()
        new_norm=-1.0*v_ave
        new_norm=Numeric.dot(TRotation, \
                         [new_norm.x(),new_norm.y(),new_norm.z()])
        new_norm=Vec3D.Vec3D(new_norm[0], new_norm[1], new_norm[2]).Normalize()

        return new_norm

######## Damage

    def PointOutside(self,qpnt,model):
        '''
        sometimes MestTools needs help answering this question if just the
        combination of the query point and poor element shape arises.
        '''

        # trying this again... 7/3/06
        # for New MeshTools.pyd (7/3/06).
        # -3  point is definitely outside (empty range tree)
        # -2  point is inside mesh
        # -1  point is outside, but can calculate shortest distance to nearest
        #     surface point 
        # >=0 id of an element for which it is not clear, FindNaturalCoords
        #     failed, i.e. exceeded max iterations
        for x in range(4): 
            (flag,dist)=model.IsPointOutsideMesh(qpnt)
            if flag == -3: return (1,0.0)
            elif flag == -2: return (0,0.0)
            elif flag == -1: return (1,dist)
            else: # slightly permute qpnt
                blah = open('Qpnt_failed.txt','a')
                blah.write(str(flag)+' '+str(qpnt.x())+' '+str(qpnt.y())+' '+ \
                           str(qpnt.z())+"\n")
                blah.close()
                qpnt=1.0001*qpnt

        # if we get here, MeshTools can't determine if the point is in or out.
        # return some bogus stuff and deal with it later.

        return ('no point','ever')

######## Damage

    def SeekLongDimension(self,vec,Rotation,model,Trans):
        '''
        a method that finds the approximate length of the structure (stored in
        model) in the direction of vec.
        vec is a vector wrt the prime basis system given by the rows of
        Rotation.  Trans is a Vec3D obj that gives the location of (0,0,0) of
        the primed basis system.

        There are two problems with this:
        1. if the corner that we are examining is less than a 90 internal angle
           we may never find a point inside the model
        2. it is concievable that we get wrong lengths if we shoot across a
           hole... 
        '''

        Max,Min=model.GetMaxDimension()
        m=1.01*((Max-Min).Magnitude())
        a=0.5*(min(self.GiveCurrent()))
        TRotation=Numeric.transpose(Rotation)

        u=Numeric.dot(TRotation,[(m*vec).x(),(m*vec).y(),(m*vec).z()])
        Gu=Trans + Vec3D.Vec3D(u[0],u[1],u[2])
        (flagm,dist)=self.PointOutside(Gu,model)
        u=Numeric.dot(TRotation,[(a*vec).x(),(a*vec).y(),(a*vec).z()])
        Gu=Trans + Vec3D.Vec3D(u[0],u[1],u[2])
        (flaga,dist)=self.PointOutside(Gu,model)

        # if for some reason neither point is in the model, or neither point is
        # outside the model, return the longest dimension of the bounding box! 
        if flaga==flagm: return m

        for i in range(20):
            am=0.5*(a+m)
            # if a and m are with 1% or their average, break... 
            if (m-a)/am < 0.01: break
            u=Numeric.dot(TRotation,[(am*vec).x(),(am*vec).y(),(am*vec).z()])
            Gu=Trans + Vec3D.Vec3D(u[0],u[1],u[2])
            (flag,dist)=self.PointOutside(Gu,model)
            if flag == flagm: m=am
            else: a=am

        return 0.5*(m+a)

######## Damage

    def Simpson(self,x,K,Cumulative_Life,dadN,tol):
        '''
        A function to integrate the K vs. a history, by simpson's rule.  x[0]
        and x[-1] are the bounds of the integral (and K matches), see .Lists()
        above.  

        x -  a list of the "x" components
        K -  NOT "y", is the SIF corresponding to x
        tol = a tolerance to compare ai to af and determine if the interval is
        effectively of zero length.  
        '''

        # sometimes when the damage transitions, ai can be larger than any a's
        # in x.  If this happens return empty lists and, in Simpson, return -1
        # when an empty list is passed in.
        if len(x) == 0: return 0.0

        # if ai = af, within some tolerance, return 0.0
        # else check to see if K(ai) > Kth, if not, return -1
        Reff=dadN.CalcReff(int(Cumulative_Life),9999)
        Kth=dadN.Calc_DKth(x[0],Reff)

        if K[0] < Kth: return -1
        elif abs(x[0] - x[-1]) < tol: return 0.0

        # Simpson's rule...
        # notice in calculation for da/dN, Life is a required arg.
        This_Life=Cumulative_Life
        count=0
        x0=x[count]
        K0=K[count]
        while (len(x)-count) > 2:
            x1=x[count+1]
            x2=x[count+2]
            K1=K[count+1]
            K2=K[count+2]
            count+=2
            # for cases where K is small (e-10), dadN can become zero
            # handle exception by returning -1 from Simpson (analgous to
            # returning -1 in Fwd_integration when K < Kth.
            # **** This may have to be rethought when variable amplitude load
            # is applied.  
            try: y0= 1.0/(dadN.Calc_dadN(K0,x0,int(This_Life)))
            except ZeroDivisionError:
                return -1

            try: y1= 1.0/(dadN.Calc_dadN(K1,x1,int(This_Life)))
            except ZeroDivisionError:
                return -1

            try: y2= 1.0/(dadN.Calc_dadN(K2,x2,int(This_Life)))
            except ZeroDivisionError:
                return -1

            # weights computed with Mathematica... 
            w0=-1.0*(2.0*x0*x0-x0*x2-x2*x2-3.0*x0*x1+3.0*x1*x2) / \
                    (6.0*(x0-x1))
            w1=(-1.0*x0*x0*x0+3.0*x0*x0*x2-3.0*x0*x2*x2+x2*x2*x2) / \
               (6.0*(x0-x1)*(x1-x2))
            w2=-1.0*(x0*x0+x0*x2-2.0*x2*x2-3.0*x0*x1+3.0*x1*x2) / \
                    (6.0*(x2-x1))

            Add=w0*y0+w1*y1+w2*y2

            # don't subtract life if da/dN comes back negative (by accident)
            if Add > 0: This_Life += Add

            # Break while loop if K1 or K2 exceed Kic. 
            if max(K1,K2) > dadN.Kic: return This_Life - Cumulative_Life 

            # change the guard... 
            x0=x2
            K0=K2

        # Get the tail end if not evenly divisible in to 3's 
        if (len(x) - count) == 2:
            x1=x[count+1]
            K1=K[count+1]
            # must calc y0 again, incase len(x) was never > 1!
            try: y0= 1.0/(dadN.Calc_dadN(K0,x0,int(This_Life)))
            except ZeroDivisionError:
                return -1
            try: y1= 1.0/(dadN.Calc_dadN(K1,x1,int(This_Life)))
            except ZeroDivisionError:
                return -1
            This_Life+=0.5*(y0+y1)*(x1-x0)

        return This_Life-Cumulative_Life

######## Damage

    def Increment(self,Klist,alist,dadNobj,r,fakeN):
        '''
        method to calculate the increment to expand ellipses for developement
        of K v. a curve.

        r = percentage of growth (like for 10% r = 1.10) 
        '''

        RANGE=range(len(Klist))

        # assumes R = 0.; if Kmax>Kic, use paris... 
        if max(Klist) > max([dadNobj[i].Kic for i in RANGE]):
            inclist = [dadNobj[i].Paris(Klist[i]) for i in RANGE]
        else:
            inclist = [dadNobj[i].Calc_dadN(Klist[i],alist[i],int(fakeN),0.0) \
                               for i in RANGE]

        if max(inclist) <= 0.0:
            inclist = [dadNobj[i].Calc_dadN(Klist[i],alist[i],0,0.0) \
                               for i in RANGE]

        maxinc=max(inclist)
        inc=[]
        for rate in inclist:
            f = rate/maxinc
            inc += [(f*(r-1.0)+1.0)] # yes. this is right

        return inc # inc[i] > 1.0...

######## Damage

    def ElliptInt(self,ko2):
        '''
        calc ellip. int. from abramowitz infinte sum.
        '''

        change1=1.0 # % error for K
        change2=1.0 # % error for E
        K_old=math.pi/2.0
        E_old=math.pi/2.0
        count=1.0
        coeff=1.0

        while abs(change1) > 1.0e-10 and abs(change2) > 1.0e-10:  
            coeff=(((2.0*count-1)/(2.0*count))**2)*coeff
            K=K_old+(math.pi/2)*coeff*(ko2**count)
            E=E_old-(math.pi/2)*coeff*(ko2**count)/(2.0*count-1)
            change1=(K-K_old)/K_old
            change2=(E-E_old)/E_old
            K_old=K
            E_old=E
            count=count+1

        return E,K

######## Damage

    def ExtendEllipse(self,r,fakeN):
        '''
        Called from BuildKva after computing K's for the ellipse.  This
        method is used to increase the dimensions of the ellipse.
        '''

        RANGE=range(self.CrackFrontPoints)

        alist = [self.a[i][0][-1] for i in RANGE]
        Klist = [self.a[i][1][-1] for i in RANGE]
        inc = self.Increment(Klist,alist,self.dadN,r,fakeN)
        anew = [inc[i]*alist[i] for i in RANGE]
        self.UpdateState(anew)

######## Damage

    def AverageSurfNorm(self,model,nid,sig1,angletol): # jme #@
        '''
        Queries the surface mesh and calculates the surface normals at the
        node, nid.  The normal is given as the average of the UNIQUE surface
        normals at a node.  The average normal is returned along with the
        standard deviation of the surface norms and a boolean that tells if
        any one of the surface norms is in the direction of sig1

        nid is the surface node we are comparing the surface normal at
        sig1 is a unit vector in the direction of the first principal stress.
        '''

        seids=model.GetAdjacentSurfElems(nid)

        # loop over the elements adjacent to nid and build a list of
        # their surface normals, norms.  
        norms=[]; compare = 0
        for seid in seids:
            n=model.SurfaceNormal(nid,seid)
            norms+=[n]
            cos=sig1*n
            if abs(cos) >= angletol: 
                compare = 1

        # Calc. average surface normal... 
        # loop over norms and keep unique normals
        zero=0.01 # to compare the dot products to
        norm=Vec3D.Vec3D(0,0,0) # to average
        count=0.0 # to keep track of how many normals are added to norm
        x_s=0.0; y_s=0.0; z_s=0.0 # for standard dev
        RANGE=range(len(norms))
        for i in RANGE:
            flag=0 # add to this if duplicate normal
            for j in range(i+1,len(norms)):
                dot=norms[i]*norms[j]
                if abs(dot) >= angletol:
                    flag+=1
            if flag == 0: # i.e. norms[i] is not the same as any other in the
                          # bottom half of the list
                norm+=norms[i]
                count+=1.0

        # mean
        norm=(norm/count)
        unorm=norm.Normalize()

        # Find the minimum angle between n's and norm to use for the
        # comparison if rotation needs to be constructed.
        COS=0.0 # i.e. angle = 90.0

        # Calc. standard deviation in surface normal... 
        # some other things used in the std dev loop... 
        count=0.0 # to keep track of how many normals are added to norm
        x_s=0.0; y_s=0.0; z_s=0.0 # for standard dev
        num=0 # a counter to loop through norms in the inner loop
        for i in RANGE:
            flag=0 # add to this if duplicate normal
            num+=1
            cos=norms[i]*unorm
            if cos > COS: COS = cos

            for j in range(num,len(norms)):
                dot=norms[i]*norms[j]
                if (1.0-abs(dot))<=zero:
                    flag+=1

            if flag == 0: # i.e. norms[i] is not the same as any other in the
                          # bottom half of the list
                count+=1
                x_s+=(norms[i].x()-norm.x())**2.0
                y_s+=(norms[i].y()-norm.y())**2.0
                z_s+=(norms[i].z()-norm.z())**2.0
        # standard dev.
        x_dev=((x_s/count))**0.5
        y_dev=((y_s/count))**0.5
        z_dev=((z_s/count))**0.5

        return unorm, max(x_dev,y_dev,z_dev),compare,COS

######## Damage

    def Intersect(self,y1a,z1a,y1b,z1b,y2a,z2a,y2b,z2b,model,a,b,TRotation, \
                  Trans): #@
        '''
        Given the coordinates of four points, returns the point where the two
        lines intersect, the "average" outward normal, the angle between the
        two lines and the parameter "check".  "check" is used to determine if
        the angle (which, by method calc'd here is always less than 180)
        is swept out within the volume (check = 1), or if it is swept out in
        thin air (check = 0).

        returns
        (y0, z0) - hopefully the intersection between external surfaces the
                   ellipse intersects
        norm - an outward normal of the intersection (in GLOBAL COORDS).  NOTE:
               this normal is necessarily in the plane of the ellipse.
        angle - the angle between the intsection lines
        check - a boolean to indicate if it is a convex or concave surface
                intersection
        '''

        # Solving this linear system as simultaneous equations give control
        # over the special cases (in theory!)... 

        # to compute the slope, m = ned/donkey
        ned1=z1a-z1b; donkey1=y1a-y1b; ned2=z2a-z2b; donkey2=y2a-y2b;

        # Define ZERO for this function.  Denominators less than this are bad!
        # also used below to compare (v1+v2).Magnitude()...should choose zero
        # relative to initial crack size
        min_crack=min(a,b)
        zero=(1e-10)*min_crack

        # if this damage element is associated with a damage origin that is a 
        # surface node, save ourselves the if/ands, set y0, z0 = 0.  This was
        # already assumed in Fellipse.__FourPoints.  Surface cracks that need 
        # to change the location of the vertix are picked off in
        # Fellipse.__FtoQ() or Hellipse.WillGrow().
        if self.is_surf:
            y0,z0=0.0,0.0

        # check the components of the two lines.
        # when both slopes are NOT infinite...
        elif abs(donkey1) > zero and abs(donkey2) > zero:
            # y-intercepts... 
            b1=-1.0*(ned1/donkey1)*y1a+z1a
            b2=-1.0*(ned2/donkey2)*y2a+z2a

            # if m1 = m2 (slopes) and b1 != b2 (my rule is kind of hoke) the 
            # lines are parallel 
            if abs((ned1/donkey1)-(ned2/donkey2)) <= zero and \
               abs(b1-b2) > (1e5)*zero:
                # the idea here is that if the lines are parallel we have a
                # through crack.  WillGrow should catch this, so return:
                #      y0 , z0,norm              ,angle,check
                return 0.0,0.0,Vec3D.Vec3D(0,0,0),  0.0,0 

            # if m1 = m2 and b1 = b2 the lines are the same line
            elif abs((ned1/donkey1)-(ned2/donkey2)) <= zero:
                #y0=0.0 # old
                #z0=0.5*(b1+b2) # old
                # return the coords of a perpendicular from this line to the
                # center of the ellipse. (jme 6/9/05)
                b=0.5*(b1+b2)
                m=0.5*((ned1/donkey1)+(ned2/donkey2))
                y0=b*m/(1.0-m*m)
                z0=(-1.0/m)*y0

            # otherwise, solve for their intersection point
            else:
                y0=(b2-b1)/((ned1/donkey1)-(ned2/donkey2))
                z0=(ned1/donkey1)*y0+b1

        # slope of line 1 is infinite 
        elif abs(donkey1) <= zero and abs(donkey2) > zero:
            b2=-1*(ned2/donkey2)*y2a+z2a
            y0=0.5*(y1a+y1b)
            z0=(ned2/donkey2)*y0+b2

        # slope of line 2 is infinite
        elif abs(donkey1) > zero and abs(donkey2) <= zero:
            b1=-1*(ned1/donkey1)*y1a+z1a
            y0=0.5*(y2a+y2b)
            z0=(ned1/donkey1)*y0+b1

        # if we've come this far and have not determined y0, and z0, then both
        # lines are parallel with z-axis.  So, two cases are possible:
        # 1.) they are the same line
        # 2.) They are parallel with each other
        
        # If 1.) is true the following elif statement should be true.  Assume
        # y0, z0 are average of 4 points.  
        elif (y1a/y2a-1.) <= 0.01 or (y1b/y2b-1.) <= 0.01:
            #y0=0.25*(y1a+y1b+y2a+y2b)
            #z0=0.25*(z1a+z1b+z2a+z2b)
            # return the coords of a perpendicular from this line to the
            # center of the ellipse. (jme 6/9/05)
            y0=0.25*(y1a+y1b+y2a+y2b)
            z0=0.0

        # If 2.) is true than that's what we have left for the else.  Raise an
        # exception that we'll catch in Fellipse.WillGrow(). 
        else:
            # the idea here is that if the lines are parallel we have a
            # through crack.  WillGrow should catch this.
            #      y0 , z0,norm              ,angle,check
            return 0.0,0.0,Vec3D.Vec3D(0,0,0),0.0  ,0

##        # if y0 and z0 are really close to zero, set them equal to zero, the
##        # center of the crack does not move.
##        if abs(y0) < 0.001:
##            y0 = 0.0
##        if abs(z0) < 0.001:
##            z0 = 0.0

        # Now compute the angle between the lines v1 and v2 point away from the
        # intersection of the line.
        # i don't want v_ave (below) to be unit length because for thin
        # structures that might put v_check (below) out of the volume.
        # So, record the lengths of v1 and v2 here and use the minimum to
        # set the length of v_ave... We want v_ave to be bisect the angle
        # between v1 & v2, consequently, after i record their lengths, i
        # normalize them.
        v1=Vec3D.Vec3D(0, y1a-y0, z1a-z0)
        v1_length=v1.Magnitude()
        v1=v1.Normalize()
        v2=Vec3D.Vec3D(0, y2a-y0, z2a-z0)
        v2_length=v2.Magnitude()
        v2=v2.Normalize()
        dot = v1*v2

        if dot < -1.0:
            dot = -1.0
        elif dot > 1.0:
            dot = 1.0
        angle=(math.acos(dot))*(180.0/math.pi)

        # sample a point between the two lines to see if it is in the volume
        # the average of the two vectors gives a point between them
        if (v1+v2).Magnitude() > zero:
            v_ave=((v1+v2).Normalize())*min(v1_length,v2_length)

        # if the v1 and v2 are equal and opposite than their sum is zero.  This
        # is problematic when finding their average, so do something else.
        else:
            v_ave=(Vec3D.Vec3D(0,v1.z(),-1*v1.y()).Normalize())* \
                   min(v1_length,v2_length)

        # translate to the intersection point in y & z (plane of ellipse)
        v_check=Vec3D.Vec3D(0,y0,z0)+v_ave

        # rotate into global coords
        u=Numeric.dot(TRotation, \
                      [v_check.x(),v_check.y(),v_check.z()])

        # translate to damage origin
        qpnt=Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        # and check it!
        (flag,dist)=self.PointOutside(qpnt,model)

        # Adopt the convention that the normal points outside the body
        if flag == 0:
            norm=-1.0*v_ave
            check=1
        else:
            norm=v_ave
            check=0

        # and rotate so is in global coord system
        norm=Numeric.dot(TRotation,[norm.x(),norm.y(),norm.z()])
        norm=Vec3D.Vec3D(norm[0], norm[1], norm[2]).Normalize()

        return y0,z0,norm,angle,check

##############################################################################
##############################################################################

class Fellipse(Damage):
    '''
    the center of the ellipse starts at the damage origin and moves
    sigma1 defines the local x axis of the ellipse, normal to the plane of the
    ellipse.  
    sigma2 defines the local y axis (major axis and a)
    sigma3 defines the local z axis (minor axis and b)
    '''

######## Fellipse

    def CalcArea(self,ab): 
        '''
        simple method to calculate the current area of the ellipse.
        '''
        a,b=0.5*(ab[0]+ab[2]),0.5*(ab[1]+ab[3])

        return math.pi*a*b

######## Fellipse

    def CalculateKi(self,scale): 
        a,b=self.GiveCurrent()
        try:
            Ki=self.__CalculateKi(a,b,self.Trans,scale)
            return Ki

        except DamErrors.FitPolyError, message:
            # for cases where the crack grows well outside the model,
            # model.GetPtStress returns 'point not in any elements'.  This is
            # flagged in __FitPoly and the K's are set to -10 (integer).  If 
            # happens, assume the crack has become unstable.  This is different
            # than the next elif where i check list because that checks if the
            # entire crack is outside the mesh.
            if self.verbose:
                print ' FitPoly dumped for', message[0], 'at doid', self.doid,\
                      'due to',message[1] 
                print ' **SIMULATION OK, means crack has out grown its welcome'
            return 4

######## Fellipse

    def ComputeFrontPoints(self,crack_step=-1):
        '''
        to compute Vec3D objects that describe the crack front in global coords
        for the state at which a,b,center and rotation describe the crack.

        crack_step indicates which crack front to return.  the default is the
        current step
        '''

        def CurrentStep():
            a,b=self.GiveCurrent()
            return self.ParentComputeFrontPoints(a,b,self.Trans,self.Rotation,\
                                             self.__yzfunc)

        # if the crakc_step possibly describes a crack state... 
        if crack_step >= 0:
            if crack_step < len(self.a[0][0]):
                a = 0.5*(self.a[0][0][crack_step] + \
                         self.a[2][0][crack_step])
                b = self.a[1][0][crack_step]
                return self.ParentComputeFrontPoints(a,b,self.Trans,\
                                               self.Rotation,self.__yzfunc)

            # if not, return the current state. 
            else: return CurrentStep()
        else: return CurrentStep()

######## Fellipse

    def __CalculateKi(self,a,b,center,scale=1.0): 
        '''
        a helper function for CalcKi that will allow computing Ki's 
        without modifying internal data (JD)
        '''

        switched=0 # a parameter to keep track of K's when b > a
                   # switched = 0 ---> a > b
                   # switched = 1 ---> b > a
        
        if a < b:
            a,b=b,a
            switched=1

        # try __FitPoly.  If no, then points lie outside the model, assume Ki
        # exceeds Kic!  Set the last Ki at 0, 90, & 180 to -10.  Look for this
        # flag in __WillGrow
        try:
            #coef=self.__FitPoly(model)
            coef=self.__FitPolynomial(a,b,center)
        except MeshTools.EmptySearchResult, message:
            raise DamErrors.FitPolyError, ('Fellipse',message)
        except LinearAlgebra.LinAlgError, message:
            raise DamErrors.FitPolyError, ('Fellipse',message)

##        print "biquadratic coefficients: y, z, y*y, y*z, z*z, 1.0"
##        print coef

        p10=coef[0]
        p01=coef[1]
        p20=coef[2]
        p11=coef[3]
        p02=coef[4]
        p00=coef[5]
        (ko2,E,K,E1,E2)=self.__GetOnce(a,b)
        phi = 0.0
        KI=[]
        # loop to move parametric angle from 0 to 90 to 180 to 270...
        for ii in xrange(4): 
            same=((b/a)**0.5)*((a*a*(Numeric.sin(phi))**2+ \
                                b*b*(Numeric.cos(phi)**2))**0.25)
            KI_0=(p00/E)*same
            I11_c=E-E1
            I11_s=E+E1
            KI_1=(2.0/3.0)*same*((p10*a*Numeric.cos(phi)/I11_c)+ \
                                 (p01*b*Numeric.sin(phi)/I11_s))
            I00_c=2.0*E
            I02_c=2.0*E1
            I22_c=E+E2
            I22_s=E-E2
            A0=((p20*a*a/4.0)+(p02*b*b/4.0))* \
               (((5*I00_c*I22_c)-3*I02_c*I02_c)/I00_c)+ \
               ((p20*a*a/2.0)-(p02*b*b/2.0)*I02_c)
            A2=(((p20*a*a/2.0)+(p02*b*b/2.0))*I02_c)+ \
                ((p20*a*a/2.0)-(p02*b*b/2.0)*I00_c)
            B2=p11*a*b/2.0*((I00_c*I22_c-I02_c*I02_c)/I22_s)
            KI_2=(8.0/(15.0*(I00_c*I22_c-I02_c*I02_c)))*same* \
                 (A0+A2*Numeric.cos(2*phi)+B2*Numeric.sin(2*phi))
            KI+=[(KI_0+KI_1+KI_2)*math.sqrt(math.pi)]
            phi=phi+(0.5*math.pi)
        if switched==1:
            KI.append(KI.pop(0))
        for i in range(len(KI)):
            if KI[i] < 0.0 : KI[i] = 0.0

        return JVT.ScalarMult(KI,scale)

######## Fellipse

    def GiveCurrent(self): 
        '''
        Returns the most current a, b.
        '''
        # we want average, and they should all be positive!  
        a = (self.a[0][0][-1] + \
             self.a[2][0][-1])/2.0
        b = (self.a[1][0][-1] + \
             self.a[3][0][-1])/2.0

        assert (a>0. and b>0.)

        return a,b

######## Fellipse

    def GrowDam(self,dN,min_inc,N,max_error,scale,integration): 
        '''
        Calculate individual growth amount for each "corner" of the ellipse
        and advance the crack.
        '''

        assert (len(self.a[1])>0) # i.e there is a SIF History... 

        # change to percentage of minimum ellipse dimension
        max_er=max_error*min(a,neg_a,b,neg_b)
        # change back to absolute tolerance
        # max_er=max_error
        max_er=min(max_er,max_error)

        # make scale an attribute for function dadN purposes, should be ok
        # because it updates it here everytime.
        self.scale=scale

        # internal function to be used with RK methods below... 
        def dAdN(N,size):
            foo = []
            abab= [max(0,size[0]),max(0,size[1]),max(0,size[2]),max(0,size[3])]
            a,b = (abab[0]+abab[2])/2.0 ,(abab[1]+abab[3])/2.0
            Xell = (abab[0]-abab[2])/2.0
            Yell = (abab[1]-abab[3])/2.0
            trans=Numeric.dot(self.TRotation, \
                          [0.0,Xell,Yell])
            localcenter=self.Trans + Vec3D.Vec3D(trans[0],trans[1],trans[2])
            KIs = self.__CalculateKi(a,b,localcenter,self.scale)
            foo.append(self.dadN[0].Calc_dadN(KIs[0], size[0],int(N)))
            foo.append(self.dadN[1].Calc_dadN(KIs[1], size[1],int(N)))
            foo.append(self.dadN[2].Calc_dadN(KIs[2], size[2],int(N)))
            foo.append(self.dadN[3].Calc_dadN(KIs[3], size[3],int(N)))
            if max(foo) <= 0.:
                raise DamErrors.dAdNError, ('Fellipse',\
                                      'All growth rates <= 0.0',foo)
            return foo

        # Get current crack lengths...
        # Don't use GiveCurrent because don't want average a and b...
        a=self.a[0][-1]
        b=b[0][-1]
        neg_a=self.na[0][-1]
        neg_b=self.nb[0][-1]

        if integration=='FWD':
            # if material contains less than 18 (req'd for NASGRO eqn), use
            # Paris
            if len(material)<10:
                CC=material[0]
                alpha=material[1]
                # calc growth rate for each "corner" of the ellipse
                # and use TWO largest (ie. a, and b)
                da_dN=CC*((self.a[1][-1])**alpha)
                neg_da_dN=CC*((self.na[1][-1])**alpha)
                db_dN=CC*((self.b[1][-1])**alpha)
                neg_db_dN=CC*((self.nb[1][-1])**alpha)
            # otherwise use NASGRO eqn.
            else:
                # create dadN model for this instance
                if self.dadN == None:
                    # note: with 0 (below), i hard code Kc = Kic in NASGRO eqn.
                    self.dadN=[dadN.dadN(material,0),dadN.dadN(material,0), \
                               dadN.dadN(material,0),dadN.dadN(material,0)]
                # calc.dadN...
                da_dN = self.dadN[0].Calc_dadN(self.a[1][-1],a,int(N))
                neg_da_dN = self.dadN[1].Calc_dadN(self.na[1][-1],neg_a,int(N))
                db_dN = self.dadN[2].Calc_dadN(self.b[1][-1],b,int(N))
                neg_db_dN = self.dadN[3].Calc_dadN(self.nb[1][-1],neg_b,int(N))

        #JohnD try RK4
        elif integration=='RK4':
            # assume always using NASGRO equations with RK4 scheme...
            junk = Integration.RK4slope_vector(dN,N,[a,b,neg_a,neg_b],dAdN)
            da_dN = junk[0]
            neg_da_dN = junk[2]
            db_dN = junk[1]
            neg_db_dN = junk[3]

        #JohnD try RK5 w/ adaptive step size
        else: # if self.IntegrationScheme=='RK5'
            # assume always using NASGRO equations with RK5 scheme...
            if (self.nextdN == 0.0): self.nextdN = dN

            # the integrations scheme can sometimes become unstable and pass
            # calcki numbers that don't make sense (too big).  Hence, we put in
            # a try and except block.  if RK5 scheme is bunk we do on small
            # step by the fwd euler method in hopes that springs us loose of
            # bad spot.  may need to add a finally or something to this in the
            # future.

            try:
                junk,dN,self.nextdN = Integration.RK5_CKslope_vector(\
                    self.nextdN,N,[a,b,neg_a,neg_b],dAdN,max_er,.05)

            except FitPolyDump, message:
                if self.verbose:
                    print ' switch to one step fwd Euler for', \
                          message[0],'at doid', self.doid, 'due to'
                    print '     ', message[1]
                junk = Integration.Eulerslope_vector(self.nextdN/1000.0,N, \
                                                     [a,b,neg_a,neg_b],dAdN)
                self.nextdN = self.nextdN/1000.0

            except dAdNError, message:
                if self.verbose:
                    print ' switch to one step fwd Euler for', \
                          message[0],'at doid', self.doid, 'due to',message[1]
                junk = Integration.Eulerslope_vector(self.nextdN/1000.0,N, \
                                                     [a,b,neg_a,neg_b],dAdN)
                self.nextdN = self.nextdN/1000.0

            except ValueError, message:
                if self.verbose:
                    print ' switch to one step fwd Euler for Fellipse', \
                          'at doid', self.doid, 'due to'
                    print '     ', message[0] 
                junk = Integration.Eulerslope_vector(self.nextdN/1000.0,N, \
                                                     [a,b,neg_a,neg_b],dAdN)
                self.nextdN = self.nextdN/1000.0

            da_dN = junk[0]
            db_dN = junk[1]
            neg_da_dN = junk[2]
            neg_db_dN = junk[3]

        # Caclulate the new a, b, a-, b-.
        # a, b and neg_a, neg_b never to take negative value (which shouldn't
        # happen) so abs() is used to insure that here. 
        inca=abs(da_dN*dN)
        incb=abs(db_dN*dN)
        neg_inca=abs(neg_da_dN*dN)
        neg_incb=abs(neg_db_dN*dN)
        big = max(inca,incb,neg_inca,neg_incb)

        if integration=='RK5': # don't mess-up the adaptivity
            anew=a+inca
            bnew=b+incb
            neg_anew=neg_a+neg_inca
            neg_bnew=neg_b+neg_incb
        else: 
            # don't grow less than min_inc...
            if big < min_inc:
                dN=min_inc/max(da_dN,db_dN,neg_da_dN,neg_db_dN)
                anew=a+abs(da_dN*dN)
                bnew=b+abs(db_dN*dN)
                neg_anew=neg_a+abs(neg_da_dN*dN)
                neg_bnew=neg_b+abs(neg_db_dN*dN)

            # don't grow more than a,b (average)
            elif max(inca,incb,neg_inca,neg_incb) > \
                 min(0.5*(a+neg_a),0.5*(b+neg_b)):
                dN=min(0.5*(a+neg_a),0.5*(b+neg_b))/ \
                       max(da_dN,db_dN,neg_da_dN,neg_db_dN)
                anew=a+abs(da_dN*dN)
                bnew=b+abs(db_dN*dN)
                neg_anew=neg_a+abs(neg_da_dN*dN)
                neg_bnew=neg_b+abs(neg_db_dN*dN)

            # otherwise, just do what seems natural... 
            else:
                anew=a+inca
                bnew=b+incb
                neg_anew=neg_a+neg_inca
                neg_bnew=neg_b+neg_incb

        self.UpdateState(dN,anew,bnew,neg_anew,neg_bnew)

######## Fellipse

    def IntegrateLife(self,ai,af,N_max,Life,tol):
        '''
        Using simpson's rule, integrate for life between ai and af.  For now,
        returns shortest life from the four possible positions.

        ai - initial crack length
        af - final crack length (usually float, but can be str()='end')
        r - increment factor (ai+1 = r * ai)
        tol - a tolerance to compare ai to af and determine if the interval is
        effectively of zero length.  
        '''

        # check the length of the SIF history... 
        if len(self.a[1]) == 0: return 0.0

        N=N_max
        Xlist,Klist=self.Lists(ai,af,self.a[0],self.a[1],tol)
        NN=self.Simpson(Xlist,Klist,Life,self.dadN[0],tol)
        if NN < N and NN >= 0:
            N = NN
        elif NN == -1:
            N = -1

        Xlist,Klist=self.Lists(ai,af,self.b[0],self.b[1],tol)
        NN=self.Simpson(Xlist,Klist,Life,self.dadN[1],tol)
        if NN < N and NN >= 0:
            N = NN
        elif NN == -1:
            N = -1

        Xlist,Klist=self.Lists(ai,af,self.na[0],self.na[1],tol)
        NN=self.Simpson(Xlist,Klist,Life,self.dadN[2],tol)
        if NN < N and NN >= 0:
            N = NN
        elif NN == -1:
            N = -1

        Xlist,Klist=self.Lists(ai,af,self.nb[0],self.nb[1],tol)
        NN=self.Simpson(Xlist,Klist,Life,self.dadN[3],tol)
        if NN < N and NN >= 0:
            N = NN
        elif NN == -1:
            N = -1

        return N

######## Fellipse

    def UpdateState(self,anew):

        # update DamHistory
        for i in range(4):
            self.a[i][0]+=[anew[i]]
        y0_new=(anew[0]-anew[2])/2.0
        z0_new=(anew[1]-anew[3])/2.0
        # update self.Trans.  Note, adding new position of center,
        # (0,y0_new,z0_new), w.r.t. ORIGINAL damage origin,
        # self.DamHistory['cent'][0]...
        trans=Numeric.dot(self.TRotation, \
                          [0.0,y0_new,z0_new])
        self.Trans=self.Trans + \
                   Vec3D.Vec3D(trans[0],trans[1],trans[2])

######## Fellipse

    def ToFile(self,path_name,did): # change for Gerd
        '''
        function to print information about Fellipse to a file.
        '''

        if len(self.DamHistory[0])>0:
            a_his = open(path_name+'''.a''', 'a')
            _a_his = open(path_name+'''.na''', 'a')
            b_his = open(path_name+'''.b''', 'a')
            _b_his = open(path_name+'''.nb''', 'a')
            dN = open(path_name+'''.dN''', 'a')
            for i in xrange(len(self.DamHistory[0])):
                a_his.write(str(did)+' '+str(self.DamHistory['ab'][i][0])+  \
                                     ' '+str(self.DamHistory[0][i])+"\n")
                _a_his.write(str(did)+' '+str(self.DamHistory['-ab'][i][0])+  \
                                     ' '+str(self.DamHistory[180][i])+"\n")
                b_his.write(str(did)+' '+str(self.DamHistory['ab'][i][1])+  \
                                     ' '+str(self.DamHistory[90][i])+"\n")
                _b_his.write(str(did)+' '+str(self.DamHistory['ab'][i][1])+  \
                                     ' '+str(self.DamHistory[270][i])+"\n")
                try: dN.write(str(did)+' '+str(self.DamHistory['dN'][i])+"\n")
                except IndexError: pass
            a_his.close()
            _a_his.close()
            b_his.close()
            _b_his.close()

        else:
            return

######## Fellipse

    def GeometryCheck(self,CMesh): 
        '''
        Determines if the Fellipse is the appropriate type of damage
        element, based on geometry.

        Returns: dam_elem, have_doubt
        
        dam_elem = 0 if no change in crack geometry is necessary
                   4 if crack has outgrown geometry (net fracture)
                   a new damage element if transition is required

        have_doubt = a parameter to indicate how much the geometry is fudged.
        For example, if the particular damage origin does not fit any type of
        damage exactly, i.e. is partway between a half ellipse and quarter
        ellipse, can use this to raise a flag saying "hey i'm not that
        confident".  I have not done that much with this yet.

        CMesh is a list an object used in GeomUtils...
        '''

        a,b=self.GiveCurrent() #@
        # __BuildPhi to replace __CheckPnts; 11/3/04 uses John D's
        # GeomUtils.pyd
        # __BuildPhi to replace __CheckPnts; 11/3/04 uses John D's
        # GeomUtils.pyd
        try: PhiList,CheckList = self.__BuildPhi(CMesh,a,b) #@
        except DamErrors.BuildPhiError, info:
            inc = 0.01 # does not need to be relative to crack size because
            # Phi is parametric
            PhiList,CheckList = self.__CheckPnts(inc)
            if self.verbose: 
                print " BuildPhi returns 1 or 3 phi's:", info
                print " CheckPnts returns:", PhiList
        # if this is a surface node and __BuildPhi returns len(PhiList)=0
        # GeomUtils may have failed to find the right points so use
        # __CheckPnts to reevaluate
        if self.is_surf:
            if len(PhiList)==0:
                inc = 0.01
                PhiList,CheckList = self.__CheckPnts(inc)

##        print PhiList, CheckList
##        print self.__CheckPnts(0.01)
####        PhiList, CheckList = self.__CheckPnts(0.01)
##        xxx

        # if the whole ellipse is within the volume len(PhiList) = 0,
        # and CheckList[0]==1 use Fellipse, calc. Ki and set will_grow...
        if len(PhiList)==0 and CheckList[0]==1: return 0,0
        # entire ellipse could fall outside (like at a corner) if so, change
        # self.Rotation so that crack is perpendicular to second principal
        # stress and recurse! 
        elif len(PhiList)==0 and CheckList[0]==0:
            self.__ChangeRotationToSecondPrincipal()
            return self.GeometryCheck(CMesh)

        # if len(PhiList) > 4 assume that means the ellipse is larger than its
        # surroundings and is therefore grown as much as possible.  Set
        # will_grow = 2 and quit.
        elif len(PhiList) > 4:
            return 4,0

        # Use PhiList and CheckList to define two lines that intersect the
        # ellipse...
        # relative to smallest dimension of ellipse
        tol=(1e-6)*min(a,b)
        # y and z in ellipse space... 
        y1a,z1a,y1b,z1b,y2a,z2a,y2b,z2b,order = \
                        self.__FourPoints(PhiList,CheckList,tol) #@

##        # to print the points in the global 
##        u=Numeric.dot(self.TRotation,[0.0,y1a,z1a])
##        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
##        print 'v', qpnt.x(), qpnt.y(), qpnt.z()
##        u=Numeric.dot(self.TRotation,[0.0,y1b,z1b])
##        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
##        print 'v', qpnt.x(), qpnt.y(), qpnt.z()
##        u=Numeric.dot(self.TRotation,[0.0,y2a,z2a])
##        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
##        print 'v', qpnt.x(), qpnt.y(), qpnt.z()
##        u=Numeric.dot(self.TRotation,[0.0,y2b,z2b])
##        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
##        print 'v', qpnt.x(), qpnt.y(), qpnt.z()
##        print y1a, z1a
##        print y1b, z1b
##        print y2a, z2a
##        print y2b, z2b
##        xx

        # Now we know two points on two lines that intersect the ellipse.  Use
        # this information to find the intersection of these two lines...
        y0_new,z0_new,norm,angle,check= \
                        self.Intersect(y1a,z1a,y1b,z1b,y2a,z2a,y2b,z2b,\
                                       self.model,a,b,self.TRotation,\
                                       self.Trans) 

##        # to print the location of the new center
##        u=Numeric.dot(self.TRotation,[0.0,y0_new,z0_new])
##        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
##        print 'v', qpnt.x(), qpnt.y(), qpnt.z()
##        qpnt=self.Trans+norm
##        print 'v', qpnt.x(), qpnt.y(), qpnt.z()
##        print "PhiList", PhiList
##        print 'y0_new, z0_new, a, b'
##        print y0_new,z0_new,a,b
##        print 'norm,                         angle, check'
##        print norm,'|',angle,check
##        xxx

        # %%%%%%%%%%%% Begin logic statements for shape change %%%%%%%%%%%%%% #
        #                                                                     #
        # now use the information about the angle between the two lines, and
        # check to determine what type of damage element is most appropriate.
        # if 0<=angle<=15 and check = 1 assume a through crack (through crack
        # checks to see if edge crack is more appropriate).
        if angle >=0. and angle <= 15. and check == 1:
            dam_el=4
            have_doubt=0

        # if 0<=angle<=15 and check = 0 assume a small chunk out of the volume
        elif angle >=0. and angle <= 15. and check == 0:
            # if intersection is out side the ellipse, switch to Hellipse
            if (((y0_new/a)**2)+((z0_new/b)**2)) > 1.:
                # I don't remember how this is possible... (10/5/05)
                if len(PhiList) == 2:
                    dam_el=self.__FtoH(y1a,z1a,y2a,z2a,PhiList[0],PhiList[1], \
                                       norm) #@
                    have_doubt=0

                # need to pass self.__FtoH() coords of two points on line
                # nearest to the center of the ellipse...
                elif ((y1a)**2+(z1a)**2) < ((y2a)**2+(z2a)**2):
                    phi_a=PhiList[order[0]]
                    phi_b=PhiList[order[1]]
                    dam_el=self.__FtoH(y1a,z1a,y1b,z1b,phi_a,phi_b,norm) 
                    have_doubt=0
                else:
                    phi_a=PhiList[order[2]]
                    phi_b=PhiList[order[3]]
                    dam_el=self.__FtoH(y2a,z2a,y2b,z2b,phi_a,phi_b,norm)
                    have_doubt=0

            # else, use Fellipse (because only a small portion of the ellipse
            # is cut out) and return 1 for have_doubt.
            else: return 0,0

        # if 15<angle<=75 and check = 1, switch to Qellipse
        # increased upper limit to 125 (from 75) 1-22-05 commented
        # __FtoQ for 75<angle<=125 below ... 
        elif angle > 15. and angle <= 125. and check == 1:
            if len(PhiList) == 4: 
                dam_el=self.__FtoQ(y1a,z1a,y2a,z2a,y0_new,z0_new, \
                               PhiList[order[0]],PhiList[order[2]],norm) 
                have_doubt=1
            else:
                dam_el=self.__FtoQ(y1a,z1a,y2a,z2a,y0_new,z0_new, \
                               PhiList[0],PhiList[1],norm) #@
                have_doubt=1

        # if 15<angle<=165 and check = 0 a slightly larger chunk out of volume
        # decreased upper limit to 145 (from 165) 5-16-05
        # from 145 to 125 (jme 10/5/05)
        elif angle > 15. and angle <= 125. and check == 0:

            # if intersection is outside the ellipse, switch to Hellipse
            if (((y0_new/a)**2)+((z0_new/b)**2)) > 1.:

                if len(PhiList) == 2:
                    dam_el=self.__FtoH(y1a,z1a,y2a,z2a,PhiList[0],PhiList[1], \
                                       norm)
                    have_doubt=0

                # need to pass self.__FtoH() coords of two points on line
                # nearest to the center of the ellipse...
                elif ((y1a)**2+(z1a)**2)>((y2a)**2+(z2a)**2):
                    phi_a=PhiList[order[0]]
                    phi_b=PhiList[order[1]]
                    dam_el=self.__FtoH(y1a,z1a,y1b,z1b,phi_a,phi_b,norm)
                    have_doubt=0
                else:
                    phi_a=PhiList[order[2]]
                    phi_b=PhiList[order[3]]
                    dam_el=self.__FtoH(y2a,z2a,y2b,z2b,phi_a,phi_b,norm)
                    have_doubt=0

            # 5/11/04 - use fellipse for this case
            else: return 0,1

        # if 165<angle<=180, switch to Hellipse
        # 5-11-04 change upper limit to 125: force to use Hellipse
        elif angle > 125. and angle <= 181.:

            if len(PhiList) == 2:
                dam_el=self.__FtoH(y1a,z1a,y2a,z2a,PhiList[0],PhiList[1], \
                                   norm)
                have_doubt=0

            # In this case it is unclear why passing self.__FtoH() coords of
            # two points on line nearest to the center of the ellipse will hurt
            # and since it seems like an easy way to pick the two points passed
            # to FtoH, why not!... THIS MAY BE DANGEROUS!
            elif ((y1a)**2+(z1a)**2)>((y2a)**2+(z2a)**2):
                phi_a=PhiList[order[0]]
                phi_b=PhiList[order[1]]
                dam_el=self.__FtoH(y1a,z1a,y1b,z1b,phi_a,phi_b,norm)
                have_doubt=0
            else:
                phi_a=PhiList[order[2]]
                phi_b=PhiList[order[3]]
                dam_el=self.__FtoH(y2a,z2a,y2b,z2b,phi_a,phi_b,norm)
                have_doubt=0

        return dam_el,have_doubt

######## Fellipse

    def __ArcLength(self,phi_1,phi_2): 
        '''
        A nine point Gauss quadrature that calculates the arc length on the
        ellipse between phi_a and phi_b.
        '''

        # get current major and minor axis dimensions
        a,b=self.GiveCurrent()

        # use GeomUtils
        x,y=self.__yzfunc(a,b,1.0,phi_1)
        if phi_1 > -1.0004 and phi_1 <= -0.5:
            theta_1=math.atan(y/x)
        elif phi_1 > -0.5 and phi_1 <= 0.:
            theta_1=math.pi+math.atan(y/x)
        elif phi_1 > 0.0 and phi_1 <= 0.5:
            theta_1=math.pi+math.atan(y/x)
        else:
            theta_1=2*math.pi+math.atan(y/x)

        x,y=self.__yzfunc(a,b,1.0,phi_2)
        if phi_2 > -1.0003125 and phi_2 <= -0.5:
            theta_2=math.atan(y/x)
        elif phi_2 > -0.5 and phi_2 <= 0.:
            theta_2=math.pi+math.atan(y/x)
        elif phi_2 > 0.0 and phi_2 <= 0.5:
            theta_2=math.pi+math.atan(y/x)
        else:
            theta_2=2*math.pi+math.atan(y/x)

        # i had code that compute the arclength numerically... can be found in
        # v1.3... 
        return GeomUtils.EllipseArcLength(a,b,theta_1,theta_2)

######## Fellipse

    def __BisectPhi(self,phi_L,phi_H,rho,tol):
        '''
        Finds the intersection of the ellipse with the surface of the FEM
        mesh.
        '''

        a,b=self.GiveCurrent()

        # phi_L
        y,z=self.__yzfunc(a,b,rho,phi_L)
        u=Numeric.dot(self.TRotation,[0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        (flag,dist)=self.PointOutside(qpnt,self.model)
        if flag == 1:
            flag_L = 0
        else:
            flag_L = 1

        # phi_H
        y,z=self.__yzfunc(a,b,rho,phi_H)
        u=Numeric.dot(self.TRotation,[0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        (flag,dist)=self.PointOutside(qpnt,self.model)
        if flag == 1:
            flag_H = 0
        else:
            flag_H = 1

        while abs(phi_H-phi_L) > tol:
            phi_M=0.5*(phi_H+phi_L)
            y,z=self.__yzfunc(a,b,rho,phi_M)
            u=Numeric.dot(self.TRotation,[0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            (flag,dist)=self.PointOutside(qpnt,self.model)
            if flag == 1:
                flag_M = 0
            else:
                flag_M = 1

            if flag_M == flag_H:
                phi_H = phi_M
            else:
                phi_L = phi_M

        # flag_H indicates whether avg(phi_L,phi_H) is a transition from or
        # into the volume.  i.e., if flag_H = 0 it is a transition OUT of the
        # volume.  Conversly, if flag_H = 1 it is a transition INTO the volume.
        return phi_L, phi_H, flag_H

######## Fellipse

    def __BisectRho(self,rho_L,rho_H,phi,tol):
        '''
        Finds the intersection of a ray of the ellipse with the surface of the 
        FEM mesh by altering rho while keeping phi constant.
        '''

        a,b=self.GiveCurrent()

        # phi_L
        y,z=self.__yzfunc(a,b,rho_L,phi)
        u=Numeric.dot(self.TRotation,[0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        (flag,dist)=self.PointOutside(qpnt,self.model)
        if flag == 1:
            flag_L = 0
        else:
            flag_L = 1

        # phi_H
        y,z=self.__yzfunc(a,b,rho_H,phi)
        u=Numeric.dot(self.TRotation,[0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        (flag,dist)=self.PointOutside(qpnt,self.model)
        if flag == 1:
            flag_H = 0
        else:
            flag_H = 1

        dist = 1e10 # initialize dist to something large.
        # tol passed in should be something relative to the minimum crack size
        while dist > tol:
            rho_M=0.5*(rho_H+rho_L)
            y,z=self.__yzfunc(a,b,rho_M,phi)
            u=Numeric.dot(self.TRotation,[0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            (flag,dist)=self.PointOutside(qpnt,self.model)
            if flag == 1:
                flag_M = 0
            else:
                flag_M = 1

            if flag_M == flag_H:
                rho_H = rho_M
            else:
                rho_L = rho_M

        # flag_H indicates whether avg(rho_L,rho_H) is a transition from or
        # into the volume.  i.e., if flag_H = 0 it is a transition OUT of the
        # volume.  Conversly, if flag_H = 1 it is a transition INTO the volume.
        return rho_L, rho_H, flag_H

######## Fellipse

    def __BuildPhi(self,CMesh,a,b): 
        '''
        use GeomUtils to build PhiList and CheckList.  (replaces __CheckPnts)

        Returns PhiList and CheckList.  PhiList is a list of parametric angle,
        phi, where the ellipse intersects a surface.  CheckList is a list of
        booleans 0, or 1.  0 means, as you travel CCW around the ellipse you're
        leaving the volume (or mesh), 1 means you're entering it.  i.e. 
        the bool reflects the answer to the question "is my next step in the
        mesh?"

        corn is list of all corner nodes coordinates (Vec3D object) of all
        surface facets.
        '''

        PhiList=[]; CheckList=[]

        # debugging code to rotate the rotation a smidgen so intersections
        # don't fall on edges, which GeomUtils.EllipseCMeshIntersections()
        # is having a hard time with.  9-20-05:
##        teta=0.5*math.pi/180.0
##        self.Rotation=Numeric.dot([[math.cos(teta), math.sin(teta),0.0],\
##                                   [-1.0*math.sin(teta),math.cos(teta),0.0],\
##                                   [0.0,0.0,1.0]],self.Rotation)

        # GeomUtils' ellipse is oriented in a more conventional manor than mine
        # i.e., x = X, y = Y, z = Z...
        rotation=[self.Rotation[1][0],self.Rotation[1][1], \
                  self.Rotation[1][2], \
                  self.Rotation[2][0],self.Rotation[2][1], \
                  self.Rotation[2][2], \
                  self.Rotation[0][0],self.Rotation[0][1], \
                  self.Rotation[0][2]]
        cent = self.Trans
##        print rotation
##        print cent.x(),',',cent.y(),',',cent.z()
##        print a,b
##        xxx
        # use GeomUtils to find the intersections of the ellipse with the
        # surfaec mesh object CMesh... 
        Thetas = GeomUtils.EllipseCMeshIntersections(CMesh,cent,rotation,a,b)
##        print 'Thetas:', Thetas

        # loop through theta and compute phi (JD uses different parametrization
        # than i do)... 
        for theta in Thetas:
            if theta < math.pi/2.0:
                PhiList += [((math.atan(((a/b)*math.tan(theta))))/math.pi)-1.0]
            elif theta > math.pi/2.0 and theta < 3.0 * math.pi/2.0:
                PhiList += [ ((math.atan(((a/b)*math.tan(theta))))/math.pi) ]
            elif theta > 3.0*math.pi/2.0:
                PhiList += [((math.atan(((a/b)*math.tan(theta))))/math.pi)+1.0]
            elif theta == math.pi/2.0:
                PhiList += [-0.5]
            else: PhiList += [0.5] # i.e. theta == 3.0*math.pi/2.0
##            elif theta == 3.0*math.pi/2.0:
##                PhiList += [0.5]
##            else: raise 
        PhiList.sort()
##        print 'PhiList:', PhiList

        # Check if phi is a transition into (check = 1) or out of the mesh...
        # changed increase of phi from: phi+0.01 to: phi+0.1 - 5.20.05 because
        # Meshtools was returning IsPointOutsideMesh=true when it should not
        # have
        for phi in PhiList:
            y,z=self.__yzfunc(a,b,1.0,phi+0.1)
            u=Numeric.dot(self.TRotation,[0.0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            (flag,dist)=self.PointOutside(qpnt,self.model)

            if flag == 1:
                CheckList+=[0]
            else:
                CheckList+=[1]

        # CheckList must be populated... so if PhiList isn't generated a sample
        # of CheckList at phi = 0.0, & 1.0.  This will be used in WillGrow to 
        # guess if the ellipse is entirely inside or possibly entirely outside 
        # the volume (like could happen at a corner). 
        if len(PhiList) == 0:
            y,z=self.__yzfunc(a,b,1.0,1.0)
            u=Numeric.dot(self.TRotation,[0.0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            # note, dialed-up the tolerance in FemModel.elements.PointInside!
            (flag1,dist1)=self.PointOutside(qpnt,self.model)
            y,z=self.__yzfunc(a,b,1.0,0.0)
            u=Numeric.dot(self.TRotation,[0.0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            # note, dialed-up the tolerance in FemModel.elements.PointInside!
            (flag2,dist2)=self.PointOutside(qpnt,self.model)
            if flag1 == 1 and flag2 == 1:
                CheckList+=[0]
            else:
                CheckList+=[1]

        # when philist has 1 or 3 values it is nonsense.  raise an exception to
        # be handled in willgrow buy using my old method of finding phi's.
        if len(PhiList) == 1 or len(PhiList) == 3:
            raise DamErrors.BuildPhiError, (PhiList)

        return PhiList,CheckList 

######## Fellipse

    def __ChangeRotationToSecondPrincipal(self):
        # query database for princ stress
        (prince,evect) = self.model.GetPtStress(self.center).PrincipalValues()
        # if already permuted, use sig1
        if self.permute: 
            # make sure right handed rule...
            e1=Vec3D.Vec3D(evect[0][0],evect[0][1],evect[0][2]).Normalize()
            e2=Vec3D.Vec3D(evect[1][0],evect[1][1],evect[1][2]).Normalize()
            e3=Vec3D.CrossProd(e1,e2).Normalize()
            self.Rotation=Numeric.array([[e1.x(),e1.y(),e1.z()], \
                                         [e2.x(),e2.y(),e2.z()], \
                                         [e3.x(),e3.y(),e3.z()]])
        # otherwise permute now
        else:
            self.permute=True
            # make sure right handed rule...
            e1=Vec3D.Vec3D(evect[1][0],evect[1][1],evect[1][2]).Normalize()
            e2=Vec3D.Vec3D(evect[2][0],evect[2][1],evect[2][2]).Normalize()
            e3=Vec3D.CrossProd(e1,e2).Normalize()
            self.Rotation=Numeric.array([[e1.x(),e1.y(),e1.z()], \
                                         [e2.x(),e2.y(),e2.z()], \
                                         [e3.x(),e3.y(),e3.z()]])

######## Fellipse

    def __CheckPnts(self,inc):
        '''
        Returns PhiList and CheckList.  PhiList is a list of parametric angle,
        phi, where the ellipse intersects a surface.  CheckList is a list of
        booleans 0, or 1.  0 means, as you travel CCW around the ellipse you're
        leaving the volume (or mesh), 1 means you're entering it.
        '''

        # begin the search with the increment handed to __CheckPnts
        for i in xrange(3):
##            print "i:", i
            PhiList=[]; CheckList=[]
            phi_L,phi_H,flag = self.__InOrOut(self.model,-1.01,inc)
##            print phi_L,phi_H
            while phi_H < 1.0 and phi_H != None:
                dummy=0.5*(phi_L+phi_H)
                PhiList+=[0.5*(phi_L+phi_H)]
                CheckList+=[flag]
                phi_L,phi_H,flag = self.__InOrOut(self.model,phi_H,inc)

            # if we found some points, call off the hounds!  
            if len(PhiList)!=0:
                break
            # if we have not found some points assume the search needs to be
            # refined... this only happens 3 times until we give up altogether.
            inc=inc/1.5

        # return the last flag even if no phi found.  This will be used in
        # WillGrow to guess if the ellipse is entirely inside or possibly
        # entirely outside the volume (like could happen at a corner).  
        if len(PhiList)==0:
            CheckList+=[flag]

        # we have intentionally overlapped phi = -1.0 and phi = 1.0 a little.
        # correct dublicate points here if necessary.
        if len(PhiList) >= 2:
            one=(PhiList[0]+2.0)/PhiList[-1]
            if abs(one-1.0)<= 2.0*inc:
                del PhiList[-1]
                del CheckList[-1]

        return PhiList,CheckList

######## Fellipse

    def __FitPolynomial(self,aa,bb,center): 
        '''
        find the biquadratic curve that best fits the stress field normal to
        the crack surface.
        '''
        # initialize least square problem:
        #           XY * coef = f
        # XY is rectangular
        ypnts=[-aa,0.0,aa]
        zpnts=[-bb,0.0,bb]

##        inc=100
##        ypnts=[aa*(float(i)*(2.0/float(inc))-1.0) for i in range(inc+1)]
##        zpnts=[bb*(float(i)*(2.0/float(inc))-1.0) for i in range(inc+1)]
##        print ypnts
##        print zpnts
        f=[]
        XY=[]
        A = self.Rotation[0][0]
        B = self.Rotation[0][1]
        C = self.Rotation[0][2]
        # get stress at each sample point

        for y in ypnts: 
            for z in zpnts:
                # u = Numeric.dot(self.TRotation,[0,y,z])
                u0=self.Rotation[1][0]*y+self.Rotation[2][0]*z
                u1=self.Rotation[1][1]*y+self.Rotation[2][1]*z
                u2=self.Rotation[1][2]*y+self.Rotation[2][2]*z
                # pnt = center + Vec3D.Vec3D(u[0],u[1],u[2])
                qpnt=center+Vec3D.Vec3D(u0,u1,u2)
                sig=self.model.GetPtStress(qpnt)

                # GetPtStress returns a ColTensor, arrange to use Numeric the 
                # global stress tensor is...
                sigG=[[sig.xx(), sig.xy(), sig.zx()],
                      [sig.xy(), sig.yy(), sig.yz()],
                      [sig.zx(), sig.yz(), sig.zz()]]

                # sigL = R*sigG*R'
                #R_sigG=Numeric.dot(self.Rotation,sigG)
                #sigL=Numeric.dot(R_sigG,self.TRotation)

                # fill matrix and right-hand-side vector
                #f+=[sigL[0][0]]
                XY+=[[y, z, y*y, y*z, z*z, 1.0]]
                f+=[A*A*sig.xx() + B*B*sig.yy() + C*C*sig.zz() \
                   +2.0*A*B*sig.xy()+2.0*A*C*sig.zx()+2.0*B*C*sig.yz()]

        # add the center of ellipse...
        sig=self.model.GetPtStress(center)

        sigG=[[sig.xx(), sig.xy(), sig.zx()],
              [sig.xy(), sig.yy(), sig.yz()],
              [sig.zx(), sig.yz(), sig.zz()]]
        #R_sigG=Numeric.dot(self.Rotation,sigG)
        #sigL=Numeric.dot(R_sigG,self.TRotation)
        #f+=[sigL[0][0]]
        f+=[A*A*sig.xx() + B*B*sig.yy() + C*C*sig.zz() \
            +2.0*A*B*sig.xy()+2.0*A*C*sig.zx()+2.0*B*C*sig.yz()]
        y0,z0=0.0,0.0
        XY+=[[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]

        # solve linear least squares problem (method of normal equations)
        XYTXY=Numeric.dot(Numeric.transpose(XY),XY)
        XYTf=Numeric.dot(Numeric.transpose(XY),f)
        coef=LinearAlgebra.solve_linear_equations(XYTXY,XYTf)

        return coef

######## Fellipse

    def __FourPoints(self,PhiList,CheckList,tol): 
        '''
        Given PhiList and CheckList, determine four points on 2 lines that make
        up the surface intesecting the ellipse.  This assumes that there are
        only two lines that make up the surface.  
        '''

        a,b=self.GiveCurrent()

        if len(PhiList)==2: # either a "corner" of the volume is within
                            # ellipse or a straight line is crossing it

            y1a,z1a=self.__yzfunc(a,b,1.0,PhiList[0])
            y2a,z2a=self.__yzfunc(a,b,1.0,PhiList[1])

            # to increment Phi in __BisectRho by some approriate amount try
            PhiInc=0.33*(PhiList[1]-PhiList[0])

            # if this damage element is associated with a damage origin that is
            # a surface node, assume that the surface intersections occur at 
            # the node. Consequently,  __BisectRho by incrementing phi does not
            # make sense and we can return (0.0,0.0) for points 1b & 2b.
            # (jme 6/8/05)
            if self.is_surf:
                y1b,z1b,y2b,z2b = 0.0,0.0,0.0,0.0
                order=[0,1]

            elif CheckList[0]==0: # first point is leaving the volume
                
                rho_L,rho_H,junk= \
                        self.__BisectRho(-1.0,1.0,PhiList[0]+PhiInc,tol) 
                rho=0.5*(rho_L+rho_H)
                y1b,z1b=self.__yzfunc(a,b,rho,PhiList[0]+0.1)
                
                rho_L,rho_H,junk= \
                        self.__BisectRho(-1.0,1.0,PhiList[1]-PhiInc,tol) 
                rho=0.5*(rho_L+rho_H)
                y2b,z2b=self.__yzfunc(a,b,rho,PhiList[1]-0.1)
                # for len(PhiList) == 2, order is meaningless
                order=[0, 1]

            else: # first point is entering volume
                rho_L,rho_H,junk= \
                        self.__BisectRho(-1.0,1.0,PhiList[0]-PhiInc,tol)
                rho=0.5*(rho_L+rho_H)
                y1b,z1b=self.__yzfunc(a,b,rho,PhiList[0]-0.1)

                rho_L,rho_H,junk= \
                        self.__BisectRho(-1.0,1.0,PhiList[1]+PhiInc,tol)
                rho=0.5*(rho_L+rho_H)
                y2b,z2b=self.__yzfunc(a,b,rho,PhiList[1]+0.1)
                # for len(PhiList) == 2, order is meaningless
                order=[0, 1]

        elif len(PhiList)==4:
            # check the mid-points between PhiList[0] and PhiList[1], and
            # PhiList[0] and PhiList[3] to determine which is a surface line.
            # First get coords of three points...
            y1a,z1a=self.__yzfunc(a,b,1.0,PhiList[0])
            y1b,z1b=self.__yzfunc(a,b,1.0,PhiList[1])
            y2b,z2b=self.__yzfunc(a,b,1.0,PhiList[3])

            # average them
            y_check_1=0.5*(y1a+y1b)
            z_check_1=0.5*(z1a+z1b)
            y_check_2=0.5*(y1a+y2b)
            z_check_2=0.5*(z1a+z2b)

            # check IsPointIn() for each
            u=Numeric.dot(self.TRotation, \
                          [0.0,y_check_1,z_check_1])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            (flag,dist1)=self.PointOutside(qpnt,self.model)
##            if flag == 1:
##                check1 = 0
##            else:
##                check1 = 1
            u=Numeric.dot(self.TRotation, \
                          [0.0,y_check_2,z_check_2])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            (flag,dist2)=self.PointOutside(qpnt,self.model)
##            if flag == 1:
##                check2 = 0
##            else:
##                check2 = 1

            # now logic statments...
            if CheckList[0] == 0:
                # if both check points are close to or within the model...
                if dist1 < tol and dist2 < tol:
                    y1a,z1a=self.__yzfunc(a,b,1.0,PhiList[0])
                    y1b,z1b=self.__yzfunc(a,b,1.0,PhiList[1])
                    y2a,z2a=self.__yzfunc(a,b,1.0,PhiList[2])
                    y2b,z2b=self.__yzfunc(a,b,1.0,PhiList[3])
                    order=[0,1,2,3]
                else:
                    y1a,z1a=self.__yzfunc(a,b,1.0,PhiList[1])
                    y1b,z1b=self.__yzfunc(a,b,1.0,PhiList[2])
                    y2a,z2a=self.__yzfunc(a,b,1.0,PhiList[0])
                    y2b,z2b=self.__yzfunc(a,b,1.0,PhiList[3])
                    order=[1,2,0,3]
            else:
                # if both check points are close to or within the model...
                if dist1 < tol and dist2 < tol:
                    y1a,z1a=self.__yzfunc(a,b,1.0,PhiList[1])
                    y1b,z1b=self.__yzfunc(a,b,1.0,PhiList[2])
                    y2a,z2a=self.__yzfunc(a,b,1.0,PhiList[0])
                    y2b,z2b=self.__yzfunc(a,b,1.0,PhiList[3])
                    order=[1,2,0,3]
                else:
                    y1a,z1a=self.__yzfunc(a,b,1.0,PhiList[0])
                    y1b,z1b=self.__yzfunc(a,b,1.0,PhiList[1])
                    y2a,z2a=self.__yzfunc(a,b,1.0,PhiList[2])
                    y2b,z2b=self.__yzfunc(a,b,1.0,PhiList[3])
                    order=[0,1,2,3]

        else:
            print ' '
            print ' exception print: PhiList:  ', PhiList
            print ' exception print: CheckList:', CheckList
            print ' '
            raise ''' * len(PhiList) = 3, or > 4. See: Fellipse.__FourPoints'''
        # note: this exception is raised if len(PhiList) == 1 too, although
        # this event should not happen!  

        return y1a,z1a,y1b,z1b,y2a,z2a,y2b,z2b,order

######## Fellipse

    def __FtoH(self,ya,za,yb,zb,phi_a,phi_b,norm): 
        '''
        This function transitions from Fellipse class to Hellipse class.
        Returns Hellipse object.
        '''

        # Get current info
        a,b=self.GiveCurrent()

        # establish the "prime" coord system:  
        # origin is at phi_b
        # y' runs FROM phi_b TO phi_a
        # z' by right hand rule
        h=math.sqrt((ya-yb)**2.+(za-zb)**2.)
        s=(za-zb)/h
        c=(ya-yb)/h
        R=[[c, s], [-s, c]] # rotation matrix
        trans=[[yb], [zb]] # translation vector

        # the origin of the unprime system must have a positive z' coordinate
        # in order for the semi-ellipse to encompass the same area.  Check
        # to see if this is true.
        check_origin=-1.0*Numeric.dot(R,[[yb], [zb]])
        switch=0 # flag to indicate if orientation of y' must be switched

        if check_origin[1]< 0.0:
            # if the unprime origin does not have a positive z' component,
            # switch and a and b so the prime origin is at phi_a (but phi_a is
            # now labeled phi_b...)
            switch=1
            ya,za,yb,zb,phi_a,phi_b=yb,zb,ya,za,phi_b,phi_a
            h=math.sqrt((ya-yb)**2.+(za-zb)**2.)
            s=(za-zb)/h
            c=(ya-yb)/h
            R=[[c, s], [-s, c]] # rotation matrix
            trans=[[yb], [zb]]

        # check 12 points along remaining elliptical crack front, convert to
        # prime coordinates and record max y', max z', and min y'.  Use these
        # to determine the new a and b for the semi-ellipse.
        max_y_prime=0.
        max_z_prime=0.
        min_y_prime=0.

        # some increments that don't need to be recomputed each time through 
        incc=(phi_b-phi_a)
        indd=(2.0-(phi_a-phi_b))

        for i in xrange(13):
            phi_check=phi_a+(incc*(float(i)/12.0))
            if switch == 1:
                # if switch, want to increment phi backwards.
                phi_check = phi_a+(indd*(float(i)/12.0))
                if phi_check > 1.0:
                    phi_check=phi_check-2.0

            # get coords at phi_check
            y,z=self.__yzfunc(a,b,1,phi_check)

            # calc. their coords in prime system (origin at phi_b)
            y_prime=Numeric.dot(R,[[y-yb],[z-zb]])

            if y_prime[0]>max_y_prime:
                max_y_prime = y_prime[0][0]
    
            if y_prime[0]<min_y_prime:
                min_y_prime = y_prime[0][0]

            if y_prime[1]>max_z_prime:
                max_z_prime = y_prime[1][0]

        a_new=0.5*(abs(max_y_prime)+abs(min_y_prime))
        b_new=abs(max_z_prime)

        # to prevent center of ellipse moving outside the volume don't move
        # center of ellipse if it is a suface node (remember this is FOR THE
        # TRANSITION FROM Fellipse TO Hellipse).
        if self.is_surf:
            y00,z00=0.0,0.0
        else:
            cent=[[min_y_prime+a_new],[0]]
            cent=[[yb],[zb]]+Numeric.dot(Numeric.transpose(R),cent)
            y00,z00=cent[0][0],cent[1][0]

        trans=Numeric.dot(self.TRotation,[0.,y00,z00])
        trans=Vec3D.Vec3D(trans[0],trans[1],trans[2])+self.Trans

        # note: pass Hellipse - self.Rotation instead of self.Evect to assure
        # that x' (direction normal to ellipse) is the same. 
        return Hellipse(self.doid,self.Trans,self.model,a_new,b_new,norm, \
                        self.material,self.Rotation,self.verbose, \
                        self.verification)

######## Fellipse

    def __FtoQ(self,ya,za,yb,zb,y00,z00,phi_a,phi_b,norm): 
        '''
        This function transitions from Fellipse class to Qellipse class.
        Returns Qellipse object.

        Global coordinate system = the coord system of the background FE model.
        Prime coord system = coord system of the Fellipse.
        Double-primed coord system = a system set up within this method to
        correspond to the new Qellipse.  Maybe changed in the Qellipse
        constructor.

        (ya,za),(yb,zb) = coords, in primed system, of phi_a & phi_b
        (y00,z00) = coords of the vertex of the new Qellipse.  May need to be
        recomputed here!!  
        phi_a,phi_b = parametric angle where Fellipse intersects the FE model
        norm = norm computed in Intersect() in Global coords

        (Last modified 12/30/05 - works well for SIPS 3002, doid 119032)
        '''

        # Get current info
        a,b=self.GiveCurrent()

        # first check to make sure (ya,za) and (yb,zb) will result in a right
        # handed rule in the order they were passed in.  If not switch them:
        A=Vec3D.Vec3D(ya-y00, za-z00, 0.)
        B=Vec3D.Vec3D(yb-y00, zb-z00, 0.)
        R=Vec3D.CrossProd(A,B)
        switch=0 # flag to indicate if orientation of y' must be switched
        if R.z() < 0.:
            switch = 1
            ya,za,yb,zb=yb,zb,ya,za

        # establish the "double-prime" coord system (primed being Fellipse
        # space, unprimed being global coordinates):  
        # origin is at (y00, z00)
        # y'' runs FROM (y00, z00) TO phi_a
        # z'' by right hand rule (general direction of (y00, z00) TO phi_b)
        h=math.sqrt((ya-y00)**2.+(za-z00)**2.)
        s=(za-z00)/h
        c=(ya-y00)/h
        R=[[c, s], [-s, c]] # rotation matrix
        trans=[[y00], [z00]] # translation vector

        # check 12 points along remaining elliptical crack front, convert to
        # double-prime coordinates and record max y'', max z''.  Use these
        # to determine the new a and b for the quarter-ellipse.
        max_y_prime=0.; max_z_prime=0.
        min_y_prime=0.; min_z_prime=0.;

        for i in xrange(13):
            phi_check=phi_a+(phi_b-phi_a)*(float(i)/12.0)
            if switch == 1:
                phi_check = phi_b+(2.0-(phi_b-phi_a))*(float(i)/12.0)
                if phi_check > 1.0:
                    phi_check=phi_check-2.0

            # get coords at phi_check
            y,z=self.__yzfunc(a,b,1,phi_check)

            # calc. their coords in double-prime system (origin at (y00,z00))
            y_prime=Numeric.dot(R,[[y-y00],[z-z00]])

            if y_prime[0][0]>max_y_prime:
                max_y_prime = y_prime[0][0]

            if y_prime[1][0]>max_z_prime:
                max_z_prime = y_prime[1][0]

            # may need to move vertex of Qellipse
            if y_prime[0][0]<min_y_prime:
                min_y_prime = y_prime[0][0]

            # can't think of why min_z_prime would ever be less than 0.0
            if y_prime[1][0]<min_z_prime:
                min_z_prime = y_prime[1][0]

        # if min_y_prime < 0.0 move center & if is_surf & if the point is
        # inside the mesh, move the center of the qellipse and adjust
        # max_y_prime to account for the negative y_prime chunk.
        if min_y_prime < 0.0 and self.is_surf:
            # get the new center's coords in primed coords
            u=Numeric.dot(Numeric.transpose(R),[min_y_prime,0.0])

            # get the new center's coords in global coords
            qpnt=Numeric.dot(self.TRotation,[0.0, u[0],u[1]])
            qpnt=self.Trans+Vec3D.Vec3D(qpnt[0],qpnt[1],qpnt[2])
            (flag,dist)=self.PointOutside(qpnt,self.model)

            if not flag or dist < min(a,b)*1e-2:
                # move center
                y00,z00=u[0],u[1]

                # adjust max_y_prime (recall, min_y_prime < 0.0)
                max_y_prime = max_y_prime - min_y_prime

                # recompute norm to account for moved ellipse vertex
                new_norm = self.NewNorm(ya,za,yb,zb,y00,z00,self.TRotation)

                # we expect new_norm to be close to norm.  check that
                dot = new_norm*norm
                if dot < 0.5:
                    # for now, just print a message and use the old norm!  
                    print "new_norm, norm:", new_norm, norm
                else:
                    del(norm)
                    norm = new_norm
                    del(new_norm)

        trans=Numeric.dot(self.TRotation,[0.,y00,z00])
        trans=Vec3D.Vec3D(trans[0],trans[1],trans[2])+self.Trans

        if max_y_prime <= 0.0:
            raise 'dumm'
        if max_z_prime <= 0.0:
            raise 'dumm'

        # pass self.Rotation instead of Evect to reflect __Rotation().  
        return Qellipse(self.doid,self.Trans,self.model,max_y_prime, \
                        max_z_prime,norm,self.material,self.Rotation, \
                        self.verbose,self.verification)

######## Fellipse

    def __GetOnce(self,a,b): 
        '''
        GetOnce(a,b): return values that are independent of the stress field
        on the ellipse.
        '''

        # take abs() because if significant digits lead to ko2 < 0 (albeit
        # really small) the elliptical integral below is zero inducing a divide
        # by zero in __CalcKi
        if a >= b:
            ko2=abs((1.0-(b*b/a/a)))
        else:
            ko2=abs((1.0-(a*a/b/b)))

        E,K = self.ElliptInt(ko2)

        # commented because no scipy for python 2.4 
##        # found scipy calcs the first and second elliptical integrals, and i
##        # assume its faster than what i do above!
##        E,K = scipy.special.ellipe(ko2),scipy.special.ellipk(ko2)

        if ko2==0:
            E1=0.0
            E2=0.0

        else:
            E1=(1/ko2/3)*((1+b*b/a/a)*E-(2*b*b/a/a*K))
            E2=(1/ko2/5)*(-4*ko2*E1+ko2*E)

        return(ko2,E,K,E1,E2)

######## Fellipse

##    def __GiveLocal(self,a,b,na,nb): 
##        return (a+na)/2,(b+nb)/2

######## Fellipse

    def __InOrOut(self,model,phi_in,inc): 
        '''
        Start at phi_in (parametric coord for ellipse) and increase phi
        until the status of the ellipse (according to FemModel.IsPointIn())
        changes.  This effectively tells you that you have moved into or out of
        the FEM mesh.

        called in __CheckPnts()
        '''

        # get current ellipse info and update self.Trans
        a,b=self.GiveCurrent()

        # check the phi_in first
        y,z=self.__yzfunc(a,b,1.0,phi_in)
        u=Numeric.dot(self.TRotation,[0.0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        # FemModel.IsPointIn() returns 1 if point is inside the mesh and 0 if
        # the point is NOT in the mesh.
        (flago,dist)=self.PointOutside(qpnt,self.model)
        if flago == 1:
            flag = 0
        else:
            flag = 1

        switch=flag

        # the while loop actually only forms a bounding interval for phi in
        # which the point of interests lies.  __BisectPhi finds the point to 
        # the gnats ass.
        while switch==flag and phi_in < 1.001:
            flag=switch
            phi_in+=inc # incrementally increase phi_in
            y,z=self.__yzfunc(a,b,1,phi_in)
            u=Numeric.dot(self.TRotation,[0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            (flago,dist)=self.PointOutside(qpnt,self.model)
            if flago == 1:
                switch = 0
            else:
                switch = 1

        if phi_in >= 1.0:
            phi_L=None
            phi_H=None

        else:
            tol=0.001 # arbitrarily chosen for now 1/13/04
            phi_L,phi_H,flag=self.__BisectPhi(phi_in-inc,phi_in,1.0,tol)

        return phi_L,phi_H,flag

######## Fellipse

##    def Intersect(self,y1a,z1a,y1b,z1b,y2a,z2a,y2b,z2b,model):
##    if something breaks here... see v1.3 for the old code...

######## Fellipse

    def __Rotation(self): 
        '''
        Sets self.Rotation by comparing the surf normal to the sigma1
        principal direction.  (rewritten 10/7/05). ** Note: only call this if
        detect node is a surface node.

        For Fellipse, the only thing we worry about is avoiding the situation
        where the ellipse lies in the plane of the surface, OR one of the
        surface elements.  
        
        surf_norm,s_dev,compare,COS = self.AverageSurfNorm()

        surf_norm = average surface normal at this surface node
        s_dev = is the maximum standard deviation of any component in surf_norm
        compare = true if any ONE of the surface elements has a surface norm 
                  that lies in the direction (+/- 5 degrees) of sig1
        COS = the cosine of the maximum angle between any ONE of the surface
              elements surface normal and surf_norm (COS kind of makes compare
              obsolete)
        '''

        # a tolerance for the angles (cos(+/-angle))
        # angletol=0.9962 # cos(5 degrees), 5 wasn't enough 
        # see tray_wall doid 44351
        # angletol=0.9659 # cos(15 degrees)
        # angletol=0.9397 # cos(20 degrees)
        angletol=0.86603 # cos(30) 7/10/06

        # First make a unit vector in the direction of the sigma1
        sig1=Vec3D.Vec3D(self.Rotation[0][0],self.Rotation[0][1], \
                         self.Rotation[0][2]).Normalize()
        sig2=Vec3D.Vec3D(self.Rotation[1][0],self.Rotation[1][1], \
                         self.Rotation[1][2]).Normalize()
        sig3=Vec3D.Vec3D(self.Rotation[2][0],self.Rotation[2][1], \
                         self.Rotation[2][2]).Normalize()

        # get the average surface norm, surf_norm, and the largest standard
        # deviation amongst its components, s_dev, and calc. angle between
        # surf_norm and the direction of sigma1
        # compare = 1 if one of the adjacent surface elements is
        # in the same direciton of sig1
        surf_norm,s_dev,compare,COS = \
                            self.AverageSurfNorm(self.model,self.doid,sig1,\
                                                 angletol)
        cos=sig1*surf_norm

        # 1.  if s_dev is small, the surface is mostly flat.  check if 
        #     surf_norm is within angletol of sig1.  If yes, permute the
        #     eig vectors.
        if s_dev<=0.05:
            if abs(cos) > angletol:
                x1=sig2
                x2=sig3
                x3=sig1
                self.permute=True 
            else:
                x1=sig1
                x2=sig2
                x3=sig3
        # 2.  try to avoid fellipse outside volume, like at a corner when sig1
        #     isn't in the direction of the average  normal.  if the sig1
        #     vector falls within theta of a surface normal and
        #     the average surface normal (+5 degrees) check if within +/- 5
        #     degrees of sig2.  If yes, punt.  If no, permute.
        #     Note: cos(a+b) = cosacosb - sinasinb...
        #     for small b... = cosacosb
        elif abs(cos) >= COS*angletol:
            x3=-1.0*surf_norm
            x2=Vec3D.CrossProd(x3,sig1).Normalize()
            x1=Vec3D.CrossProd(x2,x3).Normalize()
        # 3.  to avoid the fellipse lieing in the plane of the surface at a
        #     corner where s_dev is > 0.5.  
        elif compare:
            x3=-1.0*surf_norm
            x2=Vec3D.CrossProd(x3,sig1).Normalize()
            x1=Vec3D.CrossProd(x2,x3).Normalize()
        else: 
            x1=sig1
            x2=sig2
            x3=sig3

        self.Rotation=[[x1.x(), x1.y(), x1.z()],
                       [x2.x(), x2.y(), x2.z()],
                       [x3.x(), x3.y(), x3.z()]]

######## Fellipse

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

        RANGE=range(4)

        if Ki==4: return 4

        # store for use elsewhere... 
        for i in RANGE: self.a[i][1]+=[Ki[i]]

        #  commented out to speedup... logically shouldn't need to check this
        # here anyway?! jme 7/15/06
##        # check 4 corners of ellipse (note: 1 and -1 are same point)
##        a,b=self.GiveCurrent()
##        flist=[]
##        for i in xrange(4):
##            phi=float(i)/2.0 - 1.0
##            y,z=self.__yzfunc(a,b,1.0,phi)
##            u=Numeric.dot(self.TRotation,[0,y,z])
##            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
##            (flag,dist)=self.PointOutside(qpnt,self.model)
##            if flag == 1:
##                flist+=[0]
##            else:
##                flist+=[1]
##
##        # list contains the indicator of whether the 4 points are in model or
##        # not.  If they are all zero, assume no remaining ligament.  Set
##        # will_grow = 4.
##        if flist[0]==0 and flist[1]==0 and flist[2]==0 and flist[3]==0:
##            return 4

        Kic=[blah.Kic for blah in self.dadN]
        Reff=[self.dadN[i].CalcReff(int(N),R) for i in RANGE]
        Kth=[self.dadN[i].Calc_DKth(self.a[i][0][-1],Reff[i]) \
             for i in range(4)]
        Kmax=[Ki[i]/(1.0-R) for i in RANGE]

        # max K is negative = compression, will not grow.
        if max(Kmax) <= 0.:
            return -1

        # will not grow, even under cyclic loading if...
        elif Ki[0] <= Kth[0] and Ki[1] <= Kth[1] and \
             Ki[2] <= Kth[2] and Ki[3] <= Kth[3]: 
            return 0

        # as in unstable growth if...
        elif Kmax[0] >= Kic[0] or Kmax[1] >= Kic[1] \
             or Kmax[2] >= Kic[2] or Kmax[3] >= Kic[3]: 
            return 2

        else: return 1

##        # as in stable, fatigue growth... yes, i know this could just be an
##        # else at the end of the day... but i wanted to be explicit.
##        elif (Ki[0] > Kth[0] and Ki[0] < Kic[0]) or \
##             (Ki[1] > Kth[1] and Ki[1] < Kic[1]) or \
##             (Ki[2] > Kth[2] and Ki[2] < Kic[2]) or \
##             (Ki[3] > Kth[3] and Ki[3] < Kic[3]): 
##
##            # but if any Ki > Kic will grow unstable...
##            if Ki[0] >= Kic[0]: 
##                will_grow = 2
##            elif Ki[1] >= Kic[1]:
##                will_grow = 2
##            elif Ki[2] >= Kic[2]:
##                will_grow = 2
##            elif Ki[3] >= Kic[3]:
##                will_grow = 2
##            else:
##                will_grow = 1
##
##        # as in unstable growth if...
##        elif Ki[0] >= Kic[0] or Ki[1] >= Kic[1] \
##             or Ki[2] >= Kic[2] or Ki[3] >= Kic[3]: 
##            will_grow = 2
##
##        return will_grow,'junk'

######## Fellipse

    def __init__(self,doid,center,model,a,b,material,rotation=None, \
                 verbose=None,verification=False):
        self.a=[DamHistory.DamHistory(a),DamHistory.DamHistory(b), \
                DamHistory.DamHistory(a),DamHistory.DamHistory(b)]
        self.CrackFrontPoints=4
        self.doid = doid
        self.names = ['Fellipse','a','b','na','nb']
        self.verification = verification
        self.model = model
        self.verbose = verbose
        self.center = center # store this incase someday not a node location
        self.Trans = center # location of cntr, update as necessary. currently,
                            # original center is always retrievable by
                            # xyz,delxyz,sigxyz=model.GetNodeInfo(self.doid)
        self.dN=[] # list of time steps Numeric.sum(self.dN) = Life
        self.nextdN = 0.0

        # a boolean set to true in self.Rotation() if crack is turned normal to
        # second principal stress... see also
        # self.ChangeRotationToSecondPrincipal()
        self.permute=False

        # check if my damage origin is a surface node
        self.is_surf=model.IsSurfaceNode(doid)
        # set dadN & store material to pass to Hellipse or Qellipse
        self.material=material
        self.dadN=[dadN.dadN(material,0),dadN.dadN(material,0), \
                   dadN.dadN(material,0),dadN.dadN(material,0)]
        # initialize the plastic zone size
        for i in range(self.CrackFrontPoints): self.dadN[i].SetKol()

        if rotation:
            # make sure right handed rule...
            e1=Vec3D.Vec3D(rotation[0][0],rotation[0][1], \
                           rotation[0][2]).Normalize()
            e2=Vec3D.Vec3D(rotation[1][0],rotation[1][1], \
                           rotation[1][2]).Normalize()
            e3=Vec3D.CrossProd(e1,e2).Normalize()
            self.Rotation=Numeric.array([[e1.x(),e1.y(),e1.z()], \
                                         [e2.x(),e2.y(),e2.z()], \
                                         [e3.x(),e3.y(),e3.z()]])
        else:
            # query database for princ stress
            (prince,evect) = model.GetPtStress(center).PrincipalValues()
            # make sure right handed rule...
            e1=Vec3D.Vec3D(evect[0][0],evect[0][1],evect[0][2]).Normalize()
            e2=Vec3D.Vec3D(evect[1][0],evect[1][1],evect[1][2]).Normalize()
            e3=Vec3D.CrossProd(e1,e2).Normalize()
            self.Rotation=Numeric.array([[e1.x(),e1.y(),e1.z()], \
                                         [e2.x(),e2.y(),e2.z()], \
                                         [e3.x(),e3.y(),e3.z()]])

        # if not proceed with Rotation=Evect
        if self.is_surf:
            self.__Rotation()

##        self.Rotation=[[0.0, 1.0, 0.0],
##                       [1.0, 0.0, 0.0],
##                       [0.0, 0.0, -1.0]]

        # store transpose(self.Rotation) to reduce calls to Numeric.Transpose
        self.TRotation = Numeric.transpose(self.Rotation)

        if self.verification:
            self.verification=True
            print self
            print ' ',self.names[0], ' Rotation matrix is:'
            print repr(Numeric.array(self.Rotation))+"\n"

######## Fellipse

    def __yzfunc(self,a,b,rho,phi): 
        '''
        define mapping from rectanglular coords to the ellipse
        '''
        y=a*((rho+1.0)/2.0)*math.cos(math.pi*(phi+1.0))
        z=b*((rho+1.0)/2.0)*math.sin(math.pi*(phi+1.0))
        return y,z

######## Fellipse

    def dAdN(self,N,size,args):
        '''
        Used by GrowDam in DamMo.py
        N - current number of cycles of loading
        size - a list of ellipses dimensions: [a, b, -a, -b]
        args - tuple of supporting arguments, in this case just scale! 
        '''
        RANGE=range(4)
        # unpack args
        scale=args[0]
        growth = []
        abab= [max(0,size[i]) for i in RANGE]
        a,b = (abab[0]+abab[2])/2.0 ,(abab[1]+abab[3])/2.0
        Xell = (abab[0]-abab[2])/2.0
        Yell = (abab[1]-abab[3])/2.0
        trans=Numeric.dot(self.TRotation, \
                      [0.0,Xell,Yell])
        localcenter=self.Trans + Vec3D.Vec3D(trans[0],trans[1],trans[2])
        KIs = self.__CalculateKi(a,b,localcenter,scale)
        for i in RANGE:
            growth.append(self.dadN[i].Calc_dadN(KIs[i],size[i],int(N)))

        if max(growth) <= 0.:
            raise DamErrors.dAdNError, ('Fellipse',\
                                  'All growth rates <= 0.0',growth)
        return growth

########################################################################
########################################################################

class Hellipse(Damage):

#   geometry & numbering
#
#   DamHistory:
#                 Eg. {'aba'  : [[a, b, -a], [a, b, -a], ...]
#                      'cent': [Vec3D, Vec3D_new, ...]
#                      0     : [[Ki, Kii, Kiii], [Ki, Kii, Kiii], ...]
#                      90    : [[Ki, Kii, Kiii], [Ki, Kii, Kiii], ...]
#                      180   : [[Ki, Kii, Kiii], [Ki, Kii, Kiii], ...]
#                      dN    : [0, dN1, dN2...]}
#                 for a crack.  Where the first list of a,b,Ki is for t=0 in 
#                 the evolution.
#
#  the center of the ellipse starts at the damage origin and can move
#
#  sigma1 defines the local x-axis of the ellipse
#
#  local x-axis (major axis and a) is determined by a cross product of the
#  local z-axis with the surface normal at the damage origin
#
#  local y-axis (minor axis and b) is then obtained by e3 X e1
#
#  self.Rotation - stores the rotation matrix s.t.:
#                        v = transpose(R)*v'
#                             -OR-
#  qpnt = Numeric.dot((self.TRotation,v_prime) + self.DamOro
#

######## Hellipse

    def CalcArea(self,ab): 
        '''
        simple method to calculate the current area of the ellipse.
        '''
        a,b=0.5*(ab[0]+ab[2]),ab[1]

        return 0.5*math.pi*a*b

######## Hellipse

    def CalculateKi(self,scale):
        a,b=self.GiveCurrent()
        try:
            Ki=self.__CalculateKi(a,b,self.Trans,scale)
            return Ki

        except DamErrors.FitPolyError, message:
            # for cases where the crack grows well outside the model,
            # model.GetPtStress returns 'point not in any elements'.  This is
            # flagged in __FitPoly & the K's are set to -10 (integer). If this
            # happens, assume the crack has become unstable.  This is different
            # than the next elif where i check list because that checks if the
            # entire crack is outside the mesh.
            if self.verbose:
                print ' FitPoly dumped for', message[0], 'at doid', self.doid,\
                      'due to',message[1] 
                print ' **SIMULATION OK, means crack has out grown its welcome'
            return 4

######## Hellipse

    def ComputeFrontPoints(self,crack_step=-1):
        '''
        to compute Vec3D objects that describe the crack front in global coords
        for the state at which a,b,center and rotation describe the crack.

        crack_step indicates which crack front to return.  the default is the
        current step
        '''

        def CurrentStep():
            a,b=self.GiveCurrent()
            return self.ParentComputeFrontPoints(a,b,self.Trans,self.Rotation,\
                                             self.__yzfunc)

        # if the crakc_step possibly describes a crack state... 
        if crack_step >= 0:
            if crack_step < len(self.a[0][0]):
                a = 0.5*(self.a[0][0][crack_step] + \
                         self.a[2][0][crack_step])
                b = self.a[1][0][crack_step]
                return self.ParentComputeFrontPoints(a,b,self.Trans, \
                                                self.Rotation,self.__yzfunc)

            # if not, return the current state. 
            else: return CurrentStep()
        else: return CurrentStep()

######## Hellipse

    def __CalculateKi(self,c,a,center,scale=1.0): 
        '''
        a helper function for calcKi that will allow computing Kis 
        without modifying internal data (JD)
        '''

        #  the least squares fit in support of Raju and Newman is pretty much
        # irrelevant because the bending multiplier, H, is unity for a/t << 1

        try:
            coef=self.__FitPolynomial(c,a,center)
        except MeshTools.EmptySearchResult, message:
            raise DamErrors.FitPolyError, ('Hellipse',message)
        except LinearAlgebra.LinAlgError, message:
            raise DamErrors.FitPolyError, ('Hellipse',message)

        A=coef[0]
        B=coef[1]

        St=B+A*a
        if A >= 0.0:
            Sb = 0.0
        else:
            Sb = B - St
        KI = []

        if a <= c:
            Q=1.0+1.464*(a/c)**1.65
            M1=1.13-0.09*a/c
##            M2=-0.54+(0.89/(0.2+(a/c)))
##            M3=0.5-(1.0/(0.65+(a/c)))+14.0*((1-(a/c))**24)
            g0=1.1 # t=large! g(180)=g(0)
##            g0=1.0+(0.1+0.35*(a/t)*(a/t))
            g90=1.0 # sin(90)=1!
            f0=math.sqrt(a/c) # g(0)=g(180)
            f90=1.0
            fw=1.0 # t=large,b=large!
            KI.append((St+Sb)*math.sqrt((math.pi*a/Q))* \
                                  (M1*g0*f0))
            KI.append((St+Sb)*math.sqrt((math.pi*a/Q))* \
                                  (M1*g90*f90))
            KI.append((St+Sb)*math.sqrt((math.pi*a/Q))* \
                                   (M1*g0*f0))
        elif a > c:
            Q=1.0+1.464*(c/a)**1.65
            M1=(math.sqrt(c/a))*(1.0+0.04*c/a)
            g0=1.1 # t=large! g(0)=g(180)
            g90=1.0 # t=large!
            f0=1.0 # g(0)=g(180)
            f90=math.sqrt(c/a)
            fw=1.0 # t=large, b=large
            KI.append((St+Sb)*math.sqrt((math.pi*a/Q))* \
                                  (M1*g0*f0))
            KI.append((St+Sb)*math.sqrt((math.pi*a/Q))* \
                                  (M1*g90*f90))
            KI.append((St+Sb)*math.sqrt((math.pi*a/Q))* \
                                   (M1*g0*f0))
        return JVT.ScalarMult(KI,scale)

######## Hellipse

    def GiveCurrent(self): 
        '''
        Returns the most current a, b.
        '''
        # we want average, and they should all be positive!  
        a = (self.a[0][0][-1] + \
             self.a[2][0][-1])/2.0
        b =  self.a[1][0][-1]

        assert (a>0. and b>0.)

        return a,b

######## Hellipse

    def IntegrateLife(self,ai,af,N_max,Life,tol):
        '''
        Using simpson's rule, integrate for life between ai and af.  For now,
        returns shortest life from the four possible positions.

        ai - initial crack length
        af - final crack length (usually float, but can be str()='end')
        r - increment factor (ai+1 = r * ai)
        '''

        # check the length of the SIF history... 
        if len(self.a[1]) == 0: return 0.0

        N=N_max
        Xlist,Klist=self.Lists(ai,af,self.a[0],self.a[1],tol)
        NN=self.Simpson(Xlist,Klist,Life,self.dadN[0],tol)
        if NN < N and NN >= 0:
            N = NN
        elif NN == -1:
            N = -1

        Xlist,Klist=self.Lists(ai,af,self.b[0],self.b[1],tol)
        NN=self.Simpson(Xlist,Klist,Life,self.dadN[1],tol)
        if NN < N and NN >= 0:
            N = NN
        elif NN == -1:
            N = -1

        Xlist,Klist=self.Lists(ai,af,self.na[0],self.na[1],tol)
        NN=self.Simpson(Xlist,Klist,Life,self.dadN[2],tol)
        if NN < N and NN >= 0:
            N = NN
        elif NN == -1:
            N = -1

        return N

######## Hellipse

    def ToFile(self,path_name,did):  # change for Gerd
        '''
        Print Hellipse info to file.
        '''

        if len(self.DamHistory[0])>0:
            a_his = open(path_name+'''.a''', 'a')
            _a_his = open(path_name+'''.na''', 'a')
            b_his = open(path_name+'''.b''', 'a')
            dN = open(path_name+'''.dN''', 'a')
            for i in xrange(len(self.DamHistory[0])):
                a_his.write(str(did)+' '+str(self.DamHistory['aba'][i][0])+  \
                                     ' '+str(self.DamHistory[0][i])+"\n")
                _a_his.write(str(did)+' '+str(self.DamHistory['aba'][i][2])+  \
                                     ' '+str(self.DamHistory[180][i])+"\n")
                b_his.write(str(did)+' '+str(self.DamHistory['aba'][i][1])+  \
                                     ' '+str(self.DamHistory[90][i])+"\n")
                try: dN.write(str(did)+' '+str(self.DamHistory['dN'][i])+"\n")
                except IndexError: pass
        else:
            return


######## Hellipse

    def UpdateState(self,anew):

        # update DamHistory
        for i in range(3):
            self.a[i][0]+=[anew[i]]
        y0_new=(anew[0]-anew[2])/2.0
        # update self.Trans.  Note, adding new position of center,
        # (0,y0_new,z0_new), w.r.t. ORIGINAL damage origin,
        # self.DamHistory['cent'][0]...
        trans=Numeric.dot(self.TRotation, \
                          [0.0,y0_new,0.0])
        self.Trans=self.Trans + \
                   Vec3D.Vec3D(trans[0],trans[1],trans[2])

######## Hellipse

    def GeometryCheck(self,CMesh): 
        '''
        Determines if the Fellipse is the appropriate type of damage
        element, based on geometry.

        Returns: dam_elem, have_doubt
        
        dam_elem = 0 if no change in crack geometry is necessary
                   4 if crack has outgrown geometry (net fracture)
                   a new damage element if transition is required

        have_doubt = a parameter to indicate how much the geometry is fudged.
        For example, if the particular damage origin does not fit any type of
        damage exactly, i.e. is partway between a half ellipse and quarter
        ellipse, can use this to raise a flag saying "hey i'm not that
        confident".  I have not done that much with this yet.

        CMesh is a list an object used in GeomUtils...
        '''

        a,b=self.GiveCurrent()

        # __BuildPhi to replace __CheckPnts; 11/3/04 uses John D's
        # GeomUtils.pyd
        PhiList,CheckList = self.__BuildPhi(CMesh,a,b) #@
##        print PhiList, CheckList
##        print self.__CheckPnts(0.01)
##        xx

        # if the whole semi-ellipse is within the volume PhiList[0] = None,
        # use Hellipse, calc. Ki and set will_grow...
        if len(PhiList)==0: return 0,0

        # calculate the full arclength
        full=self.__ArcLength(-1.0,1.0) #@

        # relative to smallest dimension of ellipse
        tol=min(a,b)/100000.

        # %%%%%%%%%%%% Begin logic statements for shape change %%%%%%%%%%%%%% #
        #                                                                     #
        # if len(CheckList) == 1 --> the semi-ellipse started near a corner and
        # has grown into it.
        if len(CheckList) == 1:
            if CheckList[0] == 0: # z' points into the body
                arc=self.__ArcLength(-1.0,PhiList[0])

                # compare ratio of arc/full and make decisions
                # 5-11-04 change from 0.75 to 0.6 to avoid hybrid
                # 1-14-05 change from 0.6 to 0.9 to get Qellipse
                if arc/full <= 0.9:
                    y1a,z1a=self.__yzfunc(a,b,1.0,PhiList[0])
                    rho_L,rho_H,junk = \
                        self.__BisectRho(-1.0,1.0,PhiList[0]+0.1,tol) #@ 
                    rho=0.5*(rho_L+rho_H)
                    y1b,z1b = \
                        self.__yzfunc(a,b,rho,PhiList[0]+0.1)

                    # make dam_history
                    y0_new = 0.5*(y1a+y1b) # won't be > 0.0! 
                    z0_new = 0.0
                    # note: use (y1a,0.0) for (ya,za) so as not to over
                    # condition. i.e. in NewNorm: length(ya,za) ~ length(yb,zb)
                    norm = self.NewNorm(abs(y1a),0.0,y1a,z1a,y0_new,z0_new,\
                                        self.TRotation)
                    a_new=a-y0_new
                    b_new=b
                    trans=Numeric.dot(self.TRotation, \
                                      [0.0,y0_new,z0_new])
                    trans=self.Trans+Vec3D.Vec3D(trans[0],trans[1],trans[2])

                    dam_el = self.__HtoQ(a_new,b_new,trans,norm)
                    have_doubt=0

                else: return 0,0

            else: # i.e. if CheckList[0] == 1
                arc=self.__ArcLength(PhiList[0],1.0)

                # compare ratio of arc/full and make decisions
                # 5-11-04 change from 0.6 to 0.75 to avoid hybrid
                # 1-14-05 change from 0.75 to 0.9 to get Qellipse
                if arc/full <= 0.9:
                    # find the point at the intersection 
                    y1a,z1a=self.__yzfunc(a,b,1.0,PhiList[0])
                    # now back outside volume a little (phi-0.1) and bisect
                    # rho to find a point on the surface of the vol. inside the
                    # ellipse.  
                    rho_L,rho_H,junk = \
                        self.__BisectRho(-1.0,1.0,PhiList[0]-0.1,tol)
                    rho=0.5*(rho_L+rho_H)
                    y1b,z1b=self.__yzfunc(a,b,rho,PhiList[0]-0.1)

                    # make dam_history
                    y0_new = 0.5*(y1a+y1b) # won't be < 0.0! 
                    z0_new = 0.0
                    # note: use (y1a,0.0) for (ya,za) so as not to over
                    # condition. i.e. in NewNorm: length(ya,za) ~ length(yb,zb)
                    norm = self.NewNorm(y1a,z1a,-y1a,0.0,y0_new,z0_new,\
                                        self.TRotation)
                    a_new=b
                    b_new=a+y0_new
                    trans=Numeric.dot(self.TRotation, \
                                      [0.0,y0_new,z0_new])
                    trans=self.Trans+Vec3D.Vec3D(trans[0],trans[1],trans[2])
                    dam_el = self.__HtoQ(a_new,b_new,trans,norm)
                    have_doubt=0

                else: return 0,0

        elif len(CheckList) == 2:
            arc=self.__ArcLength(PhiList[0],PhiList[1])

            # not sure this can happen because z' for Hellipse should point
            # into the body.  
            if CheckList[0]==0: return 0,1

            else: # i.e. if CheckList[0] == 1
                # 5-11-04 change 0.6 to 0.85 to force Qellipse
                if arc/full <= 0.85:
                    y1a,z1a=self.__yzfunc(a,b,1.0,PhiList[0])
                    y1b,z1b=self.__yzfunc(a,b,1.0,PhiList[1])

                    # make Qellipse
                    a_new=math.sqrt((y1a)**2+(z1a)**2)
                    b_new=math.sqrt((y1b)**2+(z1b)**2)
                    trans=self.Trans

                    # compute the new norm
                    norm=self.NewNorm(y1a,z1a,y1b,z1b,0.0,0.0,self.TRotation) 

                    dam_el = self.__HtoQ(a_new,b_new,trans,norm)
                    have_doubt=0

                else: return 0,0

        else: return 0,1

        return dam_el,have_doubt

######## Hellipse

    def __ArcLength(self,phi_1,phi_2): 
        '''
        A nine point Gauss quadrature that calculates the arc length on the
        ellipse between phi_a and phi_b.
        '''

        # get current major and minor axis dimensions
        a,b=self.GiveCurrent()

        # define the arclength lambda
        arc=lambda x: math.sqrt(a*a*((math.sin((x+1)*math.pi/2.0))**2.0) + \
                                b*b*((math.cos((x+1)*math.pi/2.0))**2.0))

        # define gauss weights and gauss points
        gp9 = [
        -0.968160239,
        -0.836031107,
        -0.613371432,
        -0.324253423,
         0.0,
         0.324253423,
         0.613371432,
         0.836031107,
         0.968160239]

        gw9 = [
        0.081274388,
        0.180648160,
        0.260610696,
        0.312347077,
        0.330239355,
        0.312347077,
        0.260610696,
        0.180648160,
        0.081274388]

        # loop over gp9 and sum
        nsum=0

        for i in xrange(len(gp9)):
            ri = phi_1+(phi_2-phi_1)*(gp9[i]+1.0)/2.0
            wi = gw9[i]
            y = arc(ri) # arc is the lambda defined above
            nsum+=wi*y
        me=nsum*(math.pi/2.0)*(phi_2-phi_1)/2.0
        return me

######## Hellipse

    def __BisectPhi(self,phi_L,phi_H,rho,tol):
        '''
        Finds the intersection of the ellipse with the surface of the FEM
        mesh.
        '''

        a,b=self.GiveCurrent()

        # phi_L
        y,z=self.__yzfunc(a,b,rho,phi_L)
        u=Numeric.dot(self.TRotation,[0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        (flag,dist)=self.PointOutside(qpnt,self.model)
        if flag == 1:
            flag_L = 0
        else:
            flag_L = 1

        # phi_H
        y,z=self.__yzfunc(a,b,rho,phi_H)
        u=Numeric.dot(self.TRotation,[0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        (flag,dist)=self.PointOutside(qpnt,self.model)
        if flag == 1:
            flag_H = 0
        else:
            flag_H = 1

        while abs(phi_H-phi_L) > tol:
            phi_M=0.5*(phi_H+phi_L)
            y,z=self.__yzfunc(a,b,rho,phi_M)
            u=Numeric.dot(self.TRotation,[0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            (flag,dist)=self.PointOutside(qpnt,self.model)
            if flag == 1:
                flag_M = 0
            else:
                flag_M = 1

            if flag_M == flag_H:
                phi_H = phi_M
            else:
                phi_L = phi_M

        # flag_H indicates whether avg(phi_L,phi_H) is a transition from or
        # into the volume.  i.e., if flag_H = 0 it is a transition OUT of the
        # volume.  Conversly, if flag_H = 1 it is a transition INTO the volume.
        return phi_L, phi_H, flag_H

######## Hellipse

    def __BisectRho(self,rho_L,rho_H,phi,tol): 
        '''
        Finds the intersection of the ellipse with the surface of the FEM
        mesh.
       '''

        a,b=self.GiveCurrent()

        # phi_L
        y,z=self.__yzfunc(a,b,rho_L,phi)
        u=Numeric.dot(self.TRotation,[0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        (flag,dist)=self.PointOutside(qpnt,self.model)
        if flag == 1:
            flag_L = 0
        else:
            flag_L = 1

        # phi_H
        y,z=self.__yzfunc(a,b,rho_H,phi)
        u=Numeric.dot(self.TRotation,[0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        (flag,dist)=self.PointOutside(qpnt,self.model)
        if flag == 1:
            flag_H = 0
        else:
            flag_H = 1

        while abs(rho_H-rho_L) > tol:
            rho_M=0.5*(rho_H+rho_L)
            y,z=self.__yzfunc(a,b,rho_M,phi)
            u=Numeric.dot(self.TRotation,[0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            (flag,dist)=self.PointOutside(qpnt,self.model)
            if flag == 1:
                flag_M = 0
            else:
                flag_M = 1

            if flag_M == flag_H:
                rho_H = rho_M
            else:
                rho_L = rho_M

        # flag_H indicates whether avg(rho_L,rho_H) is a transition from or
        # into the volume.  i.e., if flag_H = 0 it is a transition OUT of the
        # volume.  Conversly, if flag_H = 1 it is a transition INTO the volume.
        return rho_L, rho_H, flag_H

######## Hellipse

    def __BuildPhi(self,CMesh,a,b): 
        '''
        use GeomUtils to build PhiList and CheckList.
        '''

        PhiList=[]; CheckList=[]
        # GeomUtils' ellipse is oriented in a more conventional manor than mine
        # i.e., x = X, y = Y, z = Z...
        rotation=[self.Rotation[1][0],self.Rotation[1][1], \
                  self.Rotation[1][2], \
                  self.Rotation[2][0],self.Rotation[2][1], \
                  self.Rotation[2][2], \
                  self.Rotation[0][0],self.Rotation[0][1], \
                  self.Rotation[0][2]]
        cent = self.Trans

        # use GeomUtils to find the intersections of the ellipse with the
        # surface mesh object CMesh...
        Thetas=GeomUtils.EllipseCMeshIntersections(CMesh,cent,rotation,a,b,18)

        # loop through theta and compute phi (JD uses different parametrization
        # than i do)... 
        for theta in Thetas:
            # now, if any theta is tol less than pi it is an Hellipse...
            # math.pi-0.0015 should be just inside body (if we find
            # intersection behond this it is not of interest)
            if theta >= math.pi-0.0015:
                continue
            # also, if theta is within tol of zero, we are not intersted
            # in the intersection...
            elif theta <= 0.0015:
                continue
            
            if a == b:
                e = 0.0
                r = a*((1-e*e)/(1-e*e*Numeric.cos(theta)* \
                                Numeric.cos(theta)))**.5
                phi=(2/Numeric.pi)*Numeric.arccos((r/a)* \
                                                   Numeric.cos(theta))-1.0
            elif a > b:
                e = (1-b*b/a/a)**.5
                r = a*((1-e*e)/(1-e*e*Numeric.cos(theta)* \
                                Numeric.cos(theta)))**.5
                phi=(2/Numeric.pi)*Numeric.arccos((r/a)* \
                                                   Numeric.cos(theta))-1.0
            else:
                a,b = b,a
                e = (1-b*b/a/a)**.5
                r = a*((1-e*e)/(1-e*e*Numeric.cos(theta)* \
                                Numeric.cos(theta)))**.5
                phi=(2/Numeric.pi)*Numeric.arccos((r/a)* \
                                            Numeric.cos(theta))-1.0
                a,b = b,a
            PhiList += [phi]
        PhiList.sort()

        # Check if phi is a transition into (CheckList[i] = 1) or out of the 
        # mesh... changed increase of phi from: phi+0.01 to: phi+0.1 - 5.20.05 
        # because Meshtools was returning IsPointOutsideMesh=true when it 
        # should not have.
        #
        # 6/5/05 - amount to increase phi for building CheckList should be
        # consistent with ranges of "angle" in Fellipse.WillGrow!  But that is 
        # hard to do... this is kind of hoakie and could use a better idea.
        #
        # OK here's a better idea.  The Hellipse should always be oriented such
        # that z' is into the body (by convention) so let's exploit that
        # knowledge.   For the case when phi is near 1.0, instead of moving
        # increment PLUS phi and asking if the point is in/out, move MINUS phi.
        # This check will eliminate problem that occurs when using Hellipse for
        # some slightly concave chunk.  And visa-versa for when phi is near
        # -1.0
        for phi in PhiList:
            y,z=self.__yzfunc(a,b,1.0,phi-0.01)
            u=Numeric.dot(self.TRotation,[0.0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            # note, dialed-up the tolerance in FemModel.element.PointInside!
            (flag,dist)=self.PointOutside(qpnt,self.model)
#            if phi > 0.0: 
            if flag == 1:
                CheckList+=[1] # REVERSED Because asking phi-0.01
            else:
                CheckList+=[0] 
##            else:  
##                if flag == 1:
##                    CheckList+=[0] 
##                else:
##                    CheckList+=[1] 

        return PhiList,CheckList

######## Hellipse

    def __CheckPnts(self,inc):
        '''
        Returns PhiList and CheckList.  PhiList is a list of parametric angle,
        phi, where the ellipse intersects a surface.  CheckList is a list of
        booleans 0, or 1.  0 means, as you travel CCW around the ellipse you're
        leaving the volume (or mesh), 1 means you're entering it.
        '''

        # begin the search with the increment handed to __CheckPnts
        for i in xrange(3):
            PhiList=[]; CheckList=[]
            phi_L,phi_H,flag = self.__InOrOut(self.model,-1.1,inc)

            while phi_H < 1.0 and phi_H != None:
                PhiList+=[0.5*(phi_L+phi_H)]
                CheckList+=[flag]
                phi_L,phi_H,flag = self.__InOrOut(self.model,phi_H,inc)

            # if we found some points, call off the hounds!  
            if len(PhiList)!=0:
                break
            # if we have not found some points assume the search needs to be
            # refined... this only happens 3 times until we give up altogether.
            inc=inc/1.5

        return PhiList,CheckList

######## Hellipse
    
    def __FitPolynomial(self,aa,bb,center):
        ypnts = [-aa,0.0,aa]
        zpnts = [0.0, bb/2.0, bb]
        # initialize least square problem:
        #           XY * coef = f
        # XY is rectangular
        f=[]
        XY=[]
        A = self.Rotation[0][0]
        B = self.Rotation[0][1]
        C = self.Rotation[0][2]

        # get stress at each guass point
        for y in ypnts:
            for z in zpnts:
                #u=Numeric.dot(self.TRotation,[0,y,z])
                #qpnt=center+Vec3D.Vec3D(u[0],u[1],u[2])
                u0=self.Rotation[1][0]*y+self.Rotation[2][0]*z
                u1=self.Rotation[1][1]*y+self.Rotation[2][1]*z
                u2=self.Rotation[1][2]*y+self.Rotation[2][2]*z
                qpnt=center+Vec3D.Vec3D(u0,u1,u2)
##                print aa,y,bb,z
##                print 'v',qpnt
                sig=self.model.GetPtStress(qpnt)
##                print 'out'
                # GetPtStress returns a ColTensor, arrange to use Numeric the 
                # global stress tensor is...
                sigG=[[sig.xx(), sig.xy(), sig.zx()],
                      [sig.xy(), sig.yy(), sig.yz()],
                      [sig.zx(), sig.yz(), sig.zz()]]

                # fill matrix and right-hand-side vector
                #f+=[sigL[0][0]], in place of sigL = R*sigG*R'...
                XY+=[[z, 1]]
                f+=[A*A*sig.xx() + B*B*sig.yy() + C*C*sig.zz() \
                    +2.0*A*B*sig.xy()+2.0*A*C*sig.zx()+2.0*B*C*sig.yz()]

        # add the center of ellipse...
##        print center
        sig=self.model.GetPtStress(center)
        sigG=[[sig.xx(), sig.xy(), sig.zx()],
              [sig.xy(), sig.yy(), sig.yz()],
              [sig.zx(), sig.yz(), sig.zz()]]
        #R_sigG=Numeric.dot(self.Rotation,sigG)
        #sigL=Numeric.dot(R_sigG,self.TRotation)
        #f+=[sigL[0][0]]
        f+=[A*A*sig.xx() + B*B*sig.yy() + C*C*sig.zz() \
            +2.0*A*B*sig.xy()+2.0*A*C*sig.zx()+2.0*B*C*sig.yz()]

        XY+=[[0.0, 1]]

        # Solve the least squares problem
        XYTXY=Numeric.dot(Numeric.transpose(XY),XY)
        XYTf=Numeric.dot(Numeric.transpose(XY),f)
        coef=LinearAlgebra.solve_linear_equations(XYTXY,XYTf)
        return coef

######## Hellipse

    def __HtoQ(self,a_new,b_new,trans,norm):
        '''
        method to build a Qellipse oject from given arguments.
        '''

        # return new damage element
        # note: pass Hellipse - self.Rotation instead of self.Evect
        # to assure that x' (direction normal to ellipse) is the
        # same.
        return Qellipse(self.doid,self.Trans,self.model,a_new, \
                        b_new,norm,self.material,self.Rotation,self.verbose, \
                        self.verification)

######## Hellipse

    def __InOrOut(self,phi_in,inc): #@
        '''
        Start at phi_in (parametric coord for ellipse) and increase phi
        until the status of the ellipse (according to FemModel.IsPointIn())
        changes.  This effectively tells you that you have moved into or out of
        the FEM mesh.
        '''

        # get current ellipse info
        a,b=self.GiveCurrent()

        # check the phi_in first
        y,z=self.__yzfunc(a,b,1.0,phi_in)
        u=Numeric.dot(self.TRotation,[0.0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        
        # FemModel.IsPointIn() returns 1 if point is inside the mesh and 0 if
        # the point is NOT in the mesh.
        (flago,dist)=self.PointOutside(qpnt,self.model)
        if flago == 1:
            flag = 0
        else:
            flag = 1

        switch=flag

        # the while loop actually only forms a bounding interval for phi in
        # which the point of interests lies.  __BisectPhi finds the point to 
        # the gnats ass.
        while switch==flag and phi_in < 1.001:
            flag=switch
            phi_in+=inc # incrementally increase phi_in
            y,z=self.__yzfunc(a,b,1.0,phi_in)
            u=Numeric.dot(self.TRotation,[0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            (flag,dist)=self.PointOutside(qpnt,self.model)
            if flag == 1:
                switch = 0
            else:
                switch = 1

        if phi_in >= 1.0:
            phi_L=None
            phi_H=None
            flag=None
        else:
            tol=0.001 # arbitrarily chosen for now 1/13/04
            phi_L,phi_H,flag=self.__BisectPhi(phi_in-inc,phi_in,1.0,tol)

        return phi_L,phi_H,flag

######## Hellipse

    def WillGrow(self,N,R,Ki): 
        '''
        Calculates the appropriate damage parameter,
        in this case stress intensity factor, and sets will_grow = :

        will_grow = -1 --> stress field is compressive (ie. Ki < 0.0)
        will_grow = 0  --> Ki insufficient to drive crack (Ki < K threshold)
        will_grow = 1  --> Stable, fatigue growth
        will_grow = 2  --> Unstable growth
        will_grow = 3  --> Change damage element type
        will_grow = 4  --> Assume net fracture (crack grew outside body)

        '''

        RANGE=range(3)

        if Ki == 4: return 4

        # store for use elsewhere... 
        for i in RANGE: self.a[i][1]+=[Ki[i]]

        #  commented out to speedup... logically shouldn't need to check this
        # here anyway?! jme 7/15/06
##        a=self.a[0][0][-1]
##        b=self.a[1][0][-1]
##        neg_a=self.a[2][0][-1]
##
##        # check 5 points around of semi-ellipse 
##        flist=[]
##        for i in xrange(5):
##            phi=float(i)/2.0 - 1.0
##            y,z=self.__yzfunc(a,b,1.0,phi)
##            u=Numeric.dot(self.TRotation,[0,y,z])
##            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
##            (flag,dist)=self.PointOutside(qpnt,self.model)
##            if flag == 1:
##                flist+=[0]
##            else:
##                flist+=[1]
##
##        # list contains the indicator of whether the 5 points are in model or
##        # not.  If they are all zero, assume no remaining ligament.  Set
##        # will_grow = 2.
##        if flist[0]==0 and flist[1]==0 and flist[2]==0 and \
##           flist[3]==0 and flist[4]==0:
##            return 4

        Kic=[blah.Kic for blah in self.dadN]
        Reff=[self.dadN[i].CalcReff(int(N),R) for i in RANGE]
        Kth=[self.dadN[i].Calc_DKth(self.a[i][0][-1],Reff[i]) \
             for i in RANGE]

        Kmax=[Ki[i]/(1.0-R) for i in RANGE]

        # max K is negative = compression, will not grow.
        if max(Ki) <= 0: 
            return -1 

        # will not grow, even under cyclic loading if...
        elif Ki[0] <= Kth[0] and Ki[1] <= Kth[1] and Ki[2] <= Kth[2]: 
            return 0

        # as in unstable growth if...
        elif Kmax[0] >= Kic[0] or Kmax[1] >= Kic[1] or Kmax[2] >= Kic[2]: 
            return 2

        else: return 1

######## Hellipse

    def __init__(self,doid,center,model,a,b,norm,material,rotation=None, \
                 verbose=False,verification=False):
        self.a=[DamHistory.DamHistory(a),DamHistory.DamHistory(b), \
                DamHistory.DamHistory(a)]
        self.CrackFrontPoints=3
        self.doid=doid
        self.names=['Hellipse','a','b','na']
        self.verification = verification
        self.model = model
        self.verbose=verbose
        self.Trans = center # location of cntr, update as necessary. 
                          # original center is always retrievable by
                          # xyz,delxyz,sigxyz=model.GetNodeInfo(self.doid) 
        self.dN=[] # list of time steps Numeric.sum(self.dN) = Life
        self.nextdN = 0.0

        # initiate the dadN model & store material to pass to Qellipse
        self.material = material
        self.dadN=[dadN.dadN(material,1),dadN.dadN(material,0), \
                   dadN.dadN(material,1)]
        # initialize the plastic zone size
        for i in range(self.CrackFrontPoints): self.dadN[i].SetKol()

        # ***** --> changed this 1/6/05.  occurs to me that Hellipse is never
        # instanced without norm being passed in.  furthermore, norm is always
        # calculated to be the outward surface normal in the plane of the
        # ellipse AND self.Rotation is the rotation matrix of the 
        # parent damage origin (i.e. NOT necessarily the evects).  See __FtoH()
        # and __FtoH for details --->
        self.Norm = norm

        # x1 is in direction of first principal stress
        x1=(Vec3D.Vec3D(rotation[0][0], rotation[0][1], \
                        rotation[0][2])).Normalize()
        # x2 is normal to x1 and norm
        x2=(Vec3D.CrossProd(x1,self.Norm)).Normalize()
        # x3 by cross x1 X x2...
        x3=(Vec3D.CrossProd(x1,x2)).Normalize()
        # rotation matrix then is...
##        print norm.x(), norm.y(), norm.z()
##        print x1.x(), x1.y(), x1.z()
##        print x2.x(), x2.y(), x2.z()
##        print x3.x(), x3.y(), x3.z()
##        x
        self.Rotation=[[x1.x(), x1.y(), x1.z()],
                       [x2.x(), x2.y(), x2.z()],
                       [x3.x(), x3.y(), x3.z()]]

##        self.Rotation=[[0.0, 1.0, 0.0],
##                       [1.0, 0.0, 0.0],
##                       [0.0, 0.0, -1.0]]

        # store transpose(self.Rotation) to reduce calls to Numeric.Transpose
        self.TRotation = Numeric.transpose(self.Rotation)

        if verification:
            self.verification=True
            print self
            print ' ',self.names[0], ' Rotation matrix is:'
            print repr(Numeric.array(self.Rotation))+"\n"

######## Hellipse

    def __yzfunc(self,a,b,rho,phi): 
        '''
        define mapping from rectanglular coords to the semi-ellipse
        '''
        y=a*((rho+1.0)/2.0)*math.cos(math.pi*(phi+1.0)/2.0)
        z=b*((rho+1.0)/2.0)*math.sin(math.pi*(phi+1.0)/2.0)
        return y,z

######## Hellipse

    def dAdN(self,N,size,args):
        '''
        Used by GrowDam in DamMo.py
        N - current number of cycles of loading
        size - a list of ellipses dimensions: [a, b, -a, -b]
        args - tuple of supporting arguments, in this case just scale! 
        '''
        RANGE=range(3)
        # unpack args
##        print 'in dAdN', N, size
        scale=args[0]
        growth = []
        abab= [max(0,size[i]) for i in RANGE]
        a,b = (abab[0]+abab[2])/2.0 ,abab[1]
        Xell = (abab[0]-abab[2])/2.0
        trans=Numeric.dot(self.TRotation, \
                      [0.0,Xell,0.0])
        localcenter=self.Trans + Vec3D.Vec3D(trans[0],trans[1],trans[2])
##        print 'before'
        KIs = self.__CalculateKi(a,b,localcenter,scale)
##        print 'KIs',KIs
        for i in RANGE:
            growth.append(self.dadN[i].Calc_dadN(KIs[i],size[i],int(N)))

        if max(growth) <= 0.:
            raise DamErrors.dAdNError, ('Hellipse',\
                                  'All growth rates <= 0.0',growth)
        return growth

########################################################################
########################################################################

class Qellipse(Damage):

#   geometry & numbering
#
#   DamHistory:
#                 Eg. {'ab'  : [[a, b,], [a, b,], ...]
#                      'cent': [[0, 0], [y0_new, z0_new], ...]
#                      0     : [[Ki, Kii, Kiii], [Ki, Kii, Kiii], ...]
#                      90    : [[Ki, Kii, Kiii], [Ki, Kii, Kiii], ...]
#                      dN    : [0, dN1, dN2...]}
#                 for a crack.  Where the first list of a,b,Ki is for t=0 in 
#                 the evolution.
#
#  the center of the ellipse starts at the damage origin and can move
#
#  sigma1 defines the local z axis (e3) of the ellipse
#
#  local x axis (major axis and a) is determined by a cross product of the
#  local z-axis with the surface normal at the damage origin
#
#  local y axis (minor axis and b) is then obtained by e3 X e1

######## Qellipse

    def CalcArea(self,ab): 
        '''
        simple method to calculate the current area of the ellipse.
        '''
        a,b=ab[0],ab[1]

        return 0.25*math.pi*a*b

######## Qellipse

    def CalculateKi(self,scale):
        a,b=self.GiveCurrent()
        try:
            Ki=self.__CalculateKi(a,b,self.Trans,scale)
            return Ki 

        except DamErrors.FitPolyError, message:
            # for cases where the crack grows well outside the model,
            # model.GetPtStress returns 'point not in any elements'.  This is
            # flagged in __FitPoly and the K's are set to -10 (integer).  If 
            # happens, assume the crack has become unstable.  This is different
            # than the next elif where i check list because that checks if the
            # entire crack is outside the mesh.
            if self.verbose:
                print ' FitPoly dumped for', message[0], 'at doid', self.doid,\
                      'due to',message[1] 
                print '** SIMULATION OK, means crack has out grown its welcome'
            # remove the last stuff in DamHistory because we couldn't compute
            # corresponding K's.
            return 4

######## Qellipse

    def ComputeFrontPoints(self,crack_step=-1):
        '''
        to compute Vec3D objects that describe the crack front in global coords
        for the state at which a,b,center and rotation describe the crack.

        crack_step indicates which crack front to return.  the default is the
        current step
        '''

        def CurrentStep():
            a,b=self.GiveCurrent()
            return self.ParentComputeFrontPoints(a,b,self.Trans,self.Rotation,\
                                             self.__yzfunc)

        # if the crakc_step possibly describes a crack state... 
        if crack_step >= 0:
            if crack_step < len(self.a[0][0]):
                a = self.a[0][0][crack_step]
                b = self.a[1][0][crack_step]
                return self.ParentComputeFrontPoints(a,b,self.Trans, \
                                                self.Rotation,self.__yzfunc)

            # if not, return the current state. 
            else: return CurrentStep()
        else: return CurrentStep()

######## Qellipse

    def __CalculateKi(self,c,a,center,scale=1.0):
        # try __FitPoly.  If no, then points lie outside the model, assume Ki
        # exceeds Kic!  Set the last Ki at 0, 90, & 180 to -10.  Look for this
        # flag in __WillGrow

        try:
            coef=self.__FitPolynomial(c,a,center)
        except MeshTools.EmptySearchResult, message:
            raise DamErrors.FitPolyError, ('Qellipse',message)
        except LinearAlgebra.LinAlgError, message:
            raise DamErrors.FitPolyError, ('Qellipse',message)

        # sigma_normal = A*z + B (eqn. of a line)
        A=coef[0]
        B=coef[1]
        # pure tension component of the normal stress, St.  This is
        # conservative because in pure tension (i.e. no Sb) the slope, A, would
        # be zero.  Here, if it is positive we are calculating a larger tensile
        # component than there actually is.  
        St=B+A*a

        # find bending component of normal tension, Sb.  If A > 0 means bending
        # moment is actually closing the crack.  To be conservative, set Sb = 0
        # *** NOTE: Bending multiplier, Hc = 1.0, for a/t~0...
        if A >= 0.0:
            Sb = 0.0 
        else:
            Sb = B - St
        KI=[]
        if a <= c:
            Q=1.0+1.464*(a/c)**1.65
            M1=1.08-0.03*a/c
            g0=1.08 # t=large!
            g90=1.0 # t=large!
            f0=math.sqrt(a/c) 
            f90=1.0
            fw=1.0 # t=large, b=large
            KI.append((St+Sb)*math.sqrt((math.pi*a/Q))* \
                                  (M1*g0*g0*f0*fw)) # g0*g0 on purpose!
            KI.append((St+Sb)*math.sqrt((math.pi*a/Q))* \
                                  (M1*g90*g90*f90*fw)) # g90*g90 on purpose!
            
        elif a > c:
            Q=1.0+1.464*(c/a)**1.65
            M1=(math.sqrt(c/a))*(1.08-0.03*c/a)
            g0=1.08 # assuming t=large! 
            g90=1.0 # assuming t=large!
            f0=1.0 # g(0)=g(180)
            f90=math.sqrt(c/a)
            fw=1.0 # t=large, b=large
            KI.append((St+Sb)*math.sqrt((math.pi*a/Q))* \
                                  (M1*g0*f0))
            KI.append((St+Sb)*math.sqrt((math.pi*a/Q))* \
                                  (M1*g90*f90))

        return JVT.ScalarMult(KI,scale)

######## Qellipse

    def GeometryCheck(self,CMesh):
        '''
        no further geometry to sink to...
        ''' 
        return 0,0

######## Qellipse

    def GiveCurrent(self):
        '''
        Returns the most current a, b, y0, z0.
        '''
        a = self.a[0][0][-1]
        b = self.a[1][0][-1]

        return a,b

######## Qellipse

    def IntegrateLife(self,ai,af,N_max,Life,tol):
        '''
        Using simpson's rule, integrate for life between ai and af.  For now,
        returns shortest life from the four possible positions.

        ai - initial crack length
        af - final crack length (usually float, but can be str()='end')
        r - increment factor (ai+1 = r * ai)
        '''

        # check the length of the SIF history... 
        if len(self.a[1]) == 0: return 0.0

        N=N_max
        Xlist,Klist=self.Lists(ai,af,self.a[0],self.a[1],tol)
        NN=self.Simpson(Xlist,Klist,Life,self.dadN[0],tol)
        if NN < N and NN >= 0:
            N = NN
        elif NN == -1:
            N = -1

        Xlist,Klist=self.Lists(ai,af,self.b[0],self.b[1],tol)
        NN=self.Simpson(Xlist,Klist,Life,self.dadN[1],tol)
        if NN < N and NN >= 0:
            N = NN
        elif NN == -1:
            N = -1

        return N

######## Qellipse

    def ToFile(self,path_name,did): # change for Gerd
        '''
        function to print information about Fellipse to a file.
        '''

        if len(self.DamHistory[0])>0:
            a_his = open(path_name+'''.a''', 'a')
            b_his = open(path_name+'''.b''', 'a')
            dN = open(path_name+'''.dN''', 'a')
            for i in xrange(len(self.DamHistory[0])):
                a_his.write(str(did)+' '+str(self.DamHistory['ab'][i][0])+  \
                                     ' '+str(self.DamHistory[0][i])+"\n")
                b_his.write(str(did)+' '+str(self.DamHistory['ab'][i][1])+  \
                                     ' '+str(self.DamHistory[90][i])+"\n")
                try: dN.write(str(did)+' '+str(self.DamHistory['dN'][i])+"\n")
                except IndexError: pass

        else:
            return

######## Qellipse

    def UpdateState(self,anew):

        self.a[0][0]+=[anew[0]]
        self.a[1][0]+=[anew[1]]

        # note: self.Trans can't move now... 

######## Qellipse

    def __ArcLength(self,phi_1,phi_2):
        '''
        A nine point Gauss quadrature that calculates the arc length on the
        ellipse between phi_a and phi_b.
        '''

        # get current major and minor axis dimensions
        a,b=self.GiveCurrent()

        # define the arclength lambda
        arc=lambda x: math.sqrt(a*a*((math.sin((x+1)*math.pi)/4)**2) + \
                                b*b*((math.cos((x+1)*math.pi)/4)**2))

        # define gauss weights and gauss points
        gp9 = [
        -0.968160239,
        -0.836031107,
        -0.613371432,
        -0.324253423,
         0.0,
         0.324253423,
         0.613371432,
         0.836031107,
         0.968160239]

        gw9 = [
        0.081274388,
        0.180648160,
        0.260610696,
        0.312347077,
        0.330239355,
        0.312347077,
        0.260610696,
        0.180648160,
        0.081274388]

        # loop over gp9 and sum
        nsum=0

        for i in xrange(len(gp9)):
            ri = phi_1+(phi_2-phi_1)*(gp9[i]+1.0)/2.0
            wi = (phi_2-phi_1)*(gw9[i])/2
            y = arc(ri) # arc is the lambda defined above
            nsum+=wi*y

        return nsum*math.pi     

######## Qellipse

    def __BisectPhi(self,phi_L,phi_H,rho,tol):
        '''
        Finds the intersection of the ellipse with the surface of the FEM
        mesh.
        '''

        a,b=self.GiveCurrent()

        # phi_L
        y,z=self.__yzfunc(a,b,rho,phi_L)
        u=Numeric.dot(self.TRotation,[0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        (flag,dist)=self.PointOutside(qpnt,self.model)
        if flag == 1:
            flag_L = 0
        else:
            flag_L = 1

        # phi_H
        y,z=self.__yzfunc(a,b,rho,phi_H)
        u=Numeric.dot(self.TRotation,[0,y,z])
        qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
        (flag,dist)=self.PointOutside(qpnt,self.model)
        if flag == 1:
            flag_H = 0
        else:
            flag_H = 1

        while abs(phi_H-phi_L) > tol:
            phi_M=0.5*(phi_H+phi_L)
            y,z=self.__yzfunc(a,b,rho,phi_M)
            u=Numeric.dot(self.TRotation,[0,y,z])
            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
            (flag,dist)=self.PointOutside(qpnt,self.model)
            if flag == 1:
                flag_M = 0
            else:
                flag_M = 1

            if flag_M == flag_H:
                phi_H = phi_M
            else:
                phi_L = phi_M

        # flag_H indicates whether avg(phi_L,phi_H) is a transition from or
        # into the volume.  i.e., if flag_H = 0 it is a transition OUT of the
        # volume.  Conversly, if flag_H = 1 it is a transition INTO the volume.
        return phi_L, phi_H, flag_H

######## Qellipse

    def __FitPolynomial(self,aa,bb,center): 
        ypnts=[0.0,aa/2.0,aa]
        zpnts=[0.0,bb/2.0,bb]
        f=[]
        XY=[]
        A = self.Rotation[0][0]
        B = self.Rotation[0][1]
        C = self.Rotation[0][2]
        # get stress at each guass point
        for y in ypnts:
            for z in zpnts:
                #u=Numeric.dot(self.TRotation,[0,y,z])
                u0=self.Rotation[1][0]*y+self.Rotation[2][0]*z
                u1=self.Rotation[1][1]*y+self.Rotation[2][1]*z
                u2=self.Rotation[1][2]*y+self.Rotation[2][2]*z
                #qpnt=center+Vec3D.Vec3D(u[0],u[1],u[2])
                qpnt=center+Vec3D.Vec3D(u0,u1,u2)
                sig=self.model.GetPtStress(qpnt)
                # GetPtStress returns a ColTensor, arrange to use Numeric the 
                # global stress tensor is...
                sigG=[[sig.xx(), sig.xy(), sig.zx()],
                      [sig.xy(), sig.yy(), sig.yz()],
                      [sig.zx(), sig.yz(), sig.zz()]]
                # sigL = R*sigG*R'
                #R_sigG=Numeric.dot(self.Rotation,sigG)
                #sigL=Numeric.dot(R_sigG,self.TRotation)

                # fill matrix and right-hand-side vector
                #f+=[sigL[0][0]]
                XY+=[[z, 1]]
                f+=[A*A*sig.xx() + B*B*sig.yy() + C*C*sig.zz() \
                   +2.0*A*B*sig.xy()+2.0*A*C*sig.zx()+2.0*B*C*sig.yz()]

        # add the center of ellipse...
        sig=self.model.GetPtStress(center)
        sigG=[[sig.xx(), sig.xy(), sig.zx()],
              [sig.xy(), sig.yy(), sig.yz()],
              [sig.zx(), sig.yz(), sig.zz()]]
        #R_sigG=Numeric.dot(self.Rotation,sigG)
        #sigL=Numeric.dot(R_sigG,self.TRotation)
        #f+=[sigL[0][0]]
        f+=[A*A*sig.xx() + B*B*sig.yy() + C*C*sig.zz() \
            +2.0*A*B*sig.xy()+2.0*A*C*sig.zx()+2.0*B*C*sig.yz()]
        XY+=[[0.0, 1]]

        # solve the linear least squares problem (method of normal equations)
        XYTXY=0.; XYTf=0.
        XYTXY=Numeric.dot(Numeric.transpose(XY),XY)
        XYTf=Numeric.dot(Numeric.transpose(XY),f)
        coef=LinearAlgebra.solve_linear_equations(XYTXY,XYTf)

        return coef

######## Qellipse

    def __FixRotation(self):
        '''
        to make sure b is in "thickness" direction of the plate and a is along
        "width" of the plate, per Raju & Neman's assumptions.
        '''

        # get structure's dimension in the e2 direction (what i call 'a'
        # dimension of ellipse).  
        vec=(Vec3D.Vec3D(0,1,0.1)).Normalize()
        La=self.SeekLongDimension(vec,self.Rotation,self.model,self.Trans)
        # get structure's dimension in the e3 direction
        vec=(Vec3D.Vec3D(0,0.1,1)).Normalize()
        Lb=self.SeekLongDimension(vec,self.Rotation,self.model,self.Trans)

        if Lb > La: # switch e2 & e3 and recalc Rotation
            x2=Vec3D.Vec3D(self.Rotation[2][0],self.Rotation[2][1],\
                           self.Rotation[2][2])
            x3=Vec3D.Vec3D(self.Rotation[1][0],self.Rotation[1][1],\
                           self.Rotation[1][2])
            x1=Vec3D.CrossProd(x2,x3).Normalize()
            self.Rotation=Numeric.array([[x1.x(), x1.y(), x1.z()],
                                         [x2.x(), x2.y(), x2.z()],
                                         [x3.x(), x3.y(), x3.z()]])

######## Qellipse

    def __InOrOut(self,phi_in,inc):
        '''
        Start at phi_in (parametric coord for ellipse) and increase phi
        until the status of the ellipse (according to FemModel.IsPointIn())
        changes.  This effectively tells you that you have moved into or out of
        the FEM mesh.
        '''

        # get current ellipse info
        a,b=self.GiveCurrent()

        # check the phi_in first
        y,z=self.__yzfunc(a,b,1.0,phi_in)
        u=Numeric.dot(self.TRotation,[0.0,y,z])
        qpnt=self.DamOro+Vec3D.Vec3D(u[0],u[1],u[2])

        # FemModel.IsPointIn() returns 1 if point is inside the mesh and 0 if
        # the point is NOT in the mesh.
        flag,eid,eclass,coords,dist=self.model.IsPointIn(qpnt)

        switch=flag

        # the while loop actually only forms a bounding interval for phi in
        # which the point of interests lies.  __BisectPhi finds the point to 
        # the gnats ass.
        while switch==flag and phi_in < 1.001:
            flag=switch
            phi_in+=inc # incrementally increase phi_in
            y,z=self.__yzfunc(a,b,1.0,phi_in)
            u=Numeric.dot(self.TRotation,[0,y,z])
            qpnt=self.DamOro+Vec3D.Vec3D(u[0],u[1],u[2])
            switch,eid,eclass,coords,dist=self.model.IsPointIn(qpnt)

        if phi_in >= 1.0:
            phi_L=None
            phi_H=None
            flag=None
        else:
            tol=0.001 # arbitrarily chosen for now 1/13/04
            phi_L,phi_H,flag=self.__BisectPhi(phi_in-inc,phi_in,1.0,tol)

        return phi_L,phi_H,flag

######## Qellipse

    def WillGrow(self,N,R,Ki):
        '''
        Use this if Qellipse is appropriate, it calc's Ki and compares to
        specific material properties.

        Use Kic for failure criteria since the NASGRO equation for Kc includes
        thickness, t, that will be hard to estimate here.  
        '''

        RANGE=range(2)

        if Ki == 4: return 4

        # store for use elsewhere...
        for i in RANGE: self.a[i][1]+=[Ki[i]]

        #  commented out to speedup... logically shouldn't need to check this
        # here anyway?! jme 7/15/06
##        # check 5 points around quarter-ellipse
##        a,b=self.GiveCurrent()
##        flist=[]
##        for i in range(5):
##            phi=float(i)/2.0 - 1.0
##            y,z=self.__yzfunc(a,b,1.0,phi)
##            u=Numeric.dot(self.TRotation,[0,y,z])
##            qpnt=self.Trans+Vec3D.Vec3D(u[0],u[1],u[2])
##            (flag,dist)=self.PointOutside(qpnt,self.model)
##            if flag == 1:
##                flist+=[0]
##            else:
##                flist+=[1]
##
##        # list contains the indicator of whether the 5 points are in model or
##        # not.  If they are all zero, assume no remaining ligament.  Set
##        # will_grow = 2.
##        if flist[0]==0 and flist[1]==0 and flist[2]==0 and \
##           flist[3]==0 and flist[4]==0:
##            return 4

        Kic=[self.dadN[i].Kic for i in RANGE]
        Reff=[self.dadN[i].CalcReff(int(N),R) for i in RANGE]
        Kth=[self.dadN[i].Calc_DKth(self.a[i][0][-1],Reff[i]) \
             for i in RANGE]

        Kmax=[Ki[i]/(1.0-R) for i in RANGE]

        # max K is negative = compression, will not grow.
        if max(Ki) <= 0: 
            return -1 

        # will not grow, even under cyclic loading if each Ki is less than Kth
        # at each point on ellipse (AND)...
        elif Ki[0] <= Kth[0] and Ki[1] <= Kth[1]: 
            return 0

        # as in unstable growth if...
        elif Kmax[0] >= Kic[0] or Kmax[1] >= Kic[1]: 
            return 2

        else: return 1 

######## Qellipse

    def __init__(self,doid,center,model,a,b,norm,material,Rotation=None, \
                 verbose=None,verification=False):
        self.a=[DamHistory.DamHistory(a),DamHistory.DamHistory(b)]
        self.CrackFrontPoints=2
        self.doid = doid
        self.names=['Qellipse','a','b']
        self.verification = verification
        self.model = model
        self.verbose=verbose
        self.Trans = center # location of cntr, update as necessary. 
                          # original center is always retrievable by
                          # xyz,delxyz,sigxyz=model.GetNodeInfo(self.doid) 
        self.dN=[] # list of time steps Numeric.sum(self.dN) = Life
        self.nextdN = 0.0

        # initiate dadN model ** don't bother storing material as attribute! 
        self.dadN=[dadN.dadN(material,1),dadN.dadN(material,1)]
        # initialize the plastic zone size
        for i in range(self.CrackFrontPoints): self.dadN[i].SetKol()

        # x1 is in direction of first principal stress
        x1=(Vec3D.Vec3D(Rotation[0][0], Rotation[0][1], \
                        Rotation[0][2])).Normalize()
        # x2 is normal to x1 and norm
        x2=(Vec3D.CrossProd(x1,norm)).Normalize()
        # x3 by cross x1 X x2...
        x3=(Vec3D.CrossProd(x1,x2)).Normalize()
        # rotation matix then is...
        self.Rotation=[[x1.x(), x1.y(), x1.z()],
                       [x2.x(), x2.y(), x2.z()],
                       [x3.x(), x3.y(), x3.z()]]

        # The corner is assumed to by 90 degrees so the rotation matrix based
        # on the eigen values and norm must be rotated 45 degrees...
        c=math.cos(math.pi/4.0)
        r=[[1.0, 0., 0.], [0., c, c], [0., -c, c]]
        # the matrix product of the two gives the rotation matrix from global
        # coords to the coordinate system of the quarter ellipse.
        # A'= r R A Rt rt (Rt is R transpose)
        # so self.Rotation = rR ...
        self.Rotation=Numeric.dot(r,self.Rotation)

        # Raju & Newman's solution for a quarter ellipse in a plate with 
        # dimensions w * h * t, expects 'a' to be on w face and 'b' to be in
        # the t direction.  (note R&N's solutios calls uses c for what i call
        # 'a' and uses a for what i call 'b').  The following method checks
        # changes self.Rotation so that 'b' is assigned the shortest direction
        self.__FixRotation()

        # store transpose(self.Rotation) to reduce calls to Numeric.Transpose
        self.TRotation = Numeric.transpose(self.Rotation)

        if verification:
            self.verification=True
            print self
            print ' ',self.names[0], ' Rotation matrix is:'
            print repr(Numeric.array(self.Rotation))+"\n"

######## Qellipse

    def __yzfunc(self,a,b,rho,phi):
        '''
        define mapping from rectanglular coords to the ellipse
        '''
        y=a*((rho+1)/2)*math.cos(math.pi*(phi+1)/4)
        z=b*((rho+1)/2)*math.sin(math.pi*(phi+1)/4)
        return y,z

######## Qellipse

    def dAdN(self,N,size,args):
        '''
        Used by GrowDam in DamMo.py
        N - current number of cycles of loading
        size - a list of ellipses dimensions: [a, b, -a, -b]
        args - tuple of supporting arguments, in this case just scale! 
        '''
        RANGE=range(2)
        # unpack args
        scale=args[0]
        growth = []
        abab= [max(0,size[i]) for i in RANGE]
        KIs = self.__CalculateKi(abab[0],abab[1],self.Trans,scale)
        for i in RANGE:
            growth.append(self.dadN[i].Calc_dadN(KIs[i],size[i],int(N)))

        if max(growth) <= 0.:
            raise DamErrors.dAdNError, ('Qellipse',\
                                  'All growth rates <= 0.0',growth)
        return growth

########################################################################
########################################################################

def TwoBricks():
    blah=open("check.smp","w")
    blah.write("1 Brick1IShapeFunc")
    blah.close()

    blah0=open("check.nod","w")
    blah1=open("check.con","w")
    blah2=open("check.sig","w")
    nodes={0:((0,0,0),(10,0,0,0,0,0)),1:((0,0,1),(10,0,0,0,0,0)), \
           2:((0,1,0),(10,0,0,0,0,0)),3:((0,1,1),(10,0,0,0,0,0)), \
           4:((1,0,0),(10,0,0,0,0,0)),5:((1,0,1),(10,0,0,0,0,0)), \
           6:((1,1,0),(10,0,0,0,0,0)),7:((1,1,1),(10,0,0,0,0,0)), \
           8:((0,2,0),(10,0,0,0,0,0)),9:((0,2,1),(10,0,0,0,0,0)), \
           10:((1,2,0),(10,0,0,0,0,0)),11:((1,2,1),(10,0,0,0,0,0))}
    
    for i in range(len(nodes)):
        blah0.write(str(i)+' '+str(nodes[i][0][0])+' '+str(nodes[i][0][1])+ \
                    ' '+str(nodes[i][0][2]))
        blah2.write(str(i)+' '+str(nodes[i][1][0])+' '+str(nodes[i][1][1])+ \
                    ' '+str(nodes[i][1][2])+' '+str(nodes[i][1][3])+' '+ \
                    str(nodes[i][1][4])+' '+str(nodes[i][1][5]))
        if i != len(nodes)-1:
            blah0.write("\n")
            blah2.write("\n")

    blah1.write("0 1 0 8 0 1 2 3 4 5 6 7\n")
    blah1.write("1 1 0 8 2 3 8 9 6 7 10 11")

    blah0.close()
    blah1.close()
    blah2.close()

    cwd=os.getcwd()

    print 'calling meshtools'
    model=MeshTools.MeshTools(cwd+"\\check",'RDB')
    print 'meshtools was successfull'
    os.system("del check.*")

    return model

#########################

def OneTet2I():
    blah=open("check.smp","w")
    blah.write("1 Tet2IShapeFunc \n")
    blah.close()

    blah0=open("check.nod","w")
    blah1=open("check.con","w")
    blah2=open("check.sig","w")
    blah3=open("check.edg","w")
    nodes={0:((0,0,0),(10,0,0,0,0,0)),1:((0,0,1),(10,0,0,0,0,0)), \
           2:((1,0,0),(10,0,0,0,0,0)),3:((0,1,0),(10,0,0,0,0,0)), \
           4:((0,0,0.5),(10,0,0,0,0,0)),5:((0.5,0,0.5),(10,0,0,0,0,0)), \
           6:((0.5,0,0),(10,0,0,0,0,0)),7:((0,0.5,0.5),(10,0,0,0,0,0)), \
           8:((0.5,0.5,0),(10,0,0,0,0,0)),9:((0,0.5,0),(10,0,0,0,0,0))}
    
    for i in range(len(nodes)):
        blah0.write(str(i)+' '+str(nodes[i][0][0])+' '+str(nodes[i][0][1])+ \
                    ' '+str(nodes[i][0][2]))
        blah2.write(str(i)+' '+str(nodes[i][1][0])+' '+str(nodes[i][1][1])+ \
                    ' '+str(nodes[i][1][2])+' '+str(nodes[i][1][3])+' '+ \
                    str(nodes[i][1][4])+' '+str(nodes[i][1][5]))
##        if i != len(nodes)-1:
##            blah0.write("\n")
##            blah2.write("\n")

    blah1.write("0 1 0 4 0 1 2 3\n")
    blah3.write("4 2 1 0 1\n")
    blah3.write("5 2 1 1 2\n")
    blah3.write("6 2 1 0 2\n")
    blah3.write("7 2 1 1 3\n")
    blah3.write("8 2 1 2 3\n")
    blah3.write("9 2 1 0 3\n")

    blah0.close()
    blah1.close()
    blah2.close()
    blah3.close()

    cwd=os.getcwd()

    print 'calling meshtools'
    model=MeshTools.MeshTools(cwd+"\\check",'RDB')
    print 'meshtools was successfull'
    os.system("del check.*")

    return model

if __name__=="__main__":

    
    material=[84.0, 75.0, 17.0, 12.0, 1.0, 1.0, 2.09e-005, 1.8999999999999999,\
              0.5, 1.0, 2.0, 2.0, 0.10000000000000001, 0.69999999999999996,  \
              1.8999999999999999, 0.29999999999999999, 0.0, 0.0015, \
              2.2999999999999998, 100.0]

    cwd=os.getcwd()
    model=MeshTools.MeshTools(cwd+"\\examples\\example1",'RDB')
##    model=OneTet2I()
##    print model.GetPtStress(Vec3D.Vec3D(0.25,0.25,0.25))

    blah=Fellipse(8,Vec3D.Vec3D(0,0,0),model,0.1,0.1,material)
    print blah







        
