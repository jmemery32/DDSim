import Numeric, sys, string, LinearAlgebra

def PopulateK(K,Kname):
    '''
    read Kfilename.txt and store elements of K appropriately.
    '''

    file=open(Kname,'r')
    buff=file.readline()
    while buff:
        line=string.splitfields(buff)
        i=int(line[0])
        j=int(line[1])
        K[i][j]=float(line[2])
        buff=file.readline()

########################################################

def Permute(K,d,p,name):
    '''
    read ijkgfilename.txt, build permutation matrix and permute K.
    '''

    file=open(name,'r')
    buff=file.readline()
    dof=len(K)
    gamma=Numeric.zeros([dof,dof])
    counter=0
    dofs={}
    while buff:
        line=string.splitfields(buff)
        buff=file.readline()
        dofid=[]
        for blah in line[2:]:
            i=int(blah)
            dofid+=[i]
            entry=Numeric.zeros(dof)
            entry[i]=1.0
            gamma[counter]=entry
            counter+=1

        if line:
            if line[0]=='i': dofs['i']=dofid
            elif line[0]=='j': dofs['j']=dofid
            elif line[0]=='k': dofs['k']=dofid
            elif line[0]=='gamma': dofs['g']=dofid
            else: pass

    K=Numeric.dot(gamma,K)
    K=Numeric.dot(K,Numeric.transpose(gamma))

    p=Numeric.dot(gamma,p)
    d=Numeric.dot(gamma,d)

    return K,d,p,dofs,gamma

########################################################

def G5(K,p,A,dofs,KK,pp):
    '''
    manipulate equations into G5.
    '''

    I=range(0,len(dofs['i']))
    J=range(I[-1]+1,len(dofs['i'])+len(dofs['j']))
    KS=range(J[-1]+1,len(dofs['i'])+len(dofs['j'])+len(dofs['k']))
    G=range(KS[-1]+1,len(dofs['i'])+len(dofs['j'])+len(dofs['k'])+len(dofs['g']))

    # zeros and I
    for row in KS:
        for col in range(len(K)):
            K[row][col]=0.0
    for col in KS:
        for row in range(len(K)):
            K[row][col]=0.0
    for row in KS:
        K[row][row]=1.0

    # manipulate I,J
    for i in range(len(I)):
        row=I[i]
        for j in range(len(J)):
            col=J[j]
            for k in range(len(KS)):
                colcol=KS[k]
                K[row][col]+=KK[row][colcol]*A[k][j]

    # J,J
    for i in range(len(J)):
        row=J[i]
        for j in range(len(J)):
            col=J[j]
            for k in range(len(KS)):
                colcol=KS[k]
                K[row][col]+=KK[row][colcol]*A[k][j]

    # K,J
    for i in range(len(KS)):
        row=KS[i]
        for j in range(len(J)):
            col=J[j]
            K[row][col]=-A[i][j]

    # G,J
    for i in range(len(G)):
        row=G[i]
        for j in range(len(J)):
            col=J[j]
            for k in range(len(KS)):
                colcol=KS[k]
                K[row][col]+=KK[row][colcol]*A[k][j]

    # now take care of p.  B=[0] so, only have to zero k's of p
    for k in KS:
        p[k]=0.0

    return K,p

########################################################

def G6(K,p,A,dofs,KK,pp):
    '''
    manipulate equations into G6.
    '''

    I=range(0,len(dofs['i']))
    J=range(I[-1]+1,len(dofs['i'])+len(dofs['j']))
    KS=range(J[-1]+1,len(dofs['i'])+len(dofs['j'])+len(dofs['k']))
    G=range(KS[-1]+1,len(dofs['i'])+len(dofs['j'])+len(dofs['k'])+len(dofs['g']))

    # J,I
    for i in range(len(J)):
        row=J[i]
        for j in range(len(I)):
            col=I[j]
            for k in range(len(KS)):
                rowrow=KS[k]
                K[row][col]+=A[k][i]*KK[rowrow][col]

    # J,J
    for i in range(len(J)):
        row=J[i]
        for j in range(len(J)):
            col=J[j]
            for k in range(len(KS)):
                rowrow=KS[k]
                K[row][col]+=A[k][i]*KK[rowrow][col]

    # J,K
    for i in range(len(J)):
        row=J[i]
        for j in range(len(KS)):
            col=KS[j]
            for k in range(len(KS)):
                rowrow=KS[k]
                K[row][col]+=A[k][i]*KK[rowrow][col]

    # J,G
    for i in range(len(J)):
        row=J[i]
        for j in range(len(G)):
            col=G[j]
            for k in range(len(KS)):
                rowrow=KS[k]
                K[row][col]+=A[k][i]*KK[rowrow][col]

    # now take care of p.  B=[0] so, only have to add trans(A)*Fk to j's
    for j in range(len(J)):
        prow=J[j]
        for k in range(len(KS)):
            pprow=KS[k]
            p[prow]+=A[k][j]*pp[pprow]

    return K,p

########################################################

def G7(K,p,A,dofs,KK,pp):
    '''
    manipulate equations into G7.
    '''

    I=range(0,len(dofs['i']))
    J=range(I[-1]+1,len(dofs['i'])+len(dofs['j']))
    KS=range(J[-1]+1,len(dofs['i'])+len(dofs['j'])+len(dofs['k']))
    G=range(KS[-1]+1,len(dofs['i'])+len(dofs['j'])+len(dofs['k'])+len(dofs['g']))

    # K,J
    for i in range(len(KS)):
        row=KS[i]
        for j in range(len(J)):
            col=J[j]
            sum=0.0
            for k in range(len(KS)):
                colcol=KS[k]
                sum+=KK[row][colcol]*A[k][j]
            K[row][col]=sum

    # K,K
    for i in range(len(KS)):
        row=KS[i]
        for j in range(len(KS)):
            col=KS[j]
            K[row][col] = -1.0*KK[row][col]

    # now take care of p.  B=[0] so, nothing to do here! 

    return K,p

########################################################

def G8(K,p,A,dofs,KK,pp,dd):
    '''
    manipulate equations into G8.
    '''

    I=range(0,len(dofs['i']))
    J=range(I[-1]+1,len(dofs['i'])+len(dofs['j']))
    KS=range(J[-1]+1,len(dofs['i'])+len(dofs['j'])+len(dofs['k']))
    G=range(KS[-1]+1,len(dofs['i'])+len(dofs['j'])+len(dofs['k'])+len(dofs['g']))

    # do p first because we zero important stuff in K.  
    for i in range(len(I)):
        prow=I[i]
        for g in range(len(G)):
            pprow=G[g]
            p[prow]-=KK[prow][pprow]*dd[pprow]

    for j in range(len(J)):
        prow=J[j]
        for g in range(len(G)):
            pprow=G[g]
            p[prow]-=K[prow][pprow]*dd[pprow]

    for g in range(len(G)):
        prow=G[g]
        p[prow]=KK[prow][prow]*dd[prow]

    # zeros and Gamma 
    for row in G:
        for col in range(len(K)):
            K[row][col]=0.0
    for col in G:
        for row in range(len(K)):
            K[row][col]=0.0
    for row in G:
        K[row][row]=KK[row][row]

    return K,p

########################################################

def HelpPrints():
    print " " 
    print " throw-away code to test MPCs"
    print " ---- " 
    print " Arguments are:"
    print " ...>MPCTest.py [Kfilename.txt ijkgfilename.txt DOF] (optional flags)"
    print " "
    print " Kfilename.txt has format:"
    print "     i j K(i,j)"
    print " ijkgfilename.txt has format:"
    print "     i : 3 6 8... list of free DOF's "
    print "     j : 12 13... list of master DOF's"
    print "     k : 15 16... list of constrained DOF's"
    print "     gamma : 0 1 2 4 ... list of KBC DOF's" 
    print ' '
    print ' optional flags:'
    print ' -K to print stiffness matrix as read-in from data files'
    print ' -Perm to print stiffness after permutation into i,j,k,gamma'
    print ' -G5 to print the stiffness matrix at gerds equation 5'
    print ' -G6 to print the stiffness matrix at gerds equation 6'
    print ' -G7 to print the stiffness matrix at gerds equation 7'
    print ' -G8 to print the stiffness matrix at gerds equation 8'
    print ' -v for all' 
    print ' '
    sys.exit()

########################################################

if __name__=="__main__":
    '''
    throw-away code to test MPCs
    '''

    # a helper for the keys...
    if "-help" in sys.argv or "-h" in sys.argv:
        HelpPrints()

    try: 
        Kname=sys.argv[1]
        ijkgname=sys.argv[2]
        DOF=int(sys.argv[3])
    except: HelpPrints()

    # make K, d, and p s.t. [K]{d}={p}
    K=Numeric.zeros([DOF,DOF])
    p=Numeric.zeros(DOF)
    d=p.copy()

    # HARD CODING APPLIED DISPLACEMENTS
    d[37]=1.0
    d[40]=1.0
    d[43]=1.0
    d[46]=1.0

    PopulateK(K,Kname)
    if "-K" in sys.argv or "-v" in sys.argv:
        print "K"
        print K
        print "d"
        print d
        print "p"
        print p

    # Permute K, read bc's, 
    K,d,p,dofs,gamma = Permute(K,d,p,ijkgname)
    p_orig=p.copy()
    K_orig=K.copy()
    if "-Perm" in sys.argv or "-v" in sys.argv:
        print "permuted K"
        print K
        print "permuted d"
        print d
        print "permuted p"
        print p
        print "dofs"
        print dofs

    # constraint equation (assuming B=0)
    # make A --> v_k = A v_j... needs to be more sophisticated
    # for a real problem... 
    A=Numeric.identity(len(dofs['k']))
##    print "A"
##    print A

    # G5
    K,p = G5(K,p,A,dofs,K_orig,p_orig)
    if "-G5" in sys.argv or "-v" in sys.argv:
        print "G5, K, p"
        print Numeric.dot(Numeric.transpose(gamma),Numeric.dot(K,gamma))
        print Numeric.dot(Numeric.transpose(gamma),p)

    # G6
    K,p = G6(K,p,A,dofs,K_orig,p_orig)
    if "-G6" in sys.argv or "-v" in sys.argv:
        print "G6, K, p"
        print Numeric.dot(Numeric.transpose(gamma),Numeric.dot(K,gamma))
        print Numeric.dot(Numeric.transpose(gamma),p)

    # G7
    K,p = G7(K,p,A,dofs,K_orig,p_orig)
    if "-G7" in sys.argv or "-v" in sys.argv:
        print "G7, K, p"
        print Numeric.dot(Numeric.transpose(gamma),Numeric.dot(K,gamma))
        print Numeric.dot(Numeric.transpose(gamma),p)

    # G8
    K,p = G8(K,p,A,dofs,K_orig,p_orig,d)
    if "-G8" in sys.argv or "-v" in sys.argv:
        print "G8, K, p"
        print Numeric.dot(Numeric.transpose(gamma),Numeric.dot(K,gamma))
        print Numeric.dot(Numeric.transpose(gamma),p)

    # solve the linear system
    dd=LinearAlgebra.solve_linear_equations(K,p)
    dd=Numeric.dot(Numeric.transpose(gamma),Numeric.transpose(dd))
    print " displacements per dof are:"
    for i in range(len(dd)):
        print "dof #: %3i %1.5e" % (i, dd[i])

















