# John E.
# 4/12/07
# script to read stiffness matrix & rhs and solution from FEAWD,
# compute the solution myself and compare the two.  

import sys
import math
import string
import Numeric
import LinearAlgebra

######################

def ReadVecfile(Vecfile):

    buff=Vecfile.readline() # throw-away first line
    buff=Vecfile.readline()
    vec=[]
    while buff:
        vec+=[float(buff)]
        buff=Vecfile.readline()

    return Numeric.array(vec)

######################

def ReadKfile(Kfile,DOF):

    K=Numeric.zeros((DOF,DOF))
    buff=Kfile.readline()
    while buff:
        # strip out :
        a=string.split(buff,':')

        # assign row
        row=string.split(a[0])
        row=int(row[1])

        # manipulate the rest of the line so have list of pairs 'col, float'
        b=string.strip(a[1],' (')
        b=b[:-3]
        b=string.split(b,')  (')
        for L in b:
            v=string.split(L,',')
            col=int(v[0])
            K[row][col]=float(v[1])

        buff=Kfile.readline()

    return K

######################

if __name__=="__main__":

    try:
        print "basing number of equations on length of rhs.dump file!" 
        forcefile = open('rhs.dump','r')
        force = ReadVecfile(forcefile)
        forcefile.close()
    except IOError:
        print " no mat.dump file in this directory.  Exiting"
        sys.exit()

    try:
        Kfile = open('mat.dump','r')
        K = ReadKfile(Kfile,len(force))
        Kfile.close()
    except IOError:
        print " no mat.dump file in this directory.  Exiting"
        sys.exit()

    try:
        xfile = open('x.dump','r')
        HisX = ReadVecfile(xfile)
        xfile.close()
    except IOError:
        print " no mat.dump file in this directory.  Exiting"
        sys.exit()

    print 'row 12'
    print K[12][:]
    print 'row 56'
    print K[56][:]
    print 'row 123'
    print K[123][:]

    print 'row 12'
    print force[12]
    print 'row 56'
    print force[56]
    print 'row 123'
    print force[123]

    MyX = LinearAlgebra.solve_linear_equations(K,force)

    print MyX

    print HisX-MyX

    out = open('MyX.dump','w')
    for row in MyX: out.write(str(row)+"\n")
    out.close()