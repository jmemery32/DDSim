
import Vec3D
import sys
import math


def GetFrontData(fd):
    npos = []
    Ks = []
    pts = []
    axes = []
    angles = []
    buff = fd.readline()
    buff = fd.readline()
    while len(buff) > 0:
        v = buff.split()
        npos.append(float(v[0]))
        Ks.append(Vec3D.Vec3D(float(v[1]),float(v[2]),float(v[3])))
        pts.append(Vec3D.Vec3D(float(v[4]),float(v[5]),float(v[6])))
        e1 = Vec3D.Vec3D(float(v[7]),float(v[8]),float(v[9]))
        e2 = Vec3D.Vec3D(float(v[10]),float(v[11]),float(v[12]))
        e3 = Vec3D.Vec3D(float(v[13]),float(v[14]),float(v[15]))
        axes.append((e1,e2,e3))
        angles.append(float(v[16]))
        buff = fd.readline()
    return npos,Ks,pts,axes,angles


def PrintFrontData(fd,npos,Ks,pts,axes,angles,desc,fname):

    print >> fd,"STEP_CRACK_DATA"
    print >> fd,"("
    print >> fd,"VERSION: 1"
    if desc != "": print >> fd,"DESCRIPTION:",desc
    if fname != "": print >> fd,"FILENAME:",fname
    print >> fd,"NUM_CRACK_FRONTS: 1"
    print >> fd,"  CRACK_FRONT: 0",len(npos)
    for i in xrange(len(npos)):
        print >> fd,"    ",npos[i],1,0,0
        print >> fd,"    ",Ks[i],0
        print >> fd,"    ",pts[i]
        print >> fd,"    ",axes[i][0]
        print >> fd,"    ",axes[i][1]
        print >> fd,"    ",axes[i][2]
        print >> fd,"    ",angles[i],0,0,0
        print >> fd,"    ",0,0,0
        print >> fd,"    ",0,0,0
        print >> fd,"    ",0,0,0
    print >> fd,")"


if __name__ == "__main__":

    npos,Ks,pts,axes,angles = GetFrontData(sys.stdin)

##    fd = open(fname,"w")
    PrintFrontData(sys.stdout,npos,Ks,pts,axes,angles,"","")

####    fin = sys.argv[1]
##    fd=sys.stdin
##    print fd
##    blah=fd.readline()
##    print blah
##    x
##    fd = open(fin,'r')
##    npos,Ks,pts,axes,angles = GetFrontData(fd)
##
##    fd = sys.stdout
##    PrintFrontData(fd,npos,Ks,pts,axes,angles,fname + " DATA",fname)

