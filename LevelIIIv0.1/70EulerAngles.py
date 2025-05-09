import Vec3D
import math
import string

#  little script to make the Materials file for my thesis.
# (70 grains, 70 Euler angles, RPI plasticity).. 

def ComputeEulerAngles(rot):
    '''
    compute the three euler angles equivalent to the rotation matrix, rot.

    From C++ code borrowed from Wash (via JD)
    '''

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
                    return [t1[j],\
                            t2[k],\
                            t0[i]]

    return [0.0,0.0,0.0]

# ##########################################################################

def RodriguesToMtx(r):
    '''
    from wash to take Rodrigues angles and spit out rotation matrix.
    '''

    mag = r.Magnitude() 
    ang = 2.0 * math.atan(mag) 
    n = r.Normalize() 

    mtx=[[0.0,0.0,0.0],\
         [0.0,0.0,0.0],\
         [0.0,0.0,0.0]]

    if (ang == 0.0):
        for i in range(3):
            for j in range(3):
                mtx[i][j] = 0.0 
                mtx[i][i] = 1.0 

    else:
        mtx[0][0] = (1.0-n[0]*n[0])*math.cos(ang) + n[0]*n[0] 
        mtx[0][1] = n[0]*n[1]*(1.0-math.cos(ang)) + n[2]*math.sin(ang) 
        mtx[0][2] = n[0]*n[2]*(1.0-math.cos(ang)) - n[1]*math.sin(ang) 

        mtx[1][0] = n[1]*n[0]*(1.0-math.cos(ang)) - n[2]*math.sin(ang) 
        mtx[1][1] = (1.0-n[1]*n[1])*math.cos(ang) + n[1]*n[1] 
        mtx[1][2] = n[1]*n[2]*(1-math.cos(ang))   + n[0]*math.sin(ang) 

        mtx[2][0] = n[2]*n[0]*(1.0-math.cos(ang)) + n[1]*math.sin(ang) 
        mtx[2][1] = n[2]*n[1]*(1.0-math.cos(ang)) - n[0]*math.sin(ang) 
        mtx[2][2] = (1.0-n[2]*n[2])*math.cos(ang) + n[2]*n[2] 

    return mtx

# ##########################################################################

if __name__=="__main__":

    fd=open('Z:\\research\\Dissertation\\DDSimLI\\part10000.txt','r')
    fo=open('Materials','w')

    t1=' 0 16 hardening precipitation m 0.005 g_o 31.9 gamma_dot_o 1.0 G_o 17.4 g_s_o 36.3 gamma_dot_s 5.0e+10 omega 0.0 num_slip_sys 12 mu 4105 lambda 8833 eta 740 K 11077 phi1 '
##    t1=' 0 16 hardening precipitation m 0.005 g_o 220 gamma_dot_o 1.0 G_o 120 g_s_o 250 gamma_dot_s 5.0e+10 omega 0.0 num_slip_sys 12 mu 28300 lambda 60900 eta 5100 K 76376 phi1 '
    t3=' Phi '
    t5=' phi2 '
    

    buff=fd.readline()
    for i in range(70):
        line=string.split(buff)
        buff=fd.readline()

        rod=Vec3D.Vec3D(float(line[0]),float(line[1]),float(line[2]))
        eulers=ComputeEulerAngles(RodriguesToMtx(rod))

        t0=str(i)
        t2=str(eulers[0])
        t4=str(eulers[1])
        t6=str(eulers[2])+"\n"

        fo.write(t0+t1+t2+t3+t4+t5+t6)




















