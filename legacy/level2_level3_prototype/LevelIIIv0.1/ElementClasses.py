
import sys
import Vec3D
import ColTensor


class Node:

    def __init__(self,coords,disp,stress):
        self.coords = Vec3D.Vec3D(coords[0],coords[1],coords[2])

        if disp:
            self.disp = Vec3D.Vec3D(disp[0],disp[1],disp[2])
        else:
            self.disp = None

        if stress:
            self.stress = ColTensor.ColTensor(stress[0],stress[1],
                                              stress[2],stress[3],
                                              stress[4],stress[5])
        else:
            self.stress = None

def FindNaturalCoords(self,eid,eclass,node_coords,query_pt):

    MAX_ITS = 10
    tol = 0.00001

    # set the natural coordinates to a known bad value

##        print (node_coords)

    natural = Vec3D.Vec3D(1e26,1e26,1e26)
    ncoord = eclass.CenterCoordinates()

    # newton iterations

    for num_its in xrange(MAX_ITS):

        shape = eclass.ShapeFunction(ncoord)
        deriv = eclass.ShapeDerivatives(ncoord)

        # interpolate to find the Cartesion coordinates
        # at the current trial point

        guess = Vec3D.Vec3D(0.0,0.0,0.0)
        for i in xrange(len(shape)):
            guess += shape[i] * node_coords[i]

        # compute the jacobian at this point

        jac = [[0,0,0],[0,0,0],[0,0,0]]
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(len(shape)):
                    jac[i][j] += deriv[i][k] * node_coords[k][j]

##            print (eid)
##            print (num_its)
##            print (guess)
##            print (shape)
##            print (deriv)

        # and the determinate
                    
        det = jac[0][0] * jac[1][1] * jac[2][2] + \
              jac[0][1] * jac[1][2] * jac[2][0] + \
              jac[0][2] * jac[1][0] * jac[2][1] - \
              jac[2][0] * jac[1][1] * jac[0][2] - \
              jac[2][1] * jac[1][2] * jac[0][0] - \
              jac[2][2] * jac[1][0] * jac[0][1]

##            print (det)

        # and the inverse

        inv = [[0,0,0],[0,0,0],[0,0,0]]
        inv[0][0] =  (jac[1][1]*jac[2][2] - jac[1][2]*jac[2][1]) / det
        inv[0][1] = -(jac[0][1]*jac[2][2] - jac[0][2]*jac[2][1]) / det
        inv[0][2] =  (jac[0][1]*jac[1][2] - jac[0][2]*jac[1][1]) / det
        inv[1][0] = -(jac[1][0]*jac[2][2] - jac[1][2]*jac[2][0]) / det
        inv[1][1] =  (jac[0][0]*jac[2][2] - jac[0][2]*jac[2][0]) / det
        inv[1][2] = -(jac[0][0]*jac[1][2] - jac[0][2]*jac[1][0]) / det
        inv[2][0] =  (jac[1][0]*jac[2][1] - jac[1][1]*jac[2][0]) / det
        inv[2][1] = -(jac[0][0]*jac[2][1] - jac[0][1]*jac[2][0]) / det
        inv[2][2] =  (jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0]) / det

        det = abs(det)

        # update our guessed natural coordinates based on
        # a Newton step

        delta = guess - query_pt

        update = Vec3D.Vec3D(0.0,0.0,0.0)
        for i in xrange(3):
            update[i] = inv[0][i]*delta[0] + \
                        inv[1][i]*delta[1] + \
                        inv[2][i]*delta[2]
        new_nc = ncoord - update

        # check for convergance
##            print (delta)
##            print (jac)
##            print (inv)
##            print (ncoord)
##            print (new_nc)
##            print (update)
##            print ('---')

        if abs(new_nc[0]-ncoord[0]) > tol or \
           abs(new_nc[1]-ncoord[1]) > tol or \
           abs(new_nc[2]-ncoord[2]) > tol:
            ncoord = new_nc
        else:
            natural = new_nc
            break

    return ncoord

#------------------------------------------------------
# Element Classes
#------------------------------------------------------
#
#  The Element classes support the various element types.
#  They all have the similar methods and variables:
#
#  Public Instance Variables:
#
#    order     - polynomial interpolation order, current
#                implementation supports 1 (linear) or 2 (quadratic)
#
#    mat       - element's material id
#
#    num_nodes - number of element nodes.
#
#    nodes     - list of the element's node id's
#
#  Public Methods:
#
#    "constructor" (order,mat,nodes)
#
#        Element class object constructor
#
#            order - polynomial interpolation order
#            mat   - material id
#            nodes - list of node id's
#
#
#    ShapeFunction(natural_coords)
#
#        natural_coords - Vec3D of the evaluation location
#
#        returns a list of the shape functions evaluated at the
#        input point.
#
#
#    ShapeDerivatives(natural_coords)
#
#        natural_coords - Vec3D of the evaluation location
#
#        returns a list of lists of the shape function derivatives
#        evaluated at the input point ([[dN/dr],[dN/ds],[dN/dt]]).
#
#
#    PointInside(natural_coords)
#
#        natural_coords - Vec3D of the evaluation location
#
#        returns 1 if the given natural coordinate location falls
#        within the element, 0 otherwise.
#
#
#    NearestPoint(natural_coords)
#
#        natural_coords - Vec3D of the evaluation location
#
#        returns a Vec3D of the natural coordinates for the
#        point nearest to the input that falls within or on
#        the surface of the element.
#
#
#    CenterCoordinates()
#
#        returns the natural coordinates for the point at the
#        center of the element.
#
#    GiveSurf()
#
#        returns the surfaces of an element oriented s.t. their surface
#        normal points outward.


class Tetrahedral:

#   4 & 10-noded numbering
#
#                            
#                         2 o +++
#                          / \   +++   7     
#                         /   \     ++o    
#                        /     \       +++ o 3 
#                       /       \5     ++++
#                    6 o         o ++++  +
#                     /     ++o++ \     o 8
#                    / +++++  9    \   +
#                   /++             \ +
#                  o--------o--------o  - - - - - - - - t = -1
#                   0       4        1


    def __init__(self,order,mat,nodes):
        self.order = order
        self.mat = mat
        if order == 1:
            self.num_node = 4
        else:
            self.num_node = 10
        assert len(nodes) >= self.num_node
        self.nodes = [nodes[i] for i in xrange(self.num_node)]
##        print (self.nodes)

    def ShapeFunction(self,nc):
        u = 1.0 - nc[0] - nc[1] - nc[2]
        if self.order == 1:
            shape = [nc[0],nc[1],nc[2],u]
        else:
            shape = [(2 * nc[0] - 1) * nc[0],
                     (2 * nc[1] - 1) * nc[1],
                     (2 * nc[2] - 1) * nc[2],
                     (2 * u - 1) * u,
                      4 * nc[0] * nc[1],
                      4 * nc[1] * nc[2],
                      4 * nc[0] * nc[2],
                      4 * nc[2] * u,
                      4 * nc[1] * u,
                      4 * nc[0] * u]
        return shape

    def ShapeDerivatives(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        u = 1.0 - r - s - t
        if self.order == 1:
            deriv = [
                [ 1.0, 0.0, 0.0,-1.0],
                [ 0.0, 1.0, 0.0,-1.0],
                [ 0.0, 0.0, 1.0,-1.0]]
        else:
            deriv = [
                [4.0*r - 1.0, 0.0, 0.0, 1.0 - 4.0*u, 4.0*s,
                 0.0, 4.0*t, -4.0*t, -4.0*s, 4.0 * (u - r)],
                [0.0, 4.0*s - 1.0, 0.0, 1.0 - 4.0*u, 4.0*r,
                 4.0*t, 0.0, -4.0*t, 4.0 * (u - s), -4.0*r],
                [0.0, 0.0, 4.0*t - 1.0, 1.0 - 4.0*u, 0.0,
                 4.0*s, 4.0*r, 4.0 * (u - t), -4.0*s, -4.0*r]]
        return deriv

    def PointInside(self,nc):
        u = 1.0 - nc[0] - nc[1] - nc[2]
        # jme; increased the tolerance to 0.00001 (from 0.001) to get closer
        # answers when implemented John Dailey's GeomUtils.pyd (11/9/04)
        if (nc[0] < -0.00001) or (nc[0] > 1.00001) or \
           (nc[1] < -0.00001) or (nc[1] > 1.00001) or \
           (nc[2] < -0.00001) or (nc[2] > 1.00001) or \
           (u < -0.00001) or (u > 1.00001): return 0
        return 1

    def NearestPoint(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        u = 1.0 - r - s - t
        if r < 0.0: r = 0.0
        if r > 1.0: r = 1.0
        if s < 0.0: s = 0.0
        if s > 1.0: s = 1.0
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        if u < 0.0:
            r = nc[0] - (1.0 - nc[0] - nc[1] - nc[2])/3.0
            s = nc[1] - (1.0 - nc[0] - nc[1] - nc[2])/3.0
            t = nc[2] - (1.0 - nc[0] - nc[1] - nc[2])/3.0
        return Vec3D.Vec3D(r,s,t)

    def CenterCoordinates(self):
        return Vec3D.Vec3D(0.25,0.25,0.25)

    def __repr__(self):
        return (self.nodes)

    def GiveSurf(self): #jme to build surface mesh
        surfaces = [[self.nodes[0], self.nodes[6], self.nodes[2],  \
                     self.nodes[5], self.nodes[1], self.nodes[4]], \
                    [self.nodes[0], self.nodes[9], self.nodes[3],  \
                     self.nodes[7], self.nodes[2], self.nodes[6]], \
                    [self.nodes[0], self.nodes[4], self.nodes[1],  \
                     self.nodes[8], self.nodes[3], self.nodes[9]], \
                    [self.nodes[1], self.nodes[5], self.nodes[2],  \
                     self.nodes[7], self.nodes[3], self.nodes[8]]]
        return surfaces



class Pyramid:

#
#                       
# 5 & 13-noded numbering               
#                                           
#              base                     mid             apex                   
#                         8                                           
#              0 o-------o-------o 3      12     11 
#               +               +          o-----o  
#              +               +          +     +           o 4
#           5 o               o 7        +     + 
#            +               +          o-----o 
#           +               +           9    10     
#        1 o-------o-------o 2
#                  6      
#
 
    def __init__(self,order,mat,nodes):
        self.order = order
        self.mat = mat
        if self.order == 1:
            self.num_node = 5
        else:
            self.num_node = 13
        assert len(nodes) >= self.num_node
        self.nodes = [nodes[i] for i in xrange(self.num_node)]

    def ShapeFunction(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        if self.order == 1:
            if t == 1.0:
                shape = [0.0,0.0,0.0,0.0,1.0]
            else:
                shape = [(1-(r/(1-t)))*(1-(s/(1-t)))*(1-t),
                         r*(1-(s/(1-t))),
                         r*s/(1-t),
                         (1-(r/(1-t)))*s,
                         t]
        else:
            raise 'No Quadratic Pyramids Yet'
        return shape
 
    def ShapeDerivatives(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        u = 1.0 - r - s - t
        if self.order == 1:
            if t == 1.0:
                deriv = [
                    [0.0,0.0,0.0,0.0,0.0],
                    [0.0,0.0,0.0,0.0,0.0],
                    [0.0,0.0,0.0,0.0,0.0]]
            else:
                deriv = [
                    [-(t+s-1)/(t-1), (t+s-1)/(t-1), -s/(t-1),
                     s/(t-1), 0.0],
                    [-(t+r-1)/(t-1), r/(t-1), -r/(t-1),
                     (t+r-1)/(t-1), 0.0],
                    [(r*s-1+2*t-t**2)/(t-1)**2,-r*s/(t-1)**2,
                     r*s/(t-1)**2,-r*s/(t-1)**2,1.0]]
        else:
            raise 'No Quadratic Pyramids Yet'
               
        return deriv

    def PointInside(self,nc):
        # jme; increased the tolerance to 0.00001 (from 0.001) to get closer
        # answers when implemented John Dailey's GeomUtils.pyd (11/9/04)
        r = nc[0]
        s = nc[1]
        t = nc[2]
        if self.order == 1:
            if t == 1.0:
                if r > 0.0 or s > 0.0: return 0
            if (r < -0.00001) or (r/(1-t) > 1.00001) or \
               (s < -0.00001) or (s/(1-t) > 1.00001) or \
               (t < -0.00001) or (t > 1.00001): return 0
            return 1
        else:
            raise 'No Quadratic Pyramids Yet'

    def NearestPoint(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        if t >= 1.0:
            r = 0.0
            s = 0.0
            t = 1.0
            if self.order == 1:
                if r < 0.0: r = 1.0
                if r/(1-t) > 1.0: r = 1.0 - t
                if s < 0.0: s = 0.0
                if s/(1-t) > 1.0: s = 1.0 - t
                if t < 0.0: t = 0.0
            else:
                raise 'No Quadratic Pyramids Yet'
        return Vec3D.Vec3D(r,s,t)

    def CenterCoordinates(self):
        if self.order == 1:
            return Vec3D.Vec3D(0.5,0.5,0.5)
        else:
            return Vec3D.Vec3D(0.0,0.0,0.5)

    def GiveSurf(self): #jme to build surface mesh
        surfaces = [[self.nodes[0], self.nodes[8], self.nodes[3],    \
                     self.nodes[7], self.nodes[2], self.nodes[6],    \
                     self.nodes[1], self.nodes[5]],                  \
                    [self.nodes[0], self.nodes[5], self.nodes[1],    \
                     self.nodes[9], self.nodes[4], self.nodes[12]],  \
                    [self.nodes[1], self.nodes[6], self.nodes[2],    \
                     self.nodes[10], self.nodes[4], self.nodes[9]],  \
                    [self.nodes[2], self.nodes[7], self.nodes[3],    \
                     self.nodes[11], self.nodes[4], self.nodes[10]], \
                    [self.nodes[3], self.nodes[8], self.nodes[0],    \
                     self.nodes[12], self.nodes[4], self.nodes[11]]]
        return surfaces
       


class Wedge:

#                                    5
#                                  o
#  5 &15-noded numbering         +/ \
#                            14 o/   \
#                             + /     \ 
#                          2+  o 11    o 10
#                          /\ /         \
#                         /  \     9     \
#                        /3 o-\----o------o 4 - - - - - t = +1
#                       /  +   \         +
#                    8 o +      o 7     +
#               r     / o 12     \     o 13
#               ^    / +          \   +
#               |   /+             \ +
#              --- o--------o-------o  - - - - - - - - t = -1
#                 / 0       6        1
#                /--> s

 
    def __init__(self,order,mat,nodes):
        self.order = order
        self.mat = mat
        if order == 1:
            self.num_node = 6
        else:
            self.num_node = 15
        assert len(nodes) >= self.num_node
        self.nodes = [nodes[i] for i in xrange(self.num_node)]

    def ShapeFunction(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        u = 1.0 - nc[0] - nc[1]
        if self.order == 1:
            shape = [r*(1-t), s*(1-t), u*(1-t),
                     r*t, s*t, u*t]
        else:
            shape = [0.5*u * ((2.0*u-1.0)*(1.0-t) - (1.0-t*t)),
                     0.5*s * ((2.0*s-1.0)*(1.0-t) - (1.0-t*t)),
                     0.5*r * ((2.0*r-1.0)*(1.0-t) - (1.0-t*t)),
                     0.5*u * ((2.0*u-1.0)*(1.0+t) - (1.0-t*t)),
                     0.5*s * ((2.0*s-1.0)*(1.0+t) - (1.0-t*t)),
                     0.5*r * ((2.0*r-1.0)*(1.0+t) - (1.0-t*t)),
                     2.0 * u * s * (1.0 - t),
                     2.0 * s * r * (1.0 - t),
                     2.0 * r * u * (1.0 - t),
                     2.0 * u * s * (1.0 + t),
                     2.0 * s * r * (1.0 + t),
                     2.0 * r * u * (1.0 + t),
                     u * (1.0 - t * t),
                     s * (1.0 - t * t),
                     r * (1.0 - t * t)]
        return shape

    def ShapeDerivatives(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        u = 1.0 - r - s
        if self.order == 1:
            deriv = [
                [ 1-t, 0.0, t-1, t, 0.0, -t],
                [ 0.0, 1-t, t-1, 0.0, t, -t],
                [ -r, -s, r+s-1, r, s, 1-r-s]]
        else:
            deriv = [
                [-0.5*(t-2.0+4.0*r+4.0*s)*(t-1.0),
                  0.0,
                  0.5*(t-4.0*r+2.0)*(t-1.0),
                 -0.5*(t+2.0-4.0*r-4.0*s)*(t+1.0),
                  0.0,
                  0.5*(t+4.0*r-2.0)*(t+1.0),
                  2.0*s*(t-1.0),
                 -2.0*s*(t-1.0),
                  2.0*(-1.0+2.0*r+s)*(t-1.0),
                 -2.0*s*(t+1.0),
                  2.0*s*(t+1.0),
                 -2.0*(-1.0+2.0*r+s)*(t+1.0),
                  (t-1.0)*(t+1.0),
                  0.0,
                 -(t-1.0)*(t+1.0)],
                [-0.5*(t-2.0+4.0*r+4.0*s)*(t-1.0),
                  0.5*(t-4.0*s+2)*(t-1.0),
                  0.0,
                 -0.5*(t+2.0-4.0*r-4.0*s)*(t+1.0),
                  0.5*(t+4.0*s-2.0)*(t+1.0),
                  0.0,
                  2.0*(-1.0+2.0*s+r)*(t-1.0),
                 -2.0*r*(t-1.0),
                  2.0*r*(t-1.0),
                 -2.0*(-1.0+2.0*s+r)*(t+1.0),
                  2.0*r*(t+1.0),
                 -2.0*r*(t+1.0),
                  (t-1.0)*(t+1.0),
                 -(t-1.0)*(t+1.0),
                  0.0],
                [-0.5*u*(2.0*u-1.0-2.0*t),
                 -0.5*s*(2.0*s-1.0-2.0*t),
                 -0.5*r*(2.0*r-1.0-2.0*t),
                  0.5*u*(2.0*u-1.0+2.0*t),
                  0.5*s*(2.0*s-1.0+2.0*t),
                  0.5*r*(2.0*r-1.0+2.0*t),
                 -2.0*u*s,
                 -2.0*s*r,
                 -2.0*r*u,
                  2.0*u*s,
                  2.0*s*r,
                  2.0*r*u,
                 -2.0*u*t,
                 -2.0*s*t,
                 -2.0*r*t]]
        return deriv

    def PointInside(self,nc):
        # jme; added a tolerance of 0.00001 to get closer answers when 
        # implemented John Dailey's GeomUtils.pyd (11/9/04)
        r = nc[0]
        s = nc[1]
        t = nc[2]
        u = 1.0 - r - s
        if self.order == 1:
            if (r < 0.00001) or (r > 1.00001) or \
               (s < 0.00001) or (s > 1.00001) or \
               (t < 0.00001) or (t > 1.00001) or \
               (u < 0.00001) or (u > 1.00001): return 0
            return 1
        else:
            if (r < 0.00001) or (r > 1.00001) or \
               (s < 0.00001) or (s > 1.00001) or \
               (t < -1.00001) or (t > 1.00001) or \
               (u < 0.00001) or (u > 1.00001): return 0
            return 1

    def NearestPoint(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        u = 1.0 - r - s
        if r < 0.0: r = 0.0
        if r > 1.0: r = 1.0
        if s < 0.0: s = 0.0
        if s > 1.0: s = 1.0
        if u < 0.0:
            r = (1.0 - nc[0] + nc[1]) / 2.0
            s = 1.0 - (1.0 - nc[0] + nc[1]) / 2.0
        if self.order == 1:
            if t < 0.0: t = 0.0
            if t > 1.0: t = 1.0
        else:
            if t < -1.0: t = -1.0
            if t > 1.0: t = 1.0
        return Vec3D.Vec3D(r,s,t)

    def CenterCoordinates(self):
        if self.order == 1:
            return Vec3D.Vec3D(0.333333,0.333333,0.5)
        else:
            return Vec3D.Vec3D(0.333333,0.333333,0.0)

    def GiveSurf(self): #jme to build surface mesh
        surfaces = [[self.nodes[0], self.nodes[12], self.nodes[3],  \
                     self.nodes[9], self.nodes[4], self.nodes[13],  \
                     self.nodes[1], self.nodes[6]], \
                    [self.nodes[1], self.nodes[13], self.nodes[4],  \
                     self.nodes[10], self.nodes[5], self.nodes[14], \
                     self.nodes[2], self.nodes[7]], \
                    [self.nodes[3], self.nodes[12], self.nodes[0],  \
                     self.nodes[8], self.nodes[2], self.nodes[14],  \
                     self.nodes[5], self.nodes[11]], \
                    [self.nodes[0], self.nodes[6], self.nodes[1],   \
                     self.nodes[7], self.nodes[2], self.nodes[8]],  \
                    [self.nodes[4], self.nodes[7], self.nodes[3],   \
                     self.nodes[11], self.nodes[5], self.nodes[10]]]
        return surfaces        


class Brick:

#
#                       2 o-------o-------o 6
# 8 & 20-noded numbering +|      18      +|
#                       + |             + |          s
#                    9 o  |         11 o  |          ^
#                     +   o 12        +   o 14       |
#                    +    | 19       +    |          |
#                 3 o-------o-------o 7   |          .-----> r
#                   |   0 o-------o-|-----o 4       /
#                   |    +       16 |    +         /
#                   |   +           |   +         t
#                13 o  o 8       15 o  o 10
#                   | +             | +
#                   |+              |+
#                 1 o-------o-------o 5
#                          17        
#

 
    def __init__(self,order,mat,nodes):
        self.order = order
        self.mat = mat
        if order == 1:
            self.num_node = 8
        else:
            self.num_node = 20
        assert len(nodes) >= self.num_node
        self.nodes = [nodes[i] for i in xrange(self.num_node)]

    def ShapeFunction(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        r2 = r**2
        s2 = s**2
        t2 = t**2
        if self.order == 1:
            shape = [(1.0-r)*(1.0-s)*(1.0-t),
                     (1.0-r)*(1.0-s)*t,
                     (1.0-r)*s*(1.0-t),
                     (1.0-r)*s*t,
                     r*(1.0-s)*(1.0-t),
                     r*(1.0-s)*t,
                     r*s*(1.0-t),
                     r*s*t]
        else:
            shape = [0.125*(1.0-r)*(1.0-s)*(1.0-t)*(-r-s-t-2.0),
                     0.125*(1.0-r)*(1.0-s)*(1.0+t)*(-r-s+t-2.0),
                     0.125*(1.0-r)*(1.0+s)*(1.0-t)*(-r+s-t-2.0),
                     0.125*(1.0-r)*(1.0+s)*(1.0+t)*(-r+s+t-2.0),
                     0.125*(1.0+r)*(1.0-s)*(1.0-t)*(r-s-t-2.0),
                     0.125*(1.0+r)*(1.0-s)*(1.0+t)*(r-s+t-2.0),
                     0.125*(1.0+r)*(1.0+s)*(1.0-t)*(r+s-t-2.0),
                     0.125*(1.0+r)*(1.0+s)*(1.0+t)*(r+s+t-2.0),
                     0.250*(1.0-r)*(1.0-s)*(1.0-t2),
                     0.250*(1.0-r)*(1.0+s)*(1.0-t2),
                     0.250*(1.0+r)*(1.0-s)*(1.0-t2),
                     0.250*(1.0+r)*(1.0+s)*(1.0-t2),
                     0.250*(1.0-r)*(1.0-s2)*(1.0-t),
                     0.250*(1.0-r)*(1.0-s2)*(1.0+t),
                     0.250*(1.0+r)*(1.0-s2)*(1.0-t),
                     0.250*(1.0+r)*(1.0-s2)*(1.0+t),
                     0.250*(1.0-r2)*(1.0-s)*(1.0-t),
                     0.250*(1.0-r2)*(1.0-s)*(1.0+t),
                     0.250*(1.0-r2)*(1.0+s)*(1.0-t),
                     0.250*(1.0-r2)*(1.0+s)*(1.0+t)]
        return shape

    def ShapeDerivatives(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        if self.order == 1:
            deriv = [
                [-(s-1.0)*(t-1.0),t*(s-1.0),s*(t-1.0),-s*t,
                 (s-1.0)*(t-1.0),-t*(s-1.0),-s*(t-1.0),s*t],
                [-(r-1.0)*(t-1.0),t*(r-1.0),(r-1.0)*(t-1.0),-t*(r-1.0),
                 r*(t-1.0),-r*t,-r*(t-1.0),r*t],
                [-(r-1.0)*(s-1.0),(r-1.0)*(s-1.0),s*(r-1.0),-s*(r-1.0),
                 r*(s-1.0),-r*(s-1.0),-r*s,r*s]]
        else:
            r2 = r*r
            s2 = s*s
            t2 = t*t
            deriv = [
                [-.125 * (1.0-s) * (1.0-t) * (-2.0*r-s-t-1.0),
                 -.125 * (1.0-s) * (1.0+t) * (-2.0*r-s+t-1.0),
                 -.125 * (1.0+s) * (1.0-t) * (-2.0*r+s-t-1.0),
                 -.125 * (1.0+s) * (1.0+t) * (-2.0*r+s+t-1.0),
                 0.125 * (1.0-s) * (1.0-t) * (2.0*r-s-t-1.0),
                 0.125 * (1.0-s) * (1.0+t) * (2.0*r-s+t-1.0),
                 0.125 * (1.0+s) * (1.0-t) * (2.0*r+s-t-1.0),
                 0.125 * (1.0+s) * (1.0+t) * (2.0*r+s+t-1.0),
                 -.250 * (1.0-s) * (1.0-t2),
                 -.250 * (1.0+s) * (1.0-t2),
                 0.250 * (1.0-s) * (1.0-t2),
                 0.250 * (1.0+s) * (1.0-t2),
                 -.250 * (1.0-s2) * (1.0-t),
                 -.250 * (1.0-s2) * (1.0+t),
                 0.250 * (1.0-s2) * (1.0-t),
                 0.250 * (1.0-s2) * (1.0+t),
                 -0.50 * r * (1.0-s) * (1.0-t),
                 -0.50 * r * (1.0-s) * (1.0+t),
                 -0.50 * r * (1.0+s) * (1.0-t),
                 -0.50 * r * (1.0+s) * (1.0+t)],
                [-.125 * (1.0-r) * (1.0-t) * (-r-2.0*s-t-1.0),
                 -.125 * (1.0-r) * (1.0+t) * (-r-2.0*s+t-1.0),
                 0.125 * (1.0-r) * (1.0-t) * (-r+2.0*s-t-1.0),
                 0.125 * (1.0-r) * (1.0+t) * (-r+2.0*s+t-1.0),
                 -.125 * (1.0+r) * (1.0-t) * (r-2.0*s-t-1.0),
                 -.125 * (1.0+r) * (1.0+t) * (r-2.0*s+t-1.0),
                 0.125 * (1.0+r) * (1.0-t) * (r+2.0*s-t-1.0),
                 0.125 * (1.0+r) * (1.0+t) * (r+2.0*s+t-1.0),
                 -.250 * (1.0-r) * (1-t2),
                 0.250 * (1.0-r) * (1.0-t2),
                 -.250 * (1.0+r) * (1.0-t2),
                 0.250 * (1.0+r) * (1.0-t2),
                 -.500 * (1.0-r) * s * (1.0-t),
                 -.500 * (1.0-r) * s * (1.0+t),
                 -.500 * (1.0+r) * s * (1.0-t),
                 -.500 * (1.0+r) * s * (1.0+t),
                 -.250 * (1.0-r2) * (1.0-t),
                 -.250 * (1.0-r2) * (1.0+t),
                 0.250 * (1.0-r2) * (1.0-t),
                 0.250 * (1.0-r2) * (1.0+t)],
                [-.125 * (1.0-r) * (1.0-s) * (-r-s-2.0*t-1.0),
                 0.125 * (1.0-r) * (1.0-s) * (-r-s+2.0*t-1.0),
                 -.125 * (1.0-r) * (1.0+s) * (-r+s-2.0*t-1.0),
                 0.125 * (1.0-r) * (1.0+s) * (-r+s+2.0*t-1.0),
                 -.125 * (1.0+r) * (1.0-s) * (r-s-2.0*t-1.0),
                 0.125 * (1.0+r) * (1.0-s) * (r-s+2.0*t-1.0),
                 -.125 * (1.0+r) * (1.0+s) * (r+s-2.0*t-1.0),
                 0.125 * (1.0+r) * (1.0+s) * (r+s+2.0*t-1.0),
                 -.500 * (1.0-r) * (1.0-s) * t,
                 -.500 * (1.0-r) * (1.0+s) * t,
                 -.500 * (1.0+r) * (1.0-s) * t,
                 -.500 * (1.0+r) * (1.0+s) * t,
                 -.250 * (1.0-r) * (1.0-s2),
                 0.250 * (1.0-r) * (1.0-s2),
                 -.250 * (1.0+r) * (1.0-s2),
                 0.250 * (1.0+r) * (1.0-s2),
                 -.250 * (1.0-r2) * (1.0-s),
                 0.250 * (1.0-r2) * (1.0-s),
                 -.250 * (1.0-r2) * (1.0+s),
                 0.250 * (1.0-r2) * (1.0+s)]]
        return deriv

    def PointInside(self,nc):
        # jme; increased the tolerance to 0.00001 (from 0.001) to get closer
        # answers when implemented John Dailey's GeomUtils.pyd (11/9/04)
        r = nc[0]
        s = nc[1]
        t = nc[2]
        if self.order == 1:
            if (r < -0.00001) or (r > 1.00001) or \
               (s < -0.00001) or (s > 1.00001) or \
               (t < -0.00001) or (t > 1.00001): return 0
            return 1
        else:
            if (r < -1.00001) or (r > 1.00001) or \
               (s < -1.00001) or (s > 1.00001) or \
               (t < -1.00001) or (t > 1.00001): return 0
            return 1

    def NearestPoint(self,nc):
        r = nc[0]
        s = nc[1]
        t = nc[2]
        if self.order == 1:
            if r < 0.0: r = 0.0
            if r > 1.0: r = 1.0
            if s < 0.0: s = 0.0
            if s > 1.0: s = 1.0
            if t < 0.0: t = 0.0
            if t > 1.0: t = 1.0
        else:
            if r < -1.0: r = -1.0
            if r >  1.0: r =  1.0
            if s < -1.0: s = -1.0
            if s >  1.0: s =  1.0
            if t < -1.0: t = -1.0
            if t >  1.0: t =  1.0
        return Vec3D.Vec3D(r,s,t)

    def CenterCoordinates(self):
        if self.order == 1:
            return Vec3D.Vec3D(0.5,0.5,0.5)
        else:
            return Vec3D.Vec3D(0.0,0.0,0.0)

    def GiveSurf(self): #jme to build surface mesh
        surfaces = [[self.nodes[0], self.nodes[16], self.nodes[4], \
                     self.nodes[10], self.nodes[5], self.nodes[17],\
                     self.nodes[1], self.nodes[8]] \
                    [self.nodes[0], self.nodes[8], self.nodes[1], \
                     self.nodes[13], self.nodes[3], self.nodes[9],\
                     self.nodes[2], self.nodes[12]] \
                    [self.nodes[1], self.nodes[17], self.nodes[5], \
                     self.nodes[15], self.nodes[7], self.nodes[19],\
                     self.nodes[3], self.nodes[13]] \
                    [self.nodes[5], self.nodes[10], self.nodes[4], \
                     self.nodes[14], self.nodes[6], self.nodes[11],\
                     self.nodes[7], self.nodes[15]] \
                    [self.nodes[4], self.nodes[16], self.nodes[0], \
                     self.nodes[12], self.nodes[2], self.nodes[18],\
                     self.nodes[6], self.nodes[14]] \
                    [self.nodes[2], self.nodes[9], self.nodes[3], \
                     self.nodes[19], self.nodes[7], self.nodes[11],\
                     self.nodes[6], self.nodes[18]]]
        return surfaces        



#------------------------------------------------------
# Test Driver
#------------------------------------------------------

if __name__ == "__main__":

    def Gindx(i,j,k):
        return i*16 + j*4 + k


    model = FemModel()

    # generate nodes for a 3x3x3 element block

    cur = 0 ;
    for i in xrange(4):
        for j in xrange(4):
            for k in xrange(4):
                coord = Vec3D.Vec3D(float(i),float(j),float(k))
                disp = Vec3D.Vec3D(float(i),float(j),float(k))
                stress = ColTensor.ColTensor(float(i),float(j),
                                             float(k),0,0,0)
                model.AddNode(cur,coord,disp,stress)
                cur += 1

    # create elements

    cur = 0 ;
    for i in xrange(3):
        for j in xrange(3):
            for k in xrange(3):
                nodes = [Gindx(i  ,j  ,k),Gindx(i  ,j  ,k+1),
                         Gindx(i  ,j+1,k),Gindx(i  ,j+1,k+1),
                         Gindx(i+1,j  ,k),Gindx(i+1,j  ,k+1),
                         Gindx(i+1,j+1,k),Gindx(i+1,j+1,k+1)]
                model.AddBrick(cur,1,1,nodes)
                cur += 1

    # query for the stresses at a point (x,y,z normal stress
    # values should be the same as the coordinates of the
    # query point)

    print model.GetPtStress(Vec3D.Vec3D(0.75,1.25,1.25))

