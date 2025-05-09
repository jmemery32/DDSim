#!/usr/local/bin/python

# This is from Wash originally! 

import sys
import math
import Numeric
import LinearAlgebra

# Exceptions:
##class NewtonExceptions(Exception):
##    '''
##    class of exceptions for this module
##    '''
##    def __init__(self):
##        NewtonSolve_ExceedMaxIts = "Newton.Solve() no convergence"

def PrintResid(res):
    for i in xrange(len(res)/3):
        print " %d  %11.4g %11.4g %11.4g" % \
             (i,res[i*3],res[i*3+1],res[i*3+2])

def DumpDisplacements(fem_data,disp):

    node_ids = fem_data.model.GetNodeIds()
    print "Nodal Displacements:"
    for node in node_ids:
        eqs = fem_data.eqn_nums[node]
        print " %d  %22.15e %22.15e %22.15e" % \
                (node,disp[eqs[0]],disp[eqs[1]],disp[eqs[2]])
        #print " %d  %11.4g %11.4g %11.4g" % \
        #        (node,disp[eqs[0]],disp[eqs[1]],disp[eqs[2]])

def DumpResidual(fem_data,res):

    node_ids = fem_data.model.GetNodeIds()
    print "Nodal Residuals:"
    for node in node_ids:
        eqs = fem_data.eqn_nums[node]
        print " %d  %22.15e %22.15e %22.15e" % \
                (node,res[eqs[0]],res[eqs[1]],res[eqs[2]])

def Solve(num_eqns,X,tol,EvalFunc,TangFunc,cdata,maximum_iterations=None):

    GlobalCounter = 0

##    print "X:", X

    if maximum_iterations: MAX_ITER = maximum_iterations
    else: MAX_ITER = 100

    # get memory for the solve
    del_X_itr = Numeric.zeros((num_eqns),Numeric.Float64)
    del_X_cum = Numeric.zeros((num_eqns),Numeric.Float64)
    func_val = Numeric.zeros((num_eqns),Numeric.Float64)
    jacobian = Numeric.zeros((num_eqns,num_eqns),Numeric.Float64)

    # evaluate the nonlinear function and compute the function norm
    EvalFunc(cdata,del_X_cum,X,func_val) # only change is in func_val.

    #DumpResidual(cdata,func_val)

    # if initial guess is right on, return it and be done!  
    norm = 0.0
    for i in xrange(num_eqns):
        norm += func_val[i]**2
    norm = math.sqrt(norm)
    if norm <= tol: return X

    # loop until we converge
    iter = 0
    while 1:

        iter += 1

        # solve to find the newton step

        TangFunc(cdata,del_X_cum,X,jacobian,GlobalCounter) #jacobian changed.

        #DumpResidual(cdata,func_val)

        #print "diag: "
        #for i in xrange(len(func_val)):
        #    print i,jacobian[i,i]

        # dx
        del_X_iter = LinearAlgebra.solve_linear_equations(jacobian,func_val)
        del_X_cum -= del_X_iter # minus because Dx = - F(x)/F'(x)
##        print "del_X_iter:", del_X_iter

        #DumpDisplacements(cdata,del_X_cum)
        GlobalCounter += 1

        # Be careful to sum X + del_X_xum inside EvalFunc()... 
        EvalFunc(cdata,del_X_cum,X,func_val)

        #DumpResidual(cdata,func_val)

        # Calc. norm
        norm = 0.0
        for i in xrange(num_eqns):
            norm += func_val[i]**2
        norm = math.sqrt(norm)

##        print "norm,tol : ",norm,tol
##        print "func_val : ", func_val
##        print "del_X_cum: ", del_X_cum,"\n"

        # Check against tol, if meets criteria, exit while loop
        if norm <= tol:
            # return the solution
            del_X_cum += X
            return del_X_cum

        # if exceeds iterations, raise exception and pass current estimate
        if iter > MAX_ITER and norm > tol:
            raise Exception("Newton.Solve() no convergence", del_X_cum + X)

    
# -------------------------------------------------------------


if __name__ == '__main__':

    # for a test case we try to find the inverse mapping
    # in a linear quad element

    class InverseMap:

        def __init__(self,coords,x,y):
            self.__coords = [(xt[0],xt[1]) for xt in coords]
            self.__x = x
            self.__y = y

        def EvalFunc(self,cdata,del_X,X,func_val):
            u = [X[0]+del_X[0],X[1]+del_X[1]]
            N = [(1-u[0])*(1-u[1]),u[0]*(1-u[1]),u[0]*u[1],(1-u[0])*u[1]]
            xt = N[0]*self.__coords[0][0] + \
                 N[1]*self.__coords[1][0] + \
                 N[2]*self.__coords[2][0] + \
                 N[3]*self.__coords[3][0]
            yt = N[0]*self.__coords[0][1] + \
                 N[1]*self.__coords[1][1] + \
                 N[2]*self.__coords[2][1] + \
                 N[3]*self.__coords[3][1]
            func_val[0] = self.__x - xt
            func_val[1] = self.__y - yt
##            print u[0],u[1],func_val[0],func_val[1]


        def TangFunc(self,cdata,del_X,X,jacobian,globe_counter):
            '''
            globe_counter is not used in this simple example.
            '''
            u = [X[0]+del_X[0],X[1]+del_X[1]]
            dNdr = [-(1-u[1]),1-u[1],u[1],-u[1]]
            dNds = [-(1-u[0]),-u[0],u[0],1-u[0]]
            jacobian[0,0] = dNdr[0]*self.__coords[0][0] + \
                            dNdr[1]*self.__coords[1][0] + \
                            dNdr[2]*self.__coords[2][0] + \
                            dNdr[3]*self.__coords[3][0]
            jacobian[0,1] = dNdr[0]*self.__coords[0][1] + \
                            dNdr[1]*self.__coords[1][1] + \
                            dNdr[2]*self.__coords[2][1] + \
                            dNdr[3]*self.__coords[3][1]
            jacobian[1,0] = dNds[0]*self.__coords[0][0] + \
                            dNds[1]*self.__coords[1][0] + \
                            dNds[2]*self.__coords[2][0] + \
                            dNds[3]*self.__coords[3][0]
            jacobian[1,1] = dNds[0]*self.__coords[0][1] + \
                            dNds[1]*self.__coords[1][1] + \
                            dNds[2]*self.__coords[2][1] + \
                            dNds[3]*self.__coords[3][1]

    # set up the test case

    coords = ((0,0),(1,0),(1,1),(0,1))
    pt = (0.25,0.25)
    im = InverseMap(coords,pt[0],pt[1])

    # find a tolerance based on the diagonal lengths

    dx = coords[2][0] - coords[0][0]
    dy = coords[2][1] - coords[0][1]
    dlen = math.sqrt(dx**2 + dy**2)

    dx = coords[3][0] - coords[1][0]
    dy = coords[3][1] - coords[1][1]
    dlen += math.sqrt(dx**2 + dy**2)

    tol = 0.0001 * (dlen/2)

    # solve

    delta = Solve(2,(0.5,0.5),tol,im.EvalFunc,im.TangFunc,None)

    print delta

