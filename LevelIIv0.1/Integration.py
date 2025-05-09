# a module to do numerical integration by various methods
import math

def star(vec,fac): return [vec[i]*fac for i in range(len(vec))]
def plus(v1,v2): return [v1[i]+v2[i] for i in range(len(v1))]
def minus(v1,v2): return [v1[i]-v2[i] for i in range(len(v1))]
def norm(vec):
    summation = 0.0
    for i in range(len(vec)): summation += vec[i]**2
    return summation**.5
def maxnorm(vec): return (max([vec[i]**2 for i in range(len(vec))]))**.5

def Euler(step,current_independant,current_dependant,function,args):
    '''
    Takes a single Euler integration step and returns the incremented value
    of the dependant variable
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
          current_independant- unincremented value of the independant
                               variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
                               'function(independant,dependant)'
    '''
    jump = step*function(current_independant,current_dependant,args)
    return current_dependant+jump

def Eulerslope_vector(step,current_independant,current_dependant,function, \
                      args):
    '''
    Takes a single Euler integration step and returns the slope of the
    function for that evaluation
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
          current_independant- unincremented value of the independant
                               variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
                               'function(independant,dependant)'
    '''
    return function(current_independant,current_dependant,args)

def MidPoint(step,current_independant,current_dependant,function,args):
    '''
    Takes a single midpoint integration step and returns the incremented
    value of the dependant variable
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
          current_independant- unincremented value of the independant
                               variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
                               'function(independant,dependant)'
    '''
    k1=step*function(current_independant,current_dependant,args)
    k2=step*function(current_independant+.5*step,current_dependant+.5*k1,args)
    return current_dependant+k2

def RK4(step,current_independant,current_dependant,function,args):
    '''
    Takes a single 4th order Runge-Kutta integration step and returns 
      the incremented value of the dependant variable
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
          current_independant- unincremented value of the independant
                               variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
                               'function(independant,dependant)'
    '''
    k1=step*function(current_independant,current_dependant,args)
    k2=step*function(current_independant+.5*step,current_dependant+.5*k1,args)
    k3=step*function(current_independant+.5*step,current_dependant+.5*k2,args)
    k4=step*function(current_independant+step   ,current_dependant+k3,args   )
    jump = k1/6 + k2/3 + k3/3 + k4/6
    return current_dependant + jump

def RK4slope(step,current_independant,current_dependant,function,args):
    '''
    Takes a single 4th order Runge-Kutta integration step and returns 
    the incremented value of the dependant variable
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
          current_independant- unincremented value of the independant
                               variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
                               'function(independant,dependant)'
    '''
    k1=step*function(current_independant,current_dependant,args)
    k2=step*function(current_independant+.5*step,current_dependant+.5*k1,args)
    k3=step*function(current_independant+.5*step,current_dependant+.5*k2,args)
    k4=step*function(current_independant+step   ,current_dependant+k3,args   )
    jump = k1/6 + k2/3 + k3/3 + k4/6
    return jump/step

def RK4slope_vector(step,current_independant,current_dependant,function,args):
    '''
    Takes a single 4th order Runge-Kutta integration step and returns 
    the incremented value of the dependant variable
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
                               vector
          current_independant- unincremented value of the independant
                               variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
    #                            'function(independant,dependant)'
    '''
    k1 = star(function(current_independant,current_dependant,args),step)
    k2 = star(function(current_independant+.5*step, \
                       plus(star(k1,.5),current_dependant),args), \
              step)
    k3 = star(function(current_independant+.5*step,\
                       plus(star(k2,.5),current_dependant),args), \
              step)
    k4 = star(function(current_independant+step, \
                       plus(k3,current_dependant),args), \
              step)
    jump = plus(star(k1,1.0/6),plus(star(plus(k2,k3),1.0/3),star(k4,1.0/6)))
    #print 'regular',step,k1
    return star(jump,1.0/step)

def RK4_doubling(step,current_independant,current_dependant,function,args, \
                 target_error):
    '''
    Takes a single 4th order Runge-Kutta integration step with step doubling
    adaptivity and returns the incremented value of the dependant variable,
    the size of the step taken, and the predicted size of the next step
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
          current_independant- unincremented value of the independant
          variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
                               'function(independant,dependant)'
          target_error       - the allowable truncation error
    '''
    step_2= step/2
    k1_1  = step_2*function(current_independant,current_dependant,args)
    k2_1  = step_2*function(current_independant+.25*step, \
                            current_dependant+.5*k1_1,args)
    k3_1  = step_2*function(current_independant+.25*step, \
                            current_dependant+.5*k2_1,args)
    k4_1  = step_2*function(current_independant+step_2  , \
                            current_dependant+k3_1,args)
    indep = current_independant + step_2
    jump1 = k1_1/6 + k2_1/3 + k3_1/3 + k4_1/6
    dep_1 = current_dependant + jump1
    k1_2  = step_2*function(indep, dep_1,args)
    k2_2  = step_2*function(indep+.25*step, dep_1+.5*k1_2,args)
    k3_2  = step_2*function(indep+.25*step, dep_1+.5*k2_2,args)
    k4_2  = step_2*function(indep+step_2  , dep_1+k3_2,args)
    jump2 = k1_2/6 + k2_2/3 + k3_2/3 + k4_2/6
    dep_2 = dep_1 + jump2
    k1    = k1_1*2.0
    k2    = step*function(current_independant+step_2,current_dependant+.5*k1, \
                          args)
    k3    = step*function(current_independant+step_2,current_dependant+.5*k2, \
                          args)
    k4    = step*function(current_independant+step  ,current_dependant+k3, \
                          args)
    dep   = current_dependant + k1/6 + k2/3 + k3/3 + k4/6
    error = math.fabs(dep_2-dep)
    ret_step=step
    if (error>target_error):
        return RK4_doubling(step_2,current_independant,current_dependant, \
                            function,args,target_error)
    if (error<target_error/2):
        ret_step = step*2
    return dep_2,step,ret_step

def RK4slope_vector_doubling(step,current_independant,current_dependant, \
                             function,args,target_error,errnorm=maxnorm):
    '''
    Takes a single 4th order Runge-Kutta integration step with step doubling
    adaptivity and returns the incremented value of the dependant variable,
    the size of the step taken, and the predicted size of the next step
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
                               vector
          current_independant- unincremented value of the independant
                               variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
                               'function(independant,dependant)'
      target_error       - the allowable truncation error
    '''
    step_2= step/2.0
    k1_1 = star(function(current_independant,current_dependant,args),step_2)
    k2_1 = star(function(current_independant+.25*step, \
                         plus(star(k1_1,.5),current_dependant),args),step_2)
    k3_1 = star(function(current_independant+.25*step, \
                         plus(star(k2_1,.5),current_dependant),args),step_2)
    k4_1 = star(function(current_independant+step_2, \
                         plus(k3_1,current_dependant),args),step_2)
    jump1 = plus(star(k1_1,1.0/6), \
                 plus(star(plus(k2_1,k3_1),1.0/3),star(k4_1,1.0/6)))
    indep = current_independant+step_2
    dep_1 = plus(current_dependant,jump1)
    k1_2 = star(function(indep,dep_1,args),step_2)
    k2_2 = star(function(indep+.25*step,plus(star(k1_2,.5),dep_1),args),step_2)
    k3_2 = star(function(indep+.25*step,plus(star(k2_2,.5),dep_1),args),step_2)
    k4_2 = star(function(indep+step_2,plus(k3_2,dep_1),args),step_2)
    jump2 = plus(star(k1_2,1.0/6), \
                 plus(star(plus(k2_2,k3_2),1.0/3),star(k4_2,1.0/6)))
    jumpT = plus(jump1,jump2)
    k1 = star(k1_1,2.0)
    k2 = star(function(current_independant+step_2, \
                       plus(star(k1,.5),current_dependant),args),step)
    k3 = star(function(current_independant+step_2, \
                       plus(star(k2,.5),current_dependant),args),step)
    k4 = star(function(current_independant+step, \
                       plus(k3,current_dependant),args),step)
    jump = plus(star(k1,1.0/6),plus(star(plus(k2,k3),1.0/3),star(k4,1.0/6)))
    error = errnorm(minus(jumpT,jump))
    #print 'doubling',step,k1
    ret_step = step
    if (error>target_error):
        return RK4slope_vector_doubling(step_2,current_independant, \
                                        current_dependant,function,args, \
                                        target_error)
    if (error<target_error/2):
        ret_step = step*2
    return star(jumpT,1.0/step),step,ret_step

def RK5_CK(step,current_independant,current_dependant,function,args, \
           max_error,target_error=0.9,recurse=1):
    '''
    Takes a single 5th order Runge-Kutta integration step with Cash Karp
    adaptivity and returns the incremented value of the dependant variable,
    the size of the step taken, and the predicted size of the next step
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
          current_independant- unincremented value of the independant
                               variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
                               'function(independant,dependant)'
          max_error          - the max error allowed before recursive
                               refinement
          target_error       - (optional) the target truncation error
                               expressed as a multiple of max_error, should
                               be <= 1.
          recurse            - (optional) turn recursion on/off for cases of
                               larger error
    '''
    k1 = step*function(current_independant,current_dependant,args)
    k2 = step*function(current_independant+step/5   ,\
                       current_dependant + k1/5.0,args)
    k3 = step*function(current_independant+3*step/10,\
                       current_dependant + 3*k1/40 + 9*k2/40,args)
    k4 = step*function(current_independant+3*step/5 ,\
                       current_dependant + 3*k1/10 - 9*k2/10 + 6*k3/5,args)
    k5 = step*function(current_independant+  step   ,\
             current_dependant - 11*k1/54 + 5*k2/2 - 70*k3/27 + 35*k4/27,args)
    summ =1631.*k1/55296. + 175.*k2/512.
    summ+=575.*k3/13824.  + 44275.*k4/110592.
    summ+=253.*k5/4096.
    k6 = step*function(current_independant+7*step/8 , \
                       current_dependant+summ,args)
    jump = 37*k1/378 + 250*k3/621 + 125*k4/594 + 512*k6/1771
    dep = current_dependant + jump
    error = (37.0/378-2825.0/27648)*k1 + (250.0/621-18575.0/48384)*k3 + \
            (125.0/584-13525.0/55296)*k4 - 277.0/14336*k5 + (512.0/1771-.25)*k6
    error = max((max_error*target_error/100.0),error)
    err_ratio = math.fabs(max_error*target_error/error)
    if (recurse and (err_ratio/target_error < 1.0)):
        return RK5_CK(step*err_ratio**.25,current_independant, \
                      current_dependant,function,args,max_error,target_error)
    else:
        return dep,step,step*err_ratio**.2

def RK5_CKslope_vector(step,current_independant,current_dependant,function, \
                       args,max_error,target_error=0.75,recurse=1, \
                       errnorm=maxnorm):
    '''
    Takes a single 5th order Runge-Kutta integration step with Cash Karp
    adaptivity and returns the incremented value of the dependant variable,
    the size of the step taken, and the predicted size of the next step
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
                               vector
          current_independant- unincremented value of the independant
          variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
                               'function(independant,dependant)'
          max_error          - the max error allowed before recursive
                               refinement
          target_error       - (optional) the target truncation error
                                expressed as a multiple of max_error, should
                               be <= 1.
          recurse            - (optional) turn recursion on/off for cases of
                               larger error
          errnorm            - (optional) the function to use in converting  
                               the error vector to a scalar
    '''
##    print 'k1', current_independant,current_dependant,args
    k1 = star(function(current_independant,current_dependant,args),step)
##    print 'k2', current_independant+step/5.0, \
##                       plus(current_dependant,star(k1,.2)),args
    k2 = star(function(current_independant+step/5.0, \
                       plus(current_dependant,star(k1,.2)),args),step)
##    print 'k3', current_independant+3.0*step/10, \
##                       plus(current_dependant, \
##                       plus(star(k1,3.0/40),star(k2,9.0/40))),args
    k3 = star(function(current_independant+3.0*step/10, \
                       plus(current_dependant, \
                       plus(star(k1,3.0/40),star(k2,9.0/40))),args),step)
##    print 'k4', current_independant+3.0*step/5, \
##                       plus(current_dependant, \
##                       plus(star(k1,.3),plus(star(k2,-.9),star(k3,1.2)))), \
##                       args
    k4 = star(function(current_independant+3.0*step/5, \
                       plus(current_dependant, \
                       plus(star(k1,.3),plus(star(k2,-.9),star(k3,1.2)))), \
                       args),step)
##    print 'k5', current_independant+step,\
##                       plus(current_dependant,plus(star(k1,-11.0/54),\
##                       plus(star(k2,2.5),plus(star(k3,-70.0/27), \
##                       star(k4,35.0/27))))),args
    k5 = star(function(current_independant+step,\
                       plus(current_dependant,plus(star(k1,-11.0/54),\
                       plus(star(k2,2.5),plus(star(k3,-70.0/27), \
                       star(k4,35.0/27))))),args),step)
##    print 'k6', current_independant+7.0*step/8,\
##                       plus(current_dependant, \
##                       plus(star(k1,1631.0/55296),plus(star(k2,175.0/512),\
##                       plus(star(k3,575.0/13824), \
##                       plus(star(k4,44275.0/110592),star(k5,253.0/4096)))))), \
##                       args
    k6 = star(function(current_independant+7.0*step/8,\
                       plus(current_dependant, \
                       plus(star(k1,1631.0/55296),plus(star(k2,175.0/512),\
                       plus(star(k3,575.0/13824), \
                       plus(star(k4,44275.0/110592),star(k5,253.0/4096)))))), \
                       args),step)
    jump = plus(star(k1,37.0/378),plus(star(k3,250.0/621), \
                                plus(star(k4,125.0/594),star(k6,512.0/1771))))
    err = plus(star(k1,(37.0/378-2825.0/27648)), \
               plus(star(k3,(250.0/621-18575.0/48384)), \
                    plus(star(k4,(125.0/584-13525.0/55296)),\
            plus(star(k5,277.0/14336),star(k6,(512.0/1771-.25))))))
    error = errnorm(err)
    error = max((max_error*target_error/100.0),error)
    err_ratio = math.fabs(max_error*target_error/error)
    if (recurse and (err_ratio/target_error < 1.0)):
        #print "Refining"
        return RK5_CKslope_vector(step*err_ratio**.25,current_independant, \
                                  current_dependant,function,args,max_error,\
                                  target_error,recurse,errnorm)
    else:
        return star(jump,1.0/step),step,step*err_ratio**.2

def AB3(step,current_dependant,function_eval,func_minus1,func_minus2):
    '''
    Takes a single 3rd order adams-bashforth step and returns the
    incremented value of the dependant variable
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
          function_eval      - the value of the funciton being integrated at
                               the current step
          func_minus1        - the value of the function on the previous step
          func_minus2        - the value of the function 2 steps previous
    '''
    jump = step/12.0*(23*function_eval - 16*func_minus1 + 5*func_minus2)
    return current_dependant+jump

def AM3(step,current_dependant,function_eval,func_plus1,func_minus1):
    '''
    Takes a single 3rd order adams-moulton step and returns the incremented
    value of the dependant variable
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
          function_eval      - the value of the funciton being integrated at
                               the current step
          func_plus1         - guess of the value of the function at the next
                               step
          func_minus2        - the value of the function 2 steps previous
    '''
    jump = step/12.0*(5*func_plus1 + 8*function_eval - func_minus1)
    return current_dependant+jump

def PEC_AM3(step,current_independant,current_dependant,function,func_minus1, \
            func_minus2):
    '''
    Takes a single 3rd order adams-bashforth-moulton prediction-evaluation-
    correction step and returns the incremented value of the dependant variable
    args: step               - step size to take
          current_dependant  - unincremented value of the dependant variable
          current_independant- unincremented value of the independant variable
          function           - the right hand side of the diff eq, should be
                               callable in the form
                               'function(independant,dependant)'
          func_minus1        - the value of the function on the previous step
          func_minus2        - the value of the function 2 steps previous
    '''
    eval = function(current_independant,current_dependant)
    yp_guess = AB3(step,current_dependant,eval,func_minus1,func_minus2)
    evalp = function(yp_guess,current_dependant+step)
    return AM3(step,current_dependant,eval,evalp,func_minus1)
    
