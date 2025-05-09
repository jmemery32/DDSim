import sys, time, string
import MeshTools
import DamMo as DamModel
import VarAmplitude
import Parameters
import Vec3D, ColTensor 
import numpy as Numeric
import math, os, cPickle, random
import Statistic # for plotting cdf of ai's.

# verification flag --> makes some things print automatically
#                      (1, true, makes it print)
verification = 1

def HelpPrints():
    '''
    function to print some help for the various switches.
    ''' 

    print ' '
    print '         DDSim command line arguments.'
    print ' The various command line switches are...'
    print ' -base -> follow by base file name.  Overrides further prompting'
    print '          for file name.'
    print ' -conpath -> if specify -base, include -conpath if necessary (i.e. '
    print '          if input files are not in current working directory.) One'
    print '          would never provide -conpath but not -base.  However, '
    print '          -base may be specified without -conpath'
    print '          eg.  ...>DDSim.py -base <filename> -conpath <path>'
    print '          eg.  ...>DDSim.py -base example1   -conpath examples\\'
    print ' -parpath -> same usage as above, location of *.par file'
    print ' -eg   -> process example1 (no further input necessary) '
    print ' -o    -> to be prompted for output options'
    print " -pa   -> means we just want it to print all results."
    print " -pdid -> print the damage elements (DamHistory)."
    print " -pdoid -> print the damage origin results."
    print ' -v    -> print damage results as they process'
    # -sv doesn't do anything right now... output.n and.sif are save
    # automatically with -DB option... 
##    print ' -sv   -> save output.N, output.SIF text files'
    print ' -sp   -> save pickled N file and .ori, .ai, .af files'
    print ' -sc   -> save files for viewing contours w/ MAP (*.MP)'
    print ' -si   -> save intermediate files of N for each doid as processed'
    print '         (not for parallel jobs, must include subdirectory \\inter)'
    print ' -L    -> limit the simulation spatially (will prompt for limits)'
    print ' -S    -> limit the simulation to surface nodes only'
    print ' -p    -> parallel job (include jobname.bat file in job directory)'
    print ' -usa  -> means use a user supplied ais by *.rnd file (binary no '
    print '          longer supported!)'
    print ' -DB   -> means a parallel run from SQL where the random initial '
    print '          crack sizes are supplied by the database'
    print ' -doid_list -> to specify node ids to run ddsim for.  Use #,# '
    print '               i.e. comma, no spaces! '
    print ' -pick_list -> to specify a pickled list of doids.  Follow key with'
    print '               filename.extentions'
    print ' '
    print " Integrations: Default is Constant amplitude, Adaptive time step,"
    print " RK-5" 
    print ' -Fwd  --> Constant amplitude, forward euler'
    print ' -RK4  --> Constand amplitude, 4 point RK scheme'
    print " -Simp --> Constant amplitude, Simpson's rule on generated SIF "
    print "           history"
    print " -VarAmp -> Variable amplitude loading, cycle-by-cycle integration"
    print "            Should be followed by the amplitude history file name. "
    print ' '
    print " -scale -> A multiplier for K values to scale stress if necessary."
    print "           Should be followed with an float greater than 0."
    print " -seed -> a seed to generate the monte samples with. "
    print "          Should be followed with a number."
    print " -nore -> specify to prevent the use of the Willenborg retardation"
    print "          model. i.e. growth rate is given by NASGRO eqn w/o "
    print "          retardation effects."
    print " -verify -> to print diagnostic info to the screen."
    print " -DebugGeomUtils -> to make create a file that our script to test "
    print "                  GeomUtils can read.  Follow key (-DebugGeomUtils)"
    print "                  by string to indicate name and location of file."
    print " -SVIEW <filename> -> to write an sview file.  Follow key with "
    print "           string to indicate name and location of sview file. "
    print " -cf -> prints crack front points to a file named <filename>.front"
    print "        ONLY use this for single doid runs."
    

#############

def LocalNodeInitializeStuff(model,extension,parameters):
    node_list = model.GetNodeList()

    # node_partition.%MSTI_RANK% should be there... 
    tfile=open('node_partition.'+extension,'r')
    doid_list=[]
    lines=tfile.readlines()
    for line in lines:
        doid_list+=[int(line)]
    tfile.close()

    # input.rnd should be here.  map for rid -> ai
    tfile=open("input.rnd",'r')
    ais=Statistic.Stata(parameters.sets,parameters.samples)
    for i in range(parameters.sets):
        ais.UpdateSet()
        for j in range(parameters.samples):
            l=tfile.readline()
            l=string.splitfields(l)
            ais.UpdateSamples(float(l[1]),i,int(l[0]))
    tfile.close()

    ais_map = None # remains None unless monte = 2
    if parameters.monte == 2: 
        # input.map should be here too. creates the map for nid -> rid
        tfile=open("input.map",'r')
        ais_map={}
        tfilelist=tfile.readlines()
        for line in tfilelist:
            l = string.splitfields(line)
            ais_map[int(l[0])]=map(int,l[2:])
        tfile.close()

    ncr=Numeric.arange(0,parameters.N_max*1.01,1000)
    ncr=ncr.tolist()

    return doid_list,node_list,ais,ncr,ais_map

#############

def ParallelInitializeStuff(model,doid_range,surface,seed,parameters):
    '''
    called by the master node in a parallel job to:
    1. get the list of finite element nodes
    2. randomly chop up the list into #proc's equal lists
    3. if monte = 1 (monte carlo simulation with random initial crack sizes
       coming from analytical distribution) generate the initial crack sizes
    '''
    node_list = model.GetNodeList()
    doid_list = MakeLists(node_list,model,doid_range,surface)

    # randomly chop up the doid_list for processing... 
    DoParallel(doid_list,node_list,seed)
    
    # if monte == 1, generate and write initial a's
    #  note:  if this is run from the data_base, this list will not
    #  get used because Gerd supplies the .rnd file...
    if parameters.monte == 1:
        ais = Monte(parameters.a_b,parameters.sets, \
                                   parameters.samples,seed)
        tfile=open(filename+'.rnd','w')
        s=''
        for diction in ais.rv:
            for ai in diction.keys():
                for rid in diction[ai]:
                    s+=str(rid)+' '+str(ai)+"\n"
##                print s
        s=s[0:-1]
        tfile.write(s)
        tfile.close()

#############

def DataBaseInitializeStuff(model,extension,parameters):
    '''
    from a process in a parallel job
    '''
    node_list=model.GetNodeList()

    # node_partition.%MSTI_RANK% should be there... 
    tfile=open('node_partition.'+extension,'r')
    doid_list=[]
    lines=tfile.readlines()
    for line in lines:
        doid_list+=[int(line)]
    tfile.close()

    # input.rnd should be here.  map for rid -> ai
    tfile=open("input.rnd",'r')
    ais=Statistic.Stata(parameters.sets,parameters.samples)
    for i in range(parameters.sets):
        ais.UpdateSet()
        for j in range(parameters.samples):
            l=tfile.readline()
            l=string.splitfields(l)
            ais.UpdateSamples(float(l[1]),i,int(l[0]))
    tfile.close()

    ais_map = None # remains None unless monte = 2
    if parameters.monte == 2: 
        # input.map should be here too. creates the map for doid -> rid
        tfile=open("input.map",'r')
        ais_map={}
        tfilelist=tfile.readlines()
        for line in tfilelist:
            l = string.splitfields(line)
            ais_map[int(l[0])]=map(int,l[2:])
        tfile.close()

    ncr=Numeric.arange(0,parameters.N_max*1.01,1000)
    ncr=ncr.tolist()

    return doid_list,node_list,ais,ncr,ais_map

#############

def InitializingStuff(model,user_supplied_doid,pickled_list,doid_range,\
                      surface,parameters,user_supplied_ais,parpath,filename,\
                      seed,savepickleall,default,example,base):

    # Get node_list & doid_list...
    node_list=model.GetNodeList()

    if user_supplied_doid:
        index=sys.argv.index('-doid_list')+1 # to skip '('
        doid_list=[]
        entry=sys.argv[index]
        dlist=string.splitfields(entry,',')
        for dd in dlist: doid_list+=[int(dd)]

    elif pickled_list:
        index=sys.argv.index('-pick_list')+1
        pick_list=sys.argv[index]
        tfile=open(pick_list,'r')
        # old code
##        (doid_list,node_list) = cPickle.load(tfile)
##        tfile.close()
        # node_partition way...
        doid_list=[]
        lines=tfile.readlines()
        for line in lines:
            doid_list+=[int(line)]
        tfile.close()

    else: doid_list = MakeLists(node_list,model,doid_range,surface)

    # if this is a monte carlo simulation, get the initial a's and critical
    # life range (ncr)
    ais=None; ais_map=None; ncr = None
    if parameters.monte == 1 or parameters.monte == 2:
        if parameters.monte == 1: 
            if user_supplied_ais:
                tfile=open(parpath+filename+'.rnd','r')
                ais=Statistic.Stata(parameters.sets,parameters.samples)
                for i in range(parameters.sets):
                    ais.UpdateSet()
                    for j in range(parameters.samples):
                        l=tfile.readline()
                        l=string.splitfields(l)
                        ais.UpdateSamples(float(l[1]),i,int(l[0]))
                tfile.close()
            else:
                ais = Monte(parameters.a_b,parameters.sets,parameters.samples,\
                            seed)

        else: # paramaters.monte == 2!
            # input.rnd should be here.  map for rid -> ai
            tfile=open(parpath+filename+'.rnd','r')
            ais=Statistic.Stata(parameters.sets,parameters.samples)
            for i in range(parameters.sets):
                ais.UpdateSet()
                for j in range(parameters.samples):
                    l=tfile.readline()
                    l=string.splitfields(l)
                    ais.UpdateSamples(float(l[1]),i,int(l[0]))
            tfile.close()

            ais_map = None # remains None unless monte = 2
            if parameters.monte == 2: 
                # input.map should be here too. creates the map for nid -> rid
                tfile=open(parpath+filename+".map",'r')
                ais_map={}
                tfilelist=tfile.readlines()
                for line in tfilelist:
                    l = string.splitfields(line)
                    ais_map[int(l[0])]=map(int,l[2:])
                tfile.close()

        # create the list of critical life values used to compute probabilities
        if default == 1 or example == 1 or base == 1:
            ncr=Numeric.arange(0,parameters.N_max*1.01,1000)
            ncr=ncr.tolist()
        else:
            min_ncr,max_ncr,space = GetSatisfied(GetncrRange)
            ncr=Numeric.arange(min_ncr,max_ncr,space)
            ncr=ncr.tolist()

    return doid_list,node_list,ais,ncr,ais_map

#############    

def Kva_History(ais,model,cracks,num_bcs,parameters,verbose,\
                saveintermediate,parpath,filename,extension,doid_list,\
                errfile_extension,ncr,Scale):
    '''
    Builds the K v. a curve and integrates for life prediction.
    '''

    # Get the smallest initial a
    if parameters.monte:
        initial_a=ais.rv
        init_a_sim=[]
        for set in initial_a:
            set.sort()
            init_a_sim+=set

        init_a_sim.sort()

        # begin with smallest ai in random list.
        a=init_a_sim[0] 
        b=a

        # Calc number of iterations to meet max crack size (user specified)
        # (r^n)*a=a_max
        # add 2 so that if changes damage type still gets to a_max
        iters = Numeric.log(parameters.Max_crack_size/a)/ \
                Numeric.log(parameters.r)
        iters = int(Numeric.ceil(iters))+2

    # Build and integrate the K v. a curve 
    did = 0

    for doid in doid_list: # loop over damage origins (i.e. nodes)
        if verbose: print "------ doid:", doid," --------"
        xyz,delxyz,sigxyz=model.GetNodeInfo(doid) #@
        cracks.AddDamOro(doid,xyz) #@
        # try to catch the doid causing the problem! 
##        try: 
        for j in xrange(num_bcs): # loop over boundary conditions 
            caseid=j+1
            cracks.AddDamdid(doid,did,caseid) #@
            # loop over different initial a's, and b's
            for k in xrange(len(parameters.a_b)):
                # loop over init_a_sim to get simulation for smallest ai.  
                if parameters.monte == 1: 
                    # loop over init_a_sim to get simulation for smallest
                    # ai.  
                    dam_history={'ab':[[a,b]], '-ab':[[a,b]]}
                    cracks.AddFDam(doid,did,dam_history)
                    errfile_name=parpath+filename+errfile_extension
                    cracks.BuildDkva(doid,did,parameters.material, \
                                     errfile_name, parameters.r,iters,\
                                     Scale)
                    # simulation complete.  
                    # First calculate the total life for this doid by 
                    # integrating from the smallest ai to the end.  If
                    # will_grow parameter is -1, location is in
                    # compression assume no crack will grow from there...
                    if cracks.GiveMeWG(doid) == -1:
                        N = 1.01*parameters.N_max
                        cracks.UpdateN(doid,did,a,N)
                        # now loop over the sets
                        for setid in xrange(len(initial_a)):
                            # add a blank list to set data in statistics  
                            # instances in __DamOro
                            cracks.UpdateSet(doid)  
                            for ai in initial_a[setid]:
                                cracks.UpdateSample(doid,ai,N,setid)
                            cracks.CalcStats(doid,setid,ncr)

                    else:
                        for aG in init_a_sim:
                            N_tot = cracks.IntegrateLife(aG,did,doid,\
                                                         'end', \
                                                         parameters.N_max)
                            if N_tot == -1: continue
                            else: break

                        # after determining N_tot by looping over all the
                        # initial a's, if N_tot is -1, means won't grow for
                        # even the largest initial crack, make all lives =
                        # 1.01*N_max
                        if N_tot == -1: N_tot = 1.01*parameters.N_max

                        cracks.UpdateN(doid,did,aG,N_tot)
                        # now loop over the sets
                        for setid in xrange(len(initial_a)):
                            # add a blank list to set data in statistics  
                            # instances in __DamOro
                            cracks.UpdateSet(doid)  
                            # loop over initial crack sizes (they are
                            # sorted).
                            aL=a
                            N=N_tot
                            # calc life from previous ai, aL (a left) to 
                            # this one, aR (a right) and subtract from
                            # N_tot
                            for aR in initial_a[setid]:
                                # Sometimes, we can get Nminus from cracks 
                                # with ai < aG because the interval [aL,aR]
                                # it too short to get life up in calc. of 
                                # effective R calc.  Exclude them...
                                if aR < aG:
                                    cracks.UpdateSample(doid,aR, \
                                           1.01*parameters.N_max,setid)
                                    aL=aR
                                    continue 

                                elif aR == aG:
                                    cracks.UpdateSample(doid,aR, \
                                                 N_tot,setid)
                                    aL=aR
                                    continue

                                Nminus=cracks.IntegrateLife(aG,did,doid,\
                                                          aR, \
                                                          parameters.N_max)

                                if Nminus < 0.0:
                                    cracks.UpdateSample(doid,aR, \
                                                 1.01*parameters.N_max, \
                                                 setid)

                                else:
                                    N=N_tot-Nminus
                                    cracks.UpdateSample(doid,aR,N,setid)

                                aL=aR

                            # check if total subtracted is roughly accurate
                            if verification:
                                Nminus=cracks.IntegrateLife(aG,did,doid,\
                                                  initial_a[setid][-1], \
                                                  parameters.N_max)
                                print ' '
                                print "Check if entire amount subtracted",\
                                      "is approximately equal to ",\
                                      "integrating"
                                print " interval [aG,"\
                                      "initial_a[setid][-1]]."
                                print " The integral for that limit is: ",\
                                      Nminus
                                print " N_tot = N + Nminus"
                                print " N_tot", N_tot
                                print " N + Nminus", N+Nminus
                                print ' '

                            cracks.CalcStats(doid,setid,ncr)

                    if verbose:
                        af,bf = cracks.GiveCurrentGeo(did)[1][0],\
                                cracks.GiveCurrentGeo(did)[1][1]
                        print ' '
                        print ' doid: ', doid, \
                              ' Life is:', cracks.GiveMeLife(doid), \
                              ' WillGrow:', cracks.GiveMeWG(doid), \
                              ' ai, bi: %1.4e, %1.4e' %(aG,aG),'-->',\
                              '%1.4e, %1.4e' %(af,bf)   
                        print ' '

                else:
                    # ***** this is not the most efficient way to do this 
                    # right now when len(parameters.a_b) > 1 Because this 
                    # will generate K vs. a for each a,b combo in list ****
                    # loop over init_a_sim to get sim'tion for smallest ai.
                    a = parameters.a_b[k][0]
                    b = parameters.a_b[k][1]
                    # Calc number of iterations to meet max crack size 
                    # (users specified).  (r^n)*a=a_max, add 2 so that if 
                    # changes damage type still gets to a_max.  
                    iters = Numeric.log(parameters.Max_crack_size/a)/ \
                            Numeric.log(parameters.r)
                    iters = int(Numeric.ceil(iters))+2
                    dam_history={'ab':[[a,b]], '-ab':[[a,b]]}
                    cracks.AddFDam(doid,did,dam_history)
                    errfile_name=parpath+filename+errfile_extension
                    cracks.BuildDkva(doid,did,parameters.material, \
                                     errfile_name, parameters.r,iters,\
                                     Scale)
                    # simulation complete.  
                    # calculate the total life for this doid by 
                    # integrating from the smallest ai to the end.
                    if cracks.GiveMeWG(doid) == -1: N = parameters.N_max
                    else:
                        N = cracks.IntegrateLife(a,did,doid,'end', \
                                                 parameters.N_max)
                    cracks.UpdateN(doid,did,a,N)

                    if verbose:
                        af,bf = cracks.GiveCurrentGeo(did)[1][0],\
                                cracks.GiveCurrentGeo(did)[1][1]
                        print ' '
                        print ' doid: ', doid, \
                              ' Life is:', cracks.GiveMeLife(doid), \
                              ' WillGrow:', cracks.GiveMeWG(doid), \
                              ' ai, bi: %1.4e, %1.4e' %(a,b),'-->',\
                              '%1.4e, %1.4e' %(af,bf)   
                        print ' '

                did+=1
##        # if we catch it, good, else, keep going!  
##        except:
##            errfile_name=parpath+filename+".D"+errfile_extension
##            errfile=open(errfile_name,'a')
##            errfile.write(str(doid)+"\n")
##            errfile.close()

#############

def TakeArgv():
    '''
    function to take the command arguments and pass the right stuff back.
    '''
    
    # -d default values (won't be exposed to an end user)
    default = 0
    if "-d" in sys.argv: default = 1

    # -base means the filename is passed in at command line.  need boolean to
    # avoid prompting at command line for additional parameters (like ncr)
    base = 0
    if "-base" in sys.argv: base = 1

    # -pn is flag automatically passed in for parallel jobs
    local_node = 0
    if "-pn" in sys.argv: local_node = 1

    # -eg do example1 (-d is not valid and ignored if used in conjunction
    # with -eg)
    example = 0
    if "-eg" in sys.argv: example = 1

    # -p means is a parallel job
    parallel = 0
    if "-p" in sys.argv: parallel = 1

    # -o means we want to be prompted for output options
    output = 0
    if "-o" in sys.argv: output = 1
    
    # -pa means we just want it to print all results at the end of the day.
    # -pdid print the damage elements (DamHistory)
    # -pdoid print the damage origin results
    printall = 0
    if "-pa" in sys.argv: printall = 1
    elif "-pdid" in sys.argv: printall = 2
    elif "-pdoid" in sys.argv: printall = 3

    # -v print damage results as they process
    verbose = 0
    if "-v" in sys.argv: verbose = 1

    # -sv means want to save all text files
    saveall = 0
    if "-sv" in sys.argv: saveall = 1

    # -sb means want to save various pickled files (mostly for monte).  
    savepickleall = 0
    if "-sp" in sys.argv: savepickleall = 1

    # -sc means want to save files for viewing contours (*.svw, and *.MP)
    savecontour = 0
    if "-sc" in sys.argv: savecontour = 1

    # -L means we want to limit the doid_list spatially... will prompt for
    # coordinates later.
    limit = 0
    if "-L" in sys.argv: limit = 1

    # -S means we only want surface nodes in the doid_list
    surface = 0
    if "-S" in sys.argv: surface = 1

    # -si means ...
    saveintermediate = 0
    if "-si" in sys.argv: saveintermediate = 1

    # -usa means use a user supplied ais class for the initial ai's
    user_supplied_ais = 0
    if "-usa" in sys.argv: user_supplied_ais = 1

    # -doid_list means use a user supplied list of doid's to inspect
    user_supplied_doid = 0
    if "-doid_list" in sys.argv: user_supplied_doid = 1

    # -pick_list means use a user supplied list via pickled list
    pickled_list = 0
    if "-pick_list" in sys.argv:    pickled_list = 1

    # specify the type of integration for a costant amplitude simulation
    # default --> adaptive, 5 pnt. runga-kutta scheme
    # -Fwd --> forward euler
    # -RK4 --> 4 point RK scheme
    # -Simp --> simpson's rule on generated SIF history
    Int_type = 'RK5'
    if "-Fwd" in sys.argv: Int_type = 'FWD'
    elif "-RK4" in sys.argv: Int_type = 'RK4'
    elif "-Simp" in sys.argv: Int_type = 'Simp'

    # -VarAmp means this is a variable amplitude loading case.  Should be
    # followed be the Variable amplitude input file name.  
    var_file = 0
    if "-VarAmp" in sys.argv: var_file = 1 

    # -seed provide a seed to generate the monte samples.  Should be followed
    # by a number
    seed = None
    if "-seed" in sys.argv:
        index = sys.argv.index('-seed')+1
        seed = float(sys.argv[index])

    # -scale provide a scale 0 < scale to scale K's as appropriategenerate.
    # Should be followed by a number
    scale = 1.0
    if "-scale" in sys.argv:
        index = sys.argv.index('-scale')+1
        scale = float(sys.argv[index])

    # -nore to prevent the use of the Willenborg retardation model. i.e. growth
    #       rate is given by NASGRO eqn w/o retardation effects.
    nore = 0.0
    if "-nore" in sys.argv:
        nore = 1.0

    # -verify get some diagnostic information printed up on the screen
    verify = False
    if "-verify" in sys.argv:
        verify = True

##    # to control the integration scheme in a constant amplitude loading
##    # scenario.  if -int follow by RK4 or FWD (case sensitive)
##    integration = 'RK5'
##    if "-int" in sys.argv:
##        index = sys.argv.index('-int')+1
##        integration = (sys.argv[index])

    # to make create a file that our script to test GeomUtils can read
    # follow key (-DebugGeomUtils) by string indicate name and location of file
    DebugGeomUtils=False
    if "-DebugGeomUtils" in sys.argv:
        index = sys.argv.index('-DebugGeomUtils')+1
        DebugGeomUtils = sys.argv[index]

    # to write an sview file.  follow with string file name/location... 
    SVIEW=False
    if "-SVIEW" in sys.argv:
        index = sys.argv.index('-SVIEW')+1
        SVIEW = sys.argv[index]

    data_base=0
    if "-DB" in sys.argv: data_base=1

    crack_front=0
    if "-cf" in sys.argv: crack_front=1
    

    return default,local_node,example,parallel,output,printall,verbose, \
           saveall,savepickleall,savecontour,limit,surface,saveintermediate,\
           base,user_supplied_ais,user_supplied_doid,Int_type,\
           pickled_list,var_file,seed,scale,nore,verify,DebugGeomUtils,SVIEW, \
           data_base,crack_front

#############

def Monte(a_b,sets,samples,seed=None):
    '''
    generate samples of initial crack size and stores in dictionary, initial_a
    also, record the smallest intial crack size for the life prediction
    simulation in a dictionary, smallest_a
    ***NOTE:  I am not sure that i need to generate sets of initial a values by
              monte carlo for each node, or whether i can simply do it once and
              use the same list of initial a's for each node.  Currently, i
              generate new sets of samples for each node.  This may prove to be
              to computationally demanding for larger models with more sets and
              samples. (9/04)
    
    (10/12)   Further, on that note:  It didn't take too long to become too
              large and get and overflow error.  Consequently, i am now just
              generating samples for n sets of m samples once and using them at
              each node. the code for the old way is located at the bottom of
              this script. The loop that does the actual DDSim simulation will
              be modified too but the changes should be simple enough to undo..
    '''
    # so i can easily do statistics on them and plot CDF
    ais=Statistic.Stata(sets,samples) 

    if seed: random.seed(seed)

    for set in xrange(sets):
        ais.UpdateSet()

        for samp in xrange(samples):
            # random to generate uni sample
            uni=random.uniform(0,1) 

##            # a two param weibull
##            a=a_b[0][0]*(-1.0*Numeric.log(1.0-uni))**(1.0/a_b[0][1])

            # 3 param weibull.
            gamma = 6.1 # from Gary Harlow for AA7075-T651
            samp=(a_b[0][0]*\
                  ((-1.0*Numeric.log(1.0-uni))**(1.0/a_b[0][1])))+gamma
            # samp is the TS area... a is the radius:  
            # 3.937e-5 converts from microns to inches
            a=3.937e-5*math.sqrt(samp/math.pi)

            ais.UpdateSamples(a,set)

##    temp=open('temp.ais','w')
##    for rid in ais.rid_ai.keys():
##        temp.write(str(rid)+" "+str(ais.rid_ai[rid][1])+"\n")
##    temp.close()

    return ais

#############

def MonteSimulation(ais,ais_map,cracks,parpath,filename,extension,\
                    errfile_extension,parameters,doid,xyz,verbose,\
                    ncr,scale,integration,verify,R):

    # Get the largest initial a
    if ais_map:
        local_ais=MakeLocalAis(ais,ais_map,doid)

    else: local_ais = ais

    initial_a = local_ais.rv
    init_a_sim=[]
    for set in initial_a:
        init_a_sim+=set.keys()
    init_a_sim.sort()

    # if len(initial_a_sim) == 0, then particle cracking filter said no
    # particles will crack at this doid, return N_max and exit
    if len(init_a_sim) == 0:
        # fake into having a damage element
        cracks.AddDamOro(doid,ais)
        cracks.AddFDam(doid,xyz,1.0,1.0,parameters.material,verbose,verify)
        for setid in range(parameters.sets):
            cracks.UpdateSet(doid)
            for a in ais.rv[setid]:
                cracks.UpdateSample(doid,a,parameters.N_max*1.01,setid)
            cracks.CalcStats(doid,setid,ncr)

        if verbose:
            print ' '
            print ' No particles cracked! ', \
                  ' doid: ', doid, \
                  ' Life is:', parameters.N_max*1.01, \
                  ' WillGrow:', 0, \
                  ' ai: None'
            print ' '

        return None

    cracks.AddDamOro(doid,local_ais)

    # want to skip ai's that already are too small to grow.  
    # the way to do this is, when a crack shuts down, check the
    # minor radius, b, and skip ai's <= b.  So initialize b to
    # min(init_a_sim)

    # run the simulation for the largest ai first, if it grows continue and run
    # the simulation for all initial ai's, otherwise quit.
    ai=init_a_sim[-1]
    bi = ai # so is obvious... 
    # now do the simulation
    cracks.AddFDam(doid,xyz,ai,bi,parameters.material,verbose,verify)
    errfile_name=parpath+filename+errfile_extension
    N=cracks.SimDamGrowth(doid,parameters,errfile_name,scale,R) 

    if verbose:
        af,bf = cracks.DamOro[doid].DamEl[-1].GiveCurrent()
        print ' '
        print ' For Largest ai @', \
              ' doid: ', doid, \
              ' Life is:', cracks.DamOro[doid].Life, \
              ' WillGrow:', cracks.DamOro[doid].WillGrow, \
              ' ai: %1.4e' %(ai),'-->',\
              "%1.4e, %1.4e" %(af,bf)
        print ' '

    # if the largest crack does NOT want to grow (.willgrow = 0), is growing
    # stably until Nmax is reached (.willgorw = 1) or is in compression
    # (.willgrow = -1), then do NOT bother checking all the other initial ai's.
    # This final case, .willgrow = -1, is debatable because it is conceivable
    # that a smaller crack may not be in compression.  I'm going with this for
    # now, it seems like a safe bet. 
    if cracks.DamOro[doid].WillGrow == 2 or cracks.DamOro[doid].WillGrow == 4:
        b=min(init_a_sim) 
        for ai in init_a_sim: # may use break to exit loop... 
            # compare to b
            if ai >= b:
                bi = ai
                # now do the simulation
                del cracks.DamOro[doid]
                cracks.AddDamOro(doid)
                cracks.AddFDam(doid,xyz,ai,bi,parameters.material,verbose, \
                               verify)
                errfile_name=parpath+filename+errfile_extension
                N=cracks.SimDamGrowth(doid,parameters,errfile_name,scale,R)

                if verbose:
                    af,bf = cracks.DamOro[doid].DamEl[-1].GiveCurrent()
                    print '\n', \
                          ' doid: ', doid, \
                          ' Life is:', N, \
                          ' WillGrow:', cracks.DamOro[doid].WillGrow, \
                          ' ai: %1.4e' %(ai),'-->',\
                          "%1.4e, %1.4e" %(af,bf)
                    print ' '

                # if willgrow = 2, unstable growth occurred: stop 
                # iterating over intial ai's
                # if willgrow = -1, compression
                # willgrow = 4, net fracture
                if cracks.DamOro[doid].WillGrow == 2 or \
                   cracks.DamOro[doid].WillGrow == -1 or \
                   cracks.DamOro[doid].WillGrow == 4: break

                # if willgrow = 0, crack shuts down: skip up ai 
                elif cracks.DamOro[doid].WillGrow == 0:
                    b = min(cracks.GiveCurrentGeo(doid))

                # if willgrow = 1, cracks still wants to grow 
                # stable: keep looping over ai's
                else: continue
            else: continue

        # simulation complete.  if monte==1, it is only complete 
        # for the smallest initial crack size that will grow.  use 
        # those simulation results to calculate the life for the 
        # other initial crack size here... 
        # loop over the sets...
        for setid in range(parameters.sets):
            # add a blank list to set data in statistics instances 
            # in __DamOro
            cracks.UpdateSet(doid)  
            # loop over initial crack sizes
            for a in initial_a[setid]:
                if ai==a: cracks.UpdateSample(doid,a,N,setid)
                else: 
                    Ni=cracks.InterpolateLife(a,doid,parameters.N_max)
                    cracks.UpdateSample(doid,a,Ni,setid)
            cracks.CalcStats(doid,setid,ncr)

    else: # largest crack will not grow
        for setid in range(parameters.sets):
            cracks.UpdateSet(doid)
            for a in initial_a[setid]:
                cracks.UpdateSample(doid,a,parameters.N_max*1.01,setid)
            cracks.CalcStats(doid,setid,ncr)

#############

def GetBox():
    '''
    define the bounding box for the spaitally limited simulation.  only calls
    this function if -L is in command line.
    '''
    print ' '
    print ' Must define a bounding box to limit the selection of nodes'
    print ' to be included.  This is done with two points.  One defines'
    print ' the minimum x,y,z coords, and the other defines the maximum.'
    print ' '
    print 'Enter minimum x coord'
    xmin = float(sys.stdin.readline().split()[0])
    print 'Enter minimum y coord'
    ymin = float(sys.stdin.readline().split()[0])
    print 'Enter minimum z coord'
    zmin = float(sys.stdin.readline().split()[0])
    print 'Enter maximum x coord'
    xmax = float(sys.stdin.readline().split()[0])
    print 'Enter maximum y coord'
    ymax = float(sys.stdin.readline().split()[0])
    print 'Enter maximum z coord'
    zmax = float(sys.stdin.readline().split()[0])
    return ((xmin,ymin,zmin),(xmax,ymax,zmax))

#############

def GetSetNumber():
    '''
    get the set number the user is interested in.
    '''
    print ' Enter set to display mean and variance for (-1 to exit):'
    setid = int(sys.stdin.readline().split()[0])
    return setid

#############

def GetName():
    '''
    get the input file name and location.
    ''' 
    print " Enter file name"
    filename = sys.stdin.readline().split()[0]
    print " Enter *.con file path"
    conpath = sys.stdin.readline().split()[0]
    print " Enter *.par file path"
    parpath = sys.stdin.readline().split()[0]
    return filename,conpath,parpath

#############

def GetncrRange():
    '''
    get critical life range for the statistics
    ''' 
    print " minimum acceptable crictical life "
    min_ncr = float(sys.stdin.readline().split()[0])
    print " maximum acceptable critical life"
    max_ncr = float(sys.stdin.readline().split()[0])
    print " spacing"
    space = float(sys.stdin.readline().split()[0])
    
    return min_ncr,max_ncr,space

#############

def GetSatisfied(function):
    '''
    loop over function until user has specified arguements to their
    satisfaction
    '''
    arg=function()
    satisfied = 0
    while satisfied == 0:
        print ' you entered', arg
        print " Is this correct? (y/n)"
        response = sys.stdin.readline()
        if 'y' in response or 'Y' in response:
            satisfied = 1
            continue
        else:
            arg=function()
    return arg

#############

def MakeLists(node_list,model,doid_range,surface):
    doid_list=[]
    node_list.sort()
    for node in node_list:
        if doid_range: # pick damage origins from this domain only
            xyz,delxyz,sigxyz=model.GetNodeInfo(node)
            x_val=xyz.x()
            y_val=xyz.y()
            z_val=xyz.z()
            if x_val >= doid_range[0][0] and x_val <=doid_range[1][0] and \
               y_val >= doid_range[0][1] and y_val <=doid_range[1][1] and \
               z_val >= doid_range[0][2] and z_val <=doid_range[1][2]:
                if surface:
                    if model.IsSurfaceNode(node): doid_list+=[node]
                else: doid_list+=[node]
        elif surface:
            if model.IsSurfaceNode(node): doid_list+=[node]
        else:
            doid_list+=[node]

    doid_list.sort()

    return doid_list

#############

def MakeLocalAis(ais,ais_map,doid):
    '''
    using ais (rid -> ai) and ais_map (nid -> rid) make local_ais...
    '''

    # there's only one set!!  COMPASS keeps track how rid's map to sets.
    if ais_map.has_key(doid):
        local_ais=Statistic.Stata(1,len(ais_map[doid]))
        local_ais.UpdateSet() 
        for RID in ais_map[doid]:
            local_ais.UpdateSamples(ais.rid_ai[RID][1],0,RID)
    else: local_ais = Statistic.Stata(0,0) # a dummy

    return local_ais

#############

def DoParallel(doid_list,node_list,seed=None):
    '''
    subroutine that chops up the doid_list.
    '''

    # randomly scramble the doid_list; a simple measure for speedup. 
    if seed: random.seed(seed)
    new_doidlist=[]
    while len(doid_list) > 0:
        ii=random.randint(0,len(doid_list)-1)
        new_doidlist+=[doid_list.pop(ii)]

    # loop over length of machines file to build doid_pickle.*'s
    numprocs = int(os.environ['procs'])
    old = 0
    for i in xrange(numprocs-1):
        count = (i+1)*(len(new_doidlist)/numprocs)
        templist=new_doidlist[old:count]
        templist.sort()
        tfile=open('node_partition.'+str(i),'w')
        for no in templist:
            tfile.write(str(no)+"\n")
        tfile.close()

        del(templist)
        old = count

    templist=new_doidlist[old:]
    templist.sort()

    tfile=open('node_partition.'+str(numprocs-1),'w')
    for no in templist:
        tfile.write(str(no)+"\n")
    tfile.close()

    del(templist) 

#############

def FwdDeterministic(cracks,doid,xyz,parameters,verbose,verify,parpath,\
                     filename,errfile_extension,scale):

    cracks.AddDamOro(doid)
    a=parameters.a_b[0][0] 
    b=parameters.a_b[0][1]
    # now do the simulation
    cracks.AddFDam(doid,xyz,a,b,parameters.material,verbose,verify)
    errfile_name=parpath+filename+errfile_extension
    N=cracks.SimDamGrowth(doid,parameters,errfile_name,scale,\
                        parameters.material[16])

    if verbose:
        af,bf = cracks.DamOro[doid].DamEl[-1].GiveCurrent()
        print ' '
        print ' doid: ', doid, \
              ' Life is:', N, \
              ' WillGrow:', cracks.DamOro[doid].WillGrow, \
              ' ai, bi: %1.4e, %1.4e' %(a,b),'-->',\
              '%1.4e, %1.4e' %(af,bf)   
        print ' '


#############

def Fwd_Integration(ais,model,cracks,num_bcs,parameters,verbose,\
                    saveintermediate,parpath,filename,extension,doid_list,\
                    errfile_extension,ncr,scale,integration,verify,ais_map):
    '''
    Fwd_Integration is the main processing loop in the code.  here is where i
    loop through doid_list and predict life for each node in doid_list.
    '''

    for doid in doid_list: # loop over damage origins (i.e. nodes)

        if verbose: print "------ doid:", doid," --------"
        xyz,delxyz,sigxyz=model.GetNodeInfo(doid) 

##        try: 
        if parameters.monte == 1 or parameters.monte == 2:
            MonteSimulation(ais,ais_map,cracks,parpath,filename,extension, \
                            errfile_extension,parameters,doid,xyz,\
                            verbose,ncr,scale,integration,verify,\
                            parameters.material[16])

        # if not monte carlo simulation then a_b is given and contains initial
        # crack size information; strip it out here and do life prediciton
        else: FwdDeterministic(cracks,doid,xyz,parameters,verbose,verify,\
                              parpath,filename,errfile_extension,scale)

        # ---> LIfe prediction for doid complete.  Now post process! 
        if data_base or local_node:
            nfile=open('output.N.'+extension,'a')
            cracks.NFile(doid,nfile)
            nfile.close()

            OrientFile = open(parpath+filename+'.ori.'+extension,'a')
            cracks.WriteRotations(OrientFile,doid)
            OrientFile.close()

            AiFile = open(parpath+filename+'.ai.'+extension,'a')
            cracks.WriteInitialAs(AiFile,parameters.N_max,doid)
            AiFile.close()

            AfFile = open(parpath+filename+'.af.'+extension,'a')
            cracks.WriteFinalAs(AfFile,doid)
            AfFile.close()

            # clear the stored crap that we don't need because we will get 
            # memory errors otherwise.  
            del cracks.DamOro[doid]

        elif cracks.saveall:
            nfile=open(parpath+filename+'.N','a')
            cracks.NFile(doid,nfile)
            nfile.close()

            OrientFile = open(parpath+filename+'.ori','a')
            cracks.WriteRotations(OrientFile,doid)
            OrientFile.close()

            AiFile = open(parpath+filename+'.ai','a')
            cracks.WriteInitialAs(AiFile,parameters.N_max,doid)
            AiFile.close()

            AfFile = open(parpath+filename+'.af','a')
            cracks.WriteFinalAs(AfFile,doid)
            AfFile.close()

##        except:
##            print 'exception caught in DDSim.Fwd_Integration at doid', doid
##            blah=open(parpath+filename+'.err.'+extension,'a')
##            blah.write(str(doid)+"\n")
##            blah.close()

#############

def GetFileName(example,local_node,default,parallel,base,argv,data_base):
    '''
    function to return strings:
    filename
    conpath - location where input files, *.con, *.sig, *.nod, *.edg exist
    parpath - location where input file *.par and batch file exist and where
              ouput will be written
    '''

    if example:
        filename='''example1'''
##        conpath="c:\\john\\research\\ureti\\code\\" + \
##                "python\\ddsimv1.4\\examples\\"
##        conpath="z:\\research\\ureti\\code\\python" + \
##                "\\ddsimv1.3\\examples\\"
        conpath="h:\\users\\jme32\\research\\ureti\\" + \
                "code\\python\\ddsimv1.4\\examples\\"
        parpath=conpath
    # for parallel job when executed on individual nodes.
    elif data_base:
        filename='''input'''
        conpath=''''''
        parpath=conpath
    elif local_node and default == 0: 
        filename='''input'''
        conpath=''''''
        parpath=conpath
    # for parallel job (input file names come from the batch file)
    elif parallel:
        filename=os.environ['job_name']
        conpath=os.environ['PROJECT_IN']+'\\'
        parpath=os.environ['PROJECT_OUT']+'\\'
    # a default to facilitate debugging
    elif default and parallel == 0:
##        filename='''coupon01'''
##        path='''c:\\john\\research\\SIPS\\franc3d\\coupon01\\'''
##        filename='''cube_20x20_pull'''
##        path='''c:\\john\\Gerd_FEM\\cube_20x20_pull\\'''
        filename='''Turbine_disc_15mid'''
        conpath='''Z:\\research\\URETI\\franc3d\\turbine_disc\\'''
        conpath+='''Norm_BC_15Trac_middle\\Determ_015\\'''
        parpath='''Z:\\research\\URETI\\franc3d\\turbine_disc\\'''
        parpath+='''Norm_BC_15Trac_middle\\Determ_015\\'''
    # ordinary serial run mode with file name at command line
    elif base:
        # base = 1 in TakeArgv to avoid prompts for ncr etc. 
        index = argv.index('-base')+1
        filename = argv[index]
        if "-conpath" in argv:
            index = argv.index('-conpath')+1
            conpath = argv[index]
        else:
            conpath=''''''
        if "-parpath" in argv:
            index = argv.index('-parpath')+1
            parpath = argv[index]
        else:
            parpath=''''''
    # ordinary serial run mode prompt user for file name
    else: (filename,conpath,parpath)=GetSatisfied(GetName)

    return filename,conpath,parpath

###############

def Var_Amplitude(ais,model,cracks,parameters,verbose,\
                 saveintermediate,parpath,conpath,filename,extension,\
                 doid_list,errfile_extension,ncr,Scale,nore,data_base,\
                local_node,ais_map):
    ''' 
    cycle-by-cycle variable amplitude loading.  for monte carlo simulation,
    N is first computed for the largest initial crack size.  the life for the
    other initial crack sizes is computed as the number of cycles to growth to
    the next largest crack size plus the other N's

    for parallel jobs, or database initiate jobs writes job_nam.N file that
    contains:

    doid rid N

    rules of the *.N file: 
    1.  at least one row for every doid in doid_list
    2.  If no computed N's are less than N_max, write rid = -1 and N_max
    '''

    Spec=VarAmplitude.Spectrum(conpath+filename+'.val')
    N_max = 1.01*parameters.N_max

    for doid in doid_list: # loop over damage origins (i.e. nodes)
        if verbose: print "------ doid:", doid," --------"
        xyz,delxyz,sigxyz=model.GetNodeInfo(doid)
        # try to catch the doid causing the problem! 
##        try: 
        # loop over init_a_sim to get simulation for smallest ai.  
        if parameters.monte == 1 or parameters.monte == 2:

            # Get the largest initial a
            if ais_map:
                local_ais=MakeLocalAis(ais,ais_map,doid)

            else: local_ais = ais

            initial_a = local_ais.rv

            # if len(initial_a) == 0, particle filter found no particles
            # will break at this doid set cracks.DamOro[doid].WillGrow == -1
            # (which is caught below) and don't do anything else for this doid
            if len(initial_a) == 0:
                cracks.AddDamOro(doid,ais)
                cracks.AddFDam(doid,xyz,1.0,1.0,parameters.material,verbose, \
                               verify)
                initial_a = ais.rv
                cracks.DamOro[doid].WillGrow = -1

            # otherwise, do the variable amplitude simulation for the
            # largest crack
            else: 

                init_a_sim=[]
                for set in initial_a:
                    init_a_sim+=set.keys()

                init_a_sim.sort()
                init_a_sim.reverse()
                # begin with largest ai in random list.
                ai=init_a_sim[0]

                # creat new Damage Origin... 
                cracks.AddDamOro(doid,local_ais)

                # loop over init_a_sim to get simulation for smallest
                # ai.  
                aR=parameters.Max_crack_size
                cracks.AddFDam(doid,xyz,ai,ai,parameters.material,verbose, \
                               verify)
                if verbose:
                    print ' Performing cycle-by-cycle integration for variable'
                    print '    amplitude loading for largest ai '
                N_tot,fin = cracks.VarAmp(doid,aR,parameters.N_max,\
                                          Spec,nore,Scale,parameters.r)
                akeep=ai

                if verbose:
                    af,bf,ka,kb=fin[0],fin[1],fin[2],fin[3]
                    print ' '
                    print ' For Largest ai:'
                    print ' doid: ', doid, \
                          ' Life is:', N_tot, \
                          ' WillGrow:', cracks.DamOro[doid].WillGrow
                    print ' ai, bi: %1.4e, %1.4e' %(ai,ai),'-->',\
                          ' af, bf: %1.4e, %1.4e' %(af,bf)
                    print ' K(a), K(b): %1.4e, %1.4e' %(ka,kb)
                    print ' '

            # now post process the results of variable amplitude loading from
            # the largest ai.  

            # if will_grow parameter is -1 or 0, largest crack will not grow
            # populate N's with N_max
            if cracks.DamOro[doid].WillGrow == -1 or \
               cracks.DamOro[doid].WillGrow == 0:
                cracks.DamOro[doid].Life=N_max
                # now loop over the sets
                for setid in range(parameters.sets):
                    # add a blank list to set data in statistics  
                    # instances in __DamOro
                    cracks.UpdateSet(doid)
                    for a in initial_a[setid]:
                        cracks.UpdateSample(doid,a,N_max,setid)
                    cracks.CalcStats(doid,setid,ncr)

            else:
                cracks.DamOro[doid].Life=N_tot
                for setid in range(parameters.sets):
                    cracks.UpdateSet(doid)
                    aR=akeep
                    N_cum=0.0
                    # assume that when one crack returns WillGrow==0 or -1
                    # the smaller ones will too.  Use ffwd to check it...
                    ffwd=0 
                    reversed_alist=initial_a[setid].keys()
                    reversed_alist.sort()
                    reversed_alist.reverse()
                    for a in reversed_alist:
                        if ffwd==0:
                            cracks.DamOro[doid].Flush()
                            cracks.AddFDam(doid,xyz,a,a,\
                                           parameters.material,\
                                           verbose,verify)
                            N,fin = cracks.VarAmp(doid,aR, \
                                             parameters.N_max-N_tot-N_cum,\
                                             Spec,nore,Scale,parameters.r)

                        aR=a

                        if cracks.DamOro[doid].WillGrow == 0 or \
                           cracks.DamOro[doid].WillGrow == -1:
                            Nai=N_max
                            ffwd=1

                        elif cracks.DamOro[doid].WillGrow == 1:
                            Nai=N+N_tot+N_cum
                            N_cum+=N

                        else: Nai=N

                        cracks.UpdateSample(doid,a,Nai,setid)

                        if verbose:
                            af,bf,ka,kb=fin[0],fin[1],fin[2],fin[3]
                            print ' '
                            print ' doid: ', doid, \
                                  ' Life (for this interval):', N, \
                                  ' WillGrow:', cracks.DamOro[doid].WillGrow
                            print ' ai, bi: %1.4e, %1.4e' %(a,a),'-->',\
                                  ' af, bf: %1.4e, %1.4e' %(af,bf)
                            print ' K(a), K(b): %1.4e, %1.4e' %(ka,kb)
                            print ' '

                    cracks.CalcStats(doid,setid,ncr)

        else: # as in, deterministic fatigue simulation...
            cracks.AddDamOro(doid)
            ai=parameters.a_b[0][0]
            aR=parameters.Max_crack_size
            cracks.AddFDam(doid,xyz,ai,ai,parameters.material,verbose, \
                           verify)
            N_tot,fin = cracks.VarAmp(doid,aR,parameters.N_max,\
                                      Spec,nore,Scale,parameters.r)

            if verbose:
                af,bf,ka,kb=fin[0],fin[1],fin[2],fin[3]
                print ' '
                print ' doid: ', doid, \
                      ' Life is:', N_tot, \
                      ' WillGrow:', cracks.DamOro[doid].WillGrow
                print ' ai, bi: %1.4e, %1.4e' %(ai,ai),'-->',\
                      ' af, bf: %1.4e, %1.4e' %(af,bf)
                print ' K(a), K(b): %1.4e, %1.4e' %(ka,kb)
                print ' '

        # life prediction for doid is complete ---> now postprocess.

        if data_base or local_node:
            nfile=open('output.N.'+extension,'a')
            cracks.NFile(doid,nfile)
            nfile.close()

            OrientFile = open(parpath+filename+'.ori.'+extension,'a')
            cracks.WriteRotations(OrientFile,doid)
            OrientFile.close()

            AFile = open(parpath+filename+'.ai.'+extension,'a')
            cracks.WriteInitialAs(AFile,parameters.N_max,doid)
            AFile.close()

            AfFile = open(parpath+filename+'.af.'+extension,'a')
            cracks.WriteFinalAs(AfFile,doid)
            AfFile.close()

            # clear the stored crap that we don't need because we will get 
            # memory erros otherwise.  
            del cracks.DamOro[doid]

        elif cracks.saveall:
            nfile=open(parpath+filename+'.N','a')
            cracks.NFile(doid,nfile)
            nfile.close()

            OrientFile = open(parpath+filename+'.ori','a')
            cracks.WriteRotations(OrientFile,doid)
            OrientFile.close()

            AiFile = open(parpath+filename+'.ai','a')
            cracks.WriteInitialAs(AiFile,parameters.N_max,doid)
            AiFile.close()

            AfFile = open(parpath+filename+'.af','a')
            cracks.WriteFinalAs(AfFile,doid)
            AfFile.close()

##        # if we catch it, good, else, keep going!  
##        except:
##            print "BARF--->:",doid 
##            errfile_name=parpath+filename+".D"+errfile_extension
##            errfile=open(errfile_name,'a')
##            errfile.write(str(doid)+"\n")
##            errfile.close()

############## top ##################

if __name__ == "__main__": 
    '''
    main coding stars here... 
    '''

    time.clock()

    # a helper for the keys...
    if "-help" in sys.argv or "-h" in sys.argv:
        HelpPrints()
        sys.exit()

    # read and parse the command line arguements
    default,local_node,example,parallel,output,printall,verbose, \
    saveall,savepickleall,savecontour,limit,surface,saveintermediate,\
    base,user_supplied_ais,user_supplied_doid,Int_type,pickled_list, \
    var_file,seed,scale,nore,verify,DebugGeomUtils,SVIEW,data_base,\
    crack_front = TakeArgv()

    if verbose:
        print '            ---------------------------------------------------'
        print '                    DDSim Level I - version 1.5 (M)            '
        print '             This run began at:', time.asctime(time.localtime())
        print '            ---------------------------------------------------'
        print ''

    num_bcs = 1 # hard code number of boundary condition sets to 1

    # get the filename/location if not supplied as command line arguments
    filename,conpath,parpath = GetFileName(example,local_node,default,\
                                           parallel,base,sys.argv,data_base)

    # do some processing input arguments and prompt for more information as
    # necessary
    if surface == 1:
        print '******************* surf_only = 1  *******************'

    doid_range = None
    if limit == 1 and parallel == 0 and local_node == 0:
        print '******************* limit_doid = 1 *******************'
        if default == 1:
            doid_range=((6.0,0.37,-0.05),(6.1,0.46,5.0))
        else: doid_range=GetSatisfied(GetBox)

    # read the path+filename.par file to add the parameters
    parameters = Parameters.Parameters(filename,parpath)

    # read and store the finite element model
    model = MeshTools.MeshTools(conpath+filename,'RDB')
    model.SetPointInsideTolerance(1.0e-7)

    # make file extension  
    extension = '0'; errfile_extension='0'; 
    if local_node or data_base:
        extension=os.environ['MSTI_RANK']
        errfile_extension='.err.'+extension
    else:
        errfile_extension='.err'

    # prepare the lists and initial a's...
    # in parallel, run on host node before actual crack growth simulations
    # begin.  
    if parallel: 
        ParallelInitializeStuff(model,doid_range,surface,seed,parameters)
        sys.exit()
    # if doid lists and initial crack sizes come from COMPASS... 
    elif data_base:
        doid_list,node_list,ais,ncr,ais_map = DataBaseInitializeStuff(model,\
                                                          extension,parameters)
    # local_node - run this for a local node in a parallel job... 
    elif local_node:
        doid_list,node_list,ais,ncr,ais_map = LocalNodeInitializeStuff(\
                                              model,extension,parameters)
    else: 
        doid_list,node_list,ais,ncr,ais_map = \
            InitializingStuff(model,user_supplied_doid,pickled_list, \
                              doid_range,surface,parameters, \
                              user_supplied_ais,parpath,filename,seed, \
                              savepickleall,default,example,base)

    if verbose:
        print ' arguments:', sys.argv
        print parameters
        print "Random initial a's list:"
        print ais

    # create the DamModel object
    cracks = DamModel.DamModel(model,node_list,verbose,verify,extension, \
                               saveall,parameters,DebugGeomUtils,SVIEW,\
                               Int_type)

    # now do the life prediction simulation with the appropriate ais
    if Int_type=="Simp":
        print " Using a Simpson's integration rule..."
        print ' '
        print "  --> WARNING!  This option is dead in this version. "
        print "      It is much faster and more accurate to just run the RK-5"
        print "      forward integration routine for constant amplitude "
        print "      loading.  Exiting now..."
        sys.exit()
        Kva_History(ais,model,cracks,num_bcs,parameters,verbose,\
                    saveintermediate,parpath,filename,extension,doid_list,\
                    errfile_extension,ncr,scale)
    elif var_file:
        if nore:
            print ' Variable Amplitude loading... No retardation model'
            print ' '
        else:
            print ' Variable Amplitude loading... Willenborg retardation model'
            print ' '
        Var_Amplitude(ais,model,cracks,parameters,verbose,\
                      saveintermediate,parpath,conpath,filename,extension, \
                      doid_list,errfile_extension,ncr,scale,nore,data_base, \
                      local_node,ais_map)
    else: 
        print " Using an Adaptive RK-5 Forward integration scheme..."
        print ' '
        Fwd_Integration(ais,model,cracks,num_bcs,parameters,verbose,\
                        saveintermediate,parpath,filename,extension,doid_list,\
                        errfile_extension,ncr,scale,Int_type,verify,ais_map)

    # some post-processing... 
    if printall==1: # -pa
        print ' '
        print ' The following is from DamModel.py... '
        cracks.PrintDamInfo('all','all',parameters.monte,'all')
    elif printall==2: # -pdid
        print ' '
        print ' The following is from DamModel.py... '
        cracks.PrintDamInfo('none','all',parameters.monte,'all')
    elif printall==3: # -pdoid
        print ' '
        print ' The following is from DamModel.py... '
        cracks.PrintDamInfo('all','none',parameters.monte,'all')

    if savepickleall: # -sp, binary pickled files. pass file object
        if data_base or local_node: pass
        elif parameters.monte:
            cracks.NToPickle(parpath+filename+'.stN.'+extension,'all')

            OrientFile = open(parpath+filename+'.ori.'+extension,'a')
            cracks.WriteRotations(OrientFile,'all')
            OrientFile.close()

            AFile = open(parpath+filename+'.ai.'+extension,'a')
            cracks.WriteInitialAs(AFile,parameters.N_max,'all')
            AFile.close()

            AfFile = open(parpath+filename+'.af.'+extension,'a')
            cracks.WriteFinalAs(AfFile,'all')
            AfFile.close()

        else:
            OrientFile = open(parpath+filename+'.ori','a')
            cracks.WriteRotations(OrientFile,'all')
            OrientFile.close()

            AFile = open(parpath+filename+'.ai.','a')
            cracks.WriteInitialAs(AFile,parameters.N_max,'all')
            AFile.close()

            AfFile = open(parpath+filename+'.af.','a')
            cracks.WriteFinalAs(AfFile,'all')
            AfFile.close()

    if savecontour: # -sc
        if parameters.monte:
            for i in xrange(parameters.sets):
                MAPFile = open(parpath+filename+'.MP'+str(i),'w')
                cracks.ToMAPFile(MAPFile,i)
        else:
            MAPFile = open(parpath+filename+'.MP','w')
            cracks.ToMAPFile(MAPFile)

    # print crack front points to a text file named filename.front
    if crack_front: # -cf
        frontfile=open(parpath+filename+'.front','w')
        cracks.FrontPoints(frontfile)
        frontfile.close()

    # store the total processing time
    if local_node or data_base:
        #variable extension named on line 545 above... msti_rank
        Timefile=open(parpath+filename+'.time.'+extension,'w')
        Timefile.write(str(time.clock())+"\n")
        Timefile.close()
    else:
        Timefile=open(parpath+filename+'.time','w')
        Timefile.write(str(time.clock())+"\n")
        Timefile.close()

    if verbose: 
        print ''
        print '            ---------------------------------------------------'
        print '            This run ended at:', time.asctime(time.localtime())
        print '            For a total elapsed time of (s):', time.clock()
        print '            ---------------------------------------------------'

############## bottom ##################
