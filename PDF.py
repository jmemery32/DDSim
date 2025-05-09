import MydadN, Statistic
import dadN
import pylab, math, sys, os, string, Numeric, copy
import JohnsVectorTools as JVT


########################################################

def CondProb(Nlist,n,ncrlist):
    '''
    compute P(N<=ncr|ai)

    Nlist = N[nid].rv[set].values() (a list of N's)
    n = total # of initial a's
    ncrlist = list of critical life in N <= ncr 
    '''
    m=0.0
    p=[]
    Nindex=0
    Nlist.sort()
    for ncr in ncrlist:
        for N in Nlist[int(m):]:
            if N<=ncr:
                m+=1.0
            else: pass
        p+=[float(m)/float(n)]

    return p

########################################################

def PlotReliability(N,samples,Nmax,color,lab,ncrlist,no_broken_particles):

    nodes=N.keys()
    numnodes=float(len(nodes))
    ncrlist=ncrlist.tolist()
    ncrlist[-1]=ncrlist[-1]+100.0

    total_broken_parts=0
    for statns in N.values(): total_broken_parts+=statns.samples

    P=pylab.zeros(len(ncrlist),dtype=float)
    P=P.tolist()
    one=pylab.ones(len(ncrlist),dtype=float)
    for nid in nodes:
        p=CondProb(N[nid].rv[0].values(),samples,ncrlist)
        if len(no_broken_particles) == 0:
            P=JVT.Plus(P,JVT.Star(p,JVT.ScalarMult(one,1.0/numnodes)))
        else:
            P=JVT.Plus(P,JVT.Star(p,JVT.ScalarMult(one,\
                       float(N[nid].samples)/float(total_broken_parts))))

    # compute reliability
    R=JVT.Minus(one,P)

    pylab.axes([0.15,0.2,0.7,0.7])
    pylab.plot(ncrlist,R,color,linewidth=2,label=lab)
    pylab.xticks(size=18)
    pylab.yticks(size=18)

    return P,R,copy.copy(ncrlist)

########################################################

def PlotWeibull(R,ncrlist,color,lab,Nstar):

    # can't log(zero)... so log(log(1.0)) has to be avoided
    a=Numeric.less(R,1.0)
    a=a.tolist()
    index=a.index(1)
    R=R[index:]
    ncrlist=ncrlist[index:]

    # can't log(zero)... so remove zero reliability
    a=Numeric.greater(R,0.)
    a=a.tolist()
    try:
        index=a.index(0)
        R=R[:index]
        ncrlist=ncrlist[:index]
    except ValueError: pass 

    # Weibull paper
    xo=1.0
    P=[math.log(R[i]) for i in range(len(R))]
    PPP=[math.log(-1.0*P[i]) for i in range(len(P))]
    x=[math.log(ncrlist[i]/Nstar) for i in range(len(ncrlist))]

    pylab.axes([0.15,0.2,0.7,0.7])
    pylab.plot(x,PPP,color,linewidth=2,label=lab)
    pylab.xticks(size=18)
    pylab.yticks(size=18)
    pylab.xticks(pylab.arange(10.4,11.8,0.2))

########################################################

def PlotHistOfA(ai,color,norm):
    pylab.axes([0.15,0.2,0.7,0.7])
    pylab.hist(ai.rv[0].keys(),normed=norm,align='center',bins=100,fc=color)
    pylab.xticks(size=18,rotation=45)
    pylab.yticks(size=18)

########################################################

def Populaternd(rndname):

    f = open(rndname,'r')
    b=f.readlines()
    ai = Statistic.Stata(1,len(b))
    ai.UpdateSet()
    for line in b:
        stuff=string.splitfields(line)
        if stuff: ai.UpdateSamples(float(stuff[1]),0,int(stuff[0]))

    return ai

########################################################

def PopulateN(Nname):
    '''
    The Nfile is three columns and contains the data to map a damage origin
    i.d. to the random initial crack size id to the life prediction for it. The
    file format is: 

    doid rid N_rid

    rid = -1 if no cracks grew there.  
    '''

    f = open(Nname,'r')
    b=f.readlines()
    N = {}; no_broken_particles=[]
    for line in b:
        stuff=string.splitfields(line)
        if stuff: 
            if N.has_key(int(stuff[0])):
                N[int(stuff[0])].UpdateSamples(float(stuff[2]),0,int(stuff[1]))
            else:
                N[int(stuff[0])]=Statistic.StatN(1,'dynamic')
                N[int(stuff[0])].UpdateSet()
                N[int(stuff[0])].UpdateSamples(float(stuff[2]),0,int(stuff[1]))

                if stuff[1]=='-1': no_broken_particles+=[int(stuff[0])]

    return N,no_broken_particles

########################################################

def HelpPrints():
    print " " 
    print " code to plot the cummulative probability of failure & a histogram"
    print " of the initial particles (flaw) radii"  
    print " ---- " 
    print " Arguments are:"
    print " ...>Prob_plot.py -rnd [list of rnd filenames] -N [list N filenames]"
    print "  *** the way i wrote it, the order matters... do -rnd first!"
    print " "
    print " rnd file has format:"
    print "     rid ai"
    print " N file has format:"
    print "     nid rid N"
    print "      --> if rid = -1, no cracks grew."

    sys.exit()

########################################################

if __name__=="__main__":

    # a helper for the keys...
    if "-help" in sys.argv or "-h" in sys.argv:
        HelpPrints()

    try:
        rndindex=sys.argv.index("-rnd")
        Nindex=sys.argv.index("-N")

        rndfilenames=sys.argv[rndindex+1:Nindex]
        Nfilenames=sys.argv[Nindex+1:]

    except: HelpPrints()

    print ' '

    Nmax=1e5
    Nstar=1.0 # for weibull x-axis --> Ln(x/x*) 
    colorlinestyle=['k','k-.','g','b']
    color=['g','r','b','k']
##    labels=['10','100','1,000','10,000']
    labels=['Part. Fltr.','Weibull']
    # norm = 1 normalized histogram so like a pdf... if desperately different
    # number of samples in data you're comparing, DON'T normalize... 
    norm = 1
    ncrlist=pylab.arange(0.,Nmax*1.01,Nmax*0.01) # good1
##    ncrlist=pylab.arange(36000.0,42000.0,20.0)

    for i in range(len(rndfilenames)):
        rndname=rndfilenames[i]
        Nname=Nfilenames[i]
        col=color[i]
        colline=colorlinestyle[i]
        label=labels[i]

        print rndname,Nname

        ai = Populaternd(rndname)
        samples=len(ai.rid_ai.keys())
        samples = 10000
        N,no_broken_particles = PopulateN(Nname) 
        nids=N.keys()

        # i want to compute total probs correctly.  for particle filter the
        # prob of occurence of a flaw is:
        # (# broken particles at node i)/(total # broken particles)
        # for weibull, who knows what it is, so we use 1/num_nodes.  if
        # my list of no_broken_particles is empty, then i can account for this
        # in PlotReliability and PlotWeibull... 
        if label=='Weibull': no_broken_particles=[]

        pylab.figure(1)
        PlotHistOfA(ai,col,norm)

        pylab.figure(2)
        P,R,Currentncrlist=PlotReliability(N,samples,Nmax,colline,label,ncrlist,no_broken_particles)

        pylab.figure(3)
        PlotWeibull(copy.copy(R),copy.copy(Currentncrlist),colline,label,Nstar)

        del(ai)
        del(N)
        del(no_broken_particles)

    # i want this for figure(2)... which is no longer current.
    pylab.figure(2)
##    pylab.xlim0.0,Nmax*1.01)) # good1
##    pylab.ylim((((0.975,1.002)) # good1
##    pylab.xlim((34800,36200))
    pylab.ylim((0,1.01))
    pylab.legend(loc='lower left',\
                 prop=pylab.matplotlib.font_manager.FontProperties(size=18))
    pylab.xlabel(r'$\rm{Critical}\ \rm{load}\ \rm{cycles,}\ n_{cr}$',size=18)
    pylab.ylabel(r'$\rm{Reliability,}\ P(N>n_{cr})$',size=18)

    pylab.figure(3)
    pylab.legend(loc='lower right',\
                 prop=pylab.matplotlib.font_manager.FontProperties(size=18))
    pylab.xlabel(r'$\rm{Ln}(n_{cr})$',size=18)
    pylab.ylabel(r'$\rm{Ln}(-\rm{Ln}(P(N>n_{cr})))$',size=18)

##    pylab.figure(1)
##    pylab.savefig('HistA.eps')
##    os.system('epstopdf HistA.eps')
##
##    pylab.figure(2)
##    pylab.savefig('Reliab_part_v_weib_119032.eps')
##    os.system('epstopdf Reliab_part_v_weib_119032.eps')
##
##    pylab.figure(3)
##    pylab.savefig('Weibull_part_v_weib_119032.eps')
##    os.system('epstopdf Weibull_part_v_weib_119032.eps')

    pylab.show()





















