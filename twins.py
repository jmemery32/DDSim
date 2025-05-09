import cPickle, os, sys, string, Statistic

def Ansys2Sig(ansysfile,sigfile,confile):
    '''
    converst some ansys output of stresses to a *.sig file that MeshTools can
    read.
    '''

    filelines=ansysfile.readlines()
    # build a dictionary of the corners stresses to linearly interpolate for
    # midside nodes below.
    corners={}
    for line in filelines:
        try: nid=int(line[0:8])
        except IndexError:
            continue
        except ValueError:
            continue
        sx=float(line[9:21])
        sy=float(line[21:33])
        sz=float(line[33:45])
        sxy=float(line[45:57])
        syz=float(line[57:69])
        sxz=float(line[69:81])
        sigfile.write(str(nid-1)+' '+str(sx)+' '+str(sy)+' '+str(sz)+' '+\
                      str(sxy)+' '+str(syz)+' '+str(sxz)+"\n")
        corners[nid-1]=[sx,sy,sz,sxy,syz,sxz]

    # midside nodes
    conlines=confile.readlines()
    for sline in conlines:
        line=string.splitfields(sline)
        midsideid=int(line[0])
        aid=int(line[3])
        bid=int(line[4])
        sx=0.5*(corners[aid][0]+corners[bid][0])
        sy=0.5*(corners[aid][1]+corners[bid][1])
        sz=0.5*(corners[aid][2]+corners[bid][2])
        sxy=0.5*(corners[aid][3]+corners[bid][3])
        syz=0.5*(corners[aid][4]+corners[bid][4])
        sxz=0.5*(corners[aid][5]+corners[bid][5])
        sigfile.write(str(nid)+' '+str(sx)+' '+str(sy)+' '+str(sz)+' '+\
                      str(sxy)+' '+str(syz)+' '+str(sxz)+"\n")

def CutMAPFile(oldMAPFile,newMAPFile,setid,high,low):
    '''
    Writes file for contouring in MAP with upper and lower values cutoff

    DamOro - Dictionary with key = node id, value = Statistic instance of life
             (SampMean, SampVar already populated).
    MAPFile - file object to write to
    set  - the set from which to write mean value for
    '''

    newMAPFile.write('LIFE 0'+"\n")
    #throw away the first line of oldMAPFile
    junk=oldMAPFile.readline()
    myfile = oldMAPFile.readlines()
    for line in myfile:
        x=string.splitfields(line)
        id=int(x[0])
        value=float(x[1])
        if value > high:
            value = high
        if value < low:
            value = low
        newMAPFile.write(str(id)+' '+str(value)+"\n")

def MakeDamOro(Nfile,buff,doid,setid):

    DamOro=Statistic.StatN(1,'dynamic')
    DamOro.UpdateSet()

    while buff:
        line=string.split(buff)
        newdoid=int(line[0])
        if newdoid==doid: DamOro.UpdateSamples(float(line[2]),0,int(line[1]))
        else: break
        buff=Nfile.readline()

    return DamOro,buff

def ToMAPFile(DamOro,MAPFile,set,HighLife):
    '''
    Writes file for contouring in MAP

    DamOro - Dictionary with key = node id, value = Statistic instance of life
             (SampMean, SampVar already populated).
    MAPFile - file object to write to
    set  - the set from which to write mean value for
    '''

    DamKeys = DamOro.keys()
    DamKeys.sort()

    MAPFile.write('LIFE 0'+"\n")
    for i in DamKeys: # i = doid
        if DamOro[i].SampleMean(set) == -1.0:
            MAPFile.write(str(i)+' '+str(HighLife)+"\n")
        elif len(DamOro[i].rv[set])==1 and DamOro[i].rv[set].has_key(-1):
            MAPFile.write(str(i)+' '+str(HighLife)+"\n")
        else:
            MAPFile.write(str(i)+' '+str(DamOro[i].SampleMean(set))+"\n")
##        if DamOro[i].rv[0][0] == -1.0:
##            MAPFile.write(str(i)+' '+str(HighLife)+"\n")
##        else:
##            MAPFile.write(str(i)+' '+str(DamOro[i].rv[0][0])+"\n")

def ProbFail(DamOro,MAPfile,set,ncr):
    '''
    write a MAP contour file that plots the probability of failure given a
    critical number of cycles, ncr.
    '''

    MAPFile.write("Pf_ncr 0"+"\n")
    for i in DamOro: # i = doid
        P = DamOro[i].CalcProb(set,ncr)
        MAPFile.write(str(i)+' '+str(P)+"\n")
    

def JoinPickles(newfile,num):
    '''
    Function to join pickled data from DamModel.  Assumes they are pickled in 
    binary mode.

    num - number of files being joined, usually # of nodes used in parallel 
    '''

    mydict={}
    for i in xrange(num):
        fn = open(filename+'.'+str(i),'r+b')
        addtodict = cPickle.load(fn)
        fn.close()
        mydict.update(addtodict)

    cPickle.dump(mydict,newfile,2)

def JoinASCII(f_joined,num):
    '''
    Function to join ASCII output files
    '''

    for i in xrange(num):
        f1 = open(filename+'.'+str(i),'r')
        myfile = f1.readlines()
        f1.close()
        for line in myfile:
            f_joined.write(line)

def JoinGerdASCII(f_joined,num):
    '''
    Function to join ASCII output files in the form of Gerd's new output.  He
    starts his line with <load case> <load step> (or something like that)
    '''

    f_joined.write('Stress 2'+"\n")
    for i in xrange(num):
        f1 = open(filename+'.'+str(i),'r')
        myfile = f1.readlines()
        f1.close()
        for line in myfile:
            x=string.splitfields(line)
##            x.pop(0)
##            x.pop(0)
            line=string.join(x)
            f_joined.write(line+"\n")

def PrintHelps():
    print ' '
    print ' -file <file.name>'
    print ' -num <number of processors>'
    print ' -map <integer> provied MAP contour file for given file (above) and'
    print "      set.  Actually creates from pickled N's. "
    print ' -P or -p if they are pickled files'
    print ' -A or -a if they are ASCII files'
    print " -G or -g if they are ASCII files in Gerd's new nodal stress format"
    print " -c or -C to cut chop values in a *.tab.0 file followed by:"
    print " -setid <integer> set i.d. (for cutoff map file)"
    print " -high <number> is the upper value to cut contours off at"
    print " -low <number> is the lower value to cut contours off at"
    print " -ansys to convert and ansys stress output file to *.sig (readable "
    print "        by Meshtools.pyd).  Specify ansys file with -file (above) "
    print "        but do NOT include file extension!!!!! name of ansys file"
    print "        should be filename.str, must have filename.con in same "
    print "        working directory and will get filename.sig."
    print " -report <BaseFileName> - use to report performance information.  "
    print "          Assumes already cat'd BaseFileName.stn.$MSTI_RANK$ files."
    print ' '

def ReadNFile(Nfile):
    '''
    Read the Nfile which looks like this:
            doid rid N
            doid rid N
            doid rid N...

    and store it in DamOro, a dictionary that looks like this:
            DamOro={doid:StatN,doid:StatN...}

    Currently assumes only one set!  
    '''

    LinesInNfile = Nfile.readlines()

    DamOro={}
    for line in LinesInNfile:
        sortedline=string.split(line)
        doid=int(sortedline[0])
        rid=int(sortedline[1])
        N=float(sortedline[2])
        if DamOro.has_key(doid): DamOro[doid].UpdateSamples(N,0,rid)
        else:
            DamOro[doid]=Statistic.StatN(1,'dynamic')
            DamOro[doid].UpdateSet()
            DamOro[doid].UpdateSamples(N,0,rid)

    return DamOro

def Report(report):
    Nfile=open(report+".stn",'r+b')
    DamOro=cPickle.load(Nfile)
    Nfile.close()
    HighLife=0.0
    print ' **Assuming there is at least one node that processed to N_max.'
    for i in DamOro:
        try:
            if DamOro[i].SampleMean(setid) > HighLife: 
                HighLife = DamOro[i].SampleMean(setid)
        except IndexError:
            print i
            continue

    NumNodes=len(DamOro)
    NumSets=DamOro[0].sets
    NumSamples=DamOro[0].samples
    TotalKs=0
    for node in DamOro:
        for set in DamOro[node].rv:
            for N in set.values():
                if N == HighLife: TotalKs+=1
                else: TotalKs+=N

    Time=0
    blah=open(report+'.time','r')
    times=blah.readlines()
    for time in times:
        t=float(time)
        if t>Time: Time=t

    print 'Num nodes, Num Sets, Num Samples, Total Ks computed'
    print NumNodes, NumSets, NumSamples, TotalKs, Time

if __name__ == "__main__":
    '''
    To join data files generated during parallel processing (twins, or more,
    seperated at birth!).  Data can be either:
    1. Pickled dictionaries:
       key = integer (eg. doid)
       value = object (eg. N as a Statistic.Stat() instance)
    2. ASCII text files.

    Required input:
    filename - name of files to be joined that have *.MSTI_RANK extension (eg.
               example1.MP12) **** include extension!! 
    num      - number of nodes used in processing... so upper limit to
               MSTI_RANK
    -P       - indicates that the files to be joined are pickled lists
    '''

    # flag to cut contours 
    cut = 0
    if "-c" in sys.argv or "-C" in sys.argv: cut = 1
    high = 100000.0
    if "-high" in sys.argv:
        index = sys.argv.index('-high')+1
        high = float(sys.argv[index])
    low = 0.0
    if "-low" in sys.argv:
        index = sys.argv.index('-low')+1
        low = float(sys.argv[index])
    setid=0
    if "-setid" in sys.argv:
        index = sys.argv.index('-setid')+1
        setid = int(sys.argv[index])

    if "-h" in sys.argv or "-help" in sys.argv:
        PrintHelps()
        sys.exit()

    if "-file" in sys.argv:
        index = sys.argv.index('-file')+1
        filename = sys.argv[index]

    if "-num" in sys.argv:
        index = sys.argv.index('-num')+1
        num = int(sys.argv[index])
##    elif "-list" in sys.argv:
##        index = sys.argv.index('-list')+1
##        num = string.splitfields(sys.argv[index])
##    else: num=int(os.environ['procs'])

    nfile=0
    if "-ncr" in sys.argv:
        nfile=1
        index = sys.argv.index('-ncr')+1
        ncr = float(sys.argv[index])

    ansys=0
    if "-ansys" in sys.argv: ansys=1

    map = 0
    if "-map" in sys.argv:
        map = 1
        index = sys.argv.index('-map')+1
        setid = int(sys.argv[index])

    # flag to indicate if files to be joined are pickled lists
    pickles = 0 
    if "-P" in sys.argv or "-p" in sys.argv: pickles = 1

    # flag to indicate if files to be joined are ascii text files
    ascii =0 
    if "-A" in sys.argv or "-a" in sys.argv: ascii = 1

    # flag to indicate if files to be joined are ascii text files
    gerd =0 
    if "-G" in sys.argv or "-g" in sys.argv: gerd = 1

    # flag for report performance info
    report=0
    if "-report" in sys.argv:
        index = sys.argv.index('-report')+1
        report = sys.argv[index]

# ------------------------------
    if pickles:
        newfile=open(filename,'w+b')
        JoinPickles(newfile,num)
        newfile.close()

    # new with DDSimv1.5 --> no long using pickled N file.  Now it is ASCII in
    # in the format read by the database.  
    if map: 
        Nfile=open(filename,'r')
        MAPFile=open(filename+'.'+str(setid)+'.tab.0','w')
        print ' **Assuming there is at least one node that processed to N_max.'
        HighLife=0.0
        buff=Nfile.readline()

        while buff:
            line=string.split(buff)
            doid=int(line[0])
            DamOro,buff=MakeDamOro(Nfile,buff,doid,setid)

            if DamOro.SampleMean(setid) > HighLife:
                HighLife = DamOro.SampleMean(setid)

            if DamOro.SampleMean(setid) == -1.0:
                MAPFile.write(str(doid)+' HighLife'+"\n")
            elif len(DamOro.rv[setid])==1 and DamOro.rv[setid].has_key(-1):
                MAPFile.write(str(doid)+' HighLife'+"\n")
            else:
                MAPFile.write(str(doid)+' '+str(DamOro.SampleMean(setid))+"\n")

        MAPFile.close()
        Nfile.close()

        # to sort the file.
        MAPFile=open(filename+'.'+str(setid)+'.tab.0','r')
        s={}
        buff=MAPFile.readline()
        while buff:
            line=string.split(buff)
            s[int(line[0])]=line[1]
            buff = MAPFile.readline()
        MAPFile.close()

        MAPFile=open(filename+'.'+str(setid)+'.tab.0','w')
        MAPFile.write("Life 0 \n")
        nodes=s.keys()
        nodes.sort()
        for node in nodes:
            if s[node]=="HighLife": MAPFile.write(str(node)+' '+str(HighLife)+"\n")
            else: MAPFile.write(str(node)+' '+s[node]+"\n")
        MAPFile.close()

    if ascii:
        f_joined=open(filename,'w')
        JoinASCII(f_joined,num)
        f_joined.close()

    if gerd:
        f_joined=open(filename,'w')
        JoinGerdASCII(f_joined,num)
        f_joined.close()

    if cut:
##        oldMAPFile=open(filename+'.'+str(setid)+'.tab.0','r')
##        newMAPFile=open(filename+'cut.'+str(setid)+'.tab.0','w')
        oldMAPFile=open(filename+'.tab.0','r')
        newMAPFile=open(filename+'cut.tab.0','w')
        CutMAPFile(oldMAPFile,newMAPFile,setid,high,low)
        oldMAPFile.close()
        newMAPFile.close()

    if ansys:
        ansysfile=open(filename+'.str','r')
        sigfile=open(filename+'.sig','w')
        confile=open(filename+'.con','r')
        Ansys2Sig(ansysfile,sigfile,confile)
        ansysfile.close()
        sigfile.close()
        confile.close()

    if nfile:
        Nfile=open(filename,'r+b')
        DamOro=cPickle.load(Nfile)
        Nfile.close()
        MAPFile=open(filename+'.'+str(setid)+'.Pf.tab.0','w')
        ProbFail(DamOro,MAPFile,setid,ncr)
        MAPFile.close()

    if report:
        Report(report) 


