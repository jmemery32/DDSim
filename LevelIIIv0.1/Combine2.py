
import sys, string, os

class MPCFileJoiner:
    '''
    a class to operate on global and local models to renumber the global model
    and glue them into a set of files to run with cptc.
    '''

    def __init__(self,Gfile,Lfile):
        '''
        1.  get increment by reading global files and storing 1 + max element
            id, node id, material id, and ShapeMap id
        2.  write a new file, xGfile, and add the increment to each of the
            global files.
        3.  glue them together using, e.g.:
            os.system('copy /y Lfile.con+xGfile.con LfileGfile.con')
        4.  copy to Gerd's file naming using, e.g.:
            os.system('copy /y LfileGfile.con Elements')
        '''

        self.G=Gfile # string obj that contains the global file's first name
        self.L=Lfile # "      "   "    "        "   local  "      "     "

        # 1. get increments...
        fin = open(Lfile+'.nod','r')
        self.NodeInc = self.Search(fin,0,0)+1
        print "Node increment is: ", self.NodeInc
        fin.close()

        # add extra line at end of file so when we do copy + it works... 
        fin = open(Lfile+'.nod','a')
        fin.write("\n")
        fin.close()

        fin = open(Lfile+'.con','r')
        self.ElemInc = self.Search(fin,0,0)+1
        print "Element increment is: ", self.ElemInc
        fin.close()

        fin = open(Lfile+'.Materials','r')
        self.MatInc = self.Search(fin,0,0)+1
        print "Material increment is: ", self.MatInc
        fin.close()

        fin = open(Lfile+'.ShapeMap','r')
        self.ShapeInc = self.Search(fin,0,0)+1
        print "Shape increment is: ", self.ShapeInc
        fin.close()

        # 2. write new file, incrementing the global file

        xfG = open('x'+Gfile+'.nod','w')
        fG  = open(Gfile+'.nod','r')
        self.Append(xfG,fG,self.NodeInc,0,0)
        xfG.close()
        fG.close()

        print 'xnodes done'

        xfG = open('x'+Gfile+'.con','w')
        fG  = open(Gfile+'.con','r')
        self.Append(xfG,fG,self.ElemInc,4,self.NodeInc)
        xfG.close()
        fG.close()

        print 'xcon done'

        # addjust the Matieral id and Shape Model id in the Global model
        self.Adder('x'+Gfile+'.con','int',self.MatInc,2)
        print 'xcon material update done'
        self.Adder('x'+Gfile+'.con','int',self.ShapeInc,1)
        print 'xcon ShapeFunc id update done'

        try:
            xfG = open('x'+Gfile+'.edg','w')
            fG  = open(Gfile+'.edg','r')
            self.Append(xfG,fG,self.NodeInc,3,self.NodeInc)
            fin.close()
            fG.close()
        except IOError: pass 

        print 'xedg done'

        try: 
            xfG = open('x'+Gfile+'.sig','w')
            fG  = open(Gfile+'.sig','r')
            self.Append(xfG,fG,self.NodeInc,0,0)
            fin.close()
            fG.close()
        except IOError: pass 

        print 'xsig done'

        try:
            xfG = open('x'+Gfile+'.BC1','w')
            fG  = open(Gfile+'.BC1','r')
            self.BC(xfG,fG,self.NodeInc)
            fin.close()
            fG.close()
        except IOError: pass

        print 'xBC1 done'

        try:
            xfG = open('x'+Gfile+'.ShapeMap','w')
            fG  = open(Gfile+'.ShapeMap','r')
            self.BC(xfG,fG,self.ShapeInc)
            fin.close()
            fG.close()
        except IOError: pass 

        print 'xShapeMap done'

        try:
            xfG = open('x'+Gfile+'.Materials','w')
            fG  = open(Gfile+'.Materials','r')
            self.BC(xfG,fG,self.MatInc)
            fin.close()
            fG.close()
        except IOError: pass 

        print "xMaterials done"

        try:
            xfG = open('x'+Gfile+'.BC3','w')
            fG  = open(Gfile+'.BC3','r')
            self.BC(xfG,fG,self.NodeInc)
            fin.close()
            fG.close()
        except IOError: pass 

        print 'xBC3 done'

        try:
            xfG = open('x'+Gfile+'.SrfCVecElems','w')
            fG  = open(Gfile+    '.SrfCVecElems','r')
            self.BC(xfG,fG,self.ElementInc)
            fin.close()
            fG.close()
        except IOError: pass

        print 'xSrfCVecElems done'

        try:
            xfG = open('x'+Gfile+'.NormSrfCScalElems','w')
            fG  = open(Gfile+'.NormSrfCScalElems','r')
            self.BC(xfG,fG,self.ElementInc)
            fin.close()
            fG.close()
            # because we glue them together in Glue()... 
            try:
                check=open(Lfile+'.NormSrfCScalElems','r')
                check.close()
            except IOError:
                make=open(Lfile+'.NormSrfCScalElems','w')
                make.close() 
        except IOError: pass 

        print 'xNormSrfCScalElems done'

        try:
            # first update the mpc file
            xfG = open('MPC','w')
            fin = open(Lfile+'.mpc','r')
            self.MPC(fin,xfG,self.NodeInc)
            fin.close()

            # now append with BC file 
            fuse = open('x'+Gfile+'.BC1','r')
            self.UpdateBC1(xfG,fuse)
            xfG.close()
            fuse.close()
        except IOError: pass

        print 'MPC done'

        try:
            # add local model BC1 if it exists
            localBC=open(Lfile+'.BC1','r')
            MPC=open('MPC',"a")
            self.LocalBC1(localBC,MPC)
            localBC.close()
            MPC.close()
        except IOError: pass 

        # 3. Glue together
        self.GlueMe()

################

##    def UpdateMat_Shape(self,infile):
##        con=open(infile,'r')
##        filelines=con.readlines()
##        con.close()
##
##        out=open(infile,'w')
##        line in filelines:
##            h=string.split(line)
            

################

    def Adder(self,filename,addertype,adder,column):
        blah=open(filename,'r')
        buff=blah.readline()
        s=''
        while buff:
            line=string.splitfields(buff)
            buff=blah.readline()

            if addertype=="float": new=float(line[column])+adder
            elif addertype=="int": new=int(line[column])+adder
            line[column]=str(new)
            for l in line: s+=l+' '
            s=s[0:-1]+"\n"

        blah.close()
        newblah=open(filename,'w')
        newblah.write(s[0:-1])
        newblah.close()

################

    def Append(self,afile,second,cinc,columid,inc):
        '''
        append to afile with second file incrementing
        the first column by cinc and columdid to the end of the row by inc. if
        columid == 0, only increment first column by cinc. 
        '''

        buff=second.readline()
        while buff:
            list=string.split(buff)
            buff=second.readline()

            list[0]=int(list[0])+cinc
            if columid:
                list[columid:] = \
                        [int(list[i])+inc for i in range(columid,len(list))]

            s=''
            for p in list: s+=str(p)+' '
            s=s[0:-1]+"\n"

            afile.write(s)

################

    def Search(self,fin,colid,type):
        '''
        search the file row by row and return the highest "type" in column
        "colid".

        type = 0 --> integer
        type = 1 --> float
        '''

        old_item=-1
        buff=fin.readline()
        while buff:
            list=string.split(buff)
            buff=fin.readline()

            if type: item=float(list[colid])
            else: item=int(list[colid])

            if item > old_item:
                old_item = item

        return old_item

################

    def BC(self,afile,second,cinc):
        '''
        append to afile (use open('blah','a')) with second file incrementing
        the first column by cinc and columdid to the end of the row by inc. if
        columid == 0, only increment first column by cinc. 
        '''

        buff=second.readline()
        while buff:
            list=string.split(buff)
            buff=second.readline()

            list[0]=int(list[0])+cinc

            s=''
            for p in list: s+=str(p)+' '
            s=s[0:-1]+"\n"

            afile.write(s)

################

    def MPC(self,fin,fout,inc):
        '''
        change node id's of global nodes in mpc file
        '''

        # the column in the MPC file that tells how many mpc's are coming...
        # currently column #3... but file formats are always subject to change.
        numCol=3 

        buff=fin.readline()
        while buff:
            list=string.split(buff)
            buff=fin.readline()

            number=int(list[numCol])
            for i in range(number):
                list[numCol+1+(2*i)] = int(list[numCol+1+(2*i)]) + inc

##            # add missing parameter... JD did fix...
##            list.insert(2,'0')

            s=''
            for p in list: s+=str(p)+' '

            # strip off the last space
            s=s[0:-1]+"\n"

            fout.write(s)

################

    def LocalBC1(self,localBC,MPC):
        '''
        update the MPC file for any local BC1's... no renumbering necessary.  
        '''

        buff=localBC.readline()
        while buff:
            list=string.split(buff)
            buff=localBC.readline()

            # you need the following... even when you forget you do...
            list.insert(2,'1')

            s=''
            for p in list: s+=str(p)+' '
            s=s[0:-1]+" 0 0\n"

            MPC.write(s)


################

    def UpdateBC1(self,fa,ftran):
        '''
        update the type 1 BC's in ftran.  
        '''

        buff=ftran.readline()
        while buff:
            list=string.split(buff)
            buff=ftran.readline()

            # you need the following... even when you forget you do...
            list.insert(2,'1')

            s=''
            for p in list: s+=str(p)+' '
            s=s[0:-1]+" 0 0\n"

            fa.write(s)

################

    def GlueMe(self):
        '''
        glue them together using, e.g.:
            os.system('copy /y Lfile.con+xGfile.con LfileGfile.con')
        '''

        os.system('mkdir xFEAWDx')

        osString = 'copy /Y '+self.L+".BC1+x"+self.G+".BC1"+ \
                   ' '+self.G+self.L+".BC1"
        os.system(osString)

        osString = 'copy /Y '+self.L+".con+x"+self.G+".con"+ \
                   ' '+self.G+self.L+".con"
        os.system(osString)
        osString = 'copy /Y '+self.G+self.L+".con xFEAWDx\\Elements"
        os.system(osString)

        osString = 'copy /Y '+self.L+".edg+x"+self.G+".edg"+ \
                   ' '+self.G+self.L+".edg"
        os.system(osString)
        osString = 'copy /Y '+self.G+self.L+".edg xFEAWDx\\Components"
        os.system(osString)

        osString = 'copy /Y '+self.L+".Materials+x"+self.G+".Materials"+ \
                   ' '+self.G+self.L+".Materials"
        os.system(osString)
        osString = 'copy /Y '+self.G+self.L+".Materials xFEAWDx\\Materials"
        os.system(osString)

        osString = 'copy /Y '+self.L+".nod+x"+self.G+".nod"+ \
                   ' '+self.G+self.L+".nod"
        os.system(osString)
        osString = 'copy /Y '+self.G+self.L+".nod xFEAWDx\\Nodes"
        os.system(osString)

        osString = 'copy /Y '+self.L+".NormSrfCScalElems+x"+self.G+".NormSrfCScalElems"+ \
                   ' '+self.G+self.L+".NormSrfCScalElems"
        os.system(osString)
        osString = 'copy /Y '+self.G+self.L+".NormSrfCScalElems xFEAWDx\\NormSrfCScalElems"
        os.system(osString)

        osString = 'copy /Y '+self.L+".ShapeMap+x"+self.G+".ShapeMap"+ \
                   ' '+self.G+self.L+".ShapeMap"
        os.system(osString)
        osString = 'copy /Y '+self.G+self.L+".ShapeMap xFEAWDx\\ShapeMap"
        os.system(osString)

        os.system('copy /Y MPC xFEAWDx\\MPC')

        os.system('del x'+self.G+".BC1")
        os.system('del x'+self.G+".con")
        os.system('del x'+self.G+".edg")
        os.system('del x'+self.G+".Materials")
        os.system('del x'+self.G+".nod")
        os.system('del x'+self.G+".NormSrfCScalElems")
        os.system('del x'+self.G+".ShapeMap")
        os.system('del x'+self.G+".SrfCVecElems")
        os.system('del x'+self.G+".sig")
        os.system('del x'+self.G+".BC3")

        print "I'm not doing the MaterialModelMap file... just do it yourself"
        print " oh, and by the way, you need to switch the Components file column."

################

def HelpPrints():
    print " " 
    print " This script reads two finite element models, the Global and the "
    print " Local, in RDB format and the Local.mpc file.  Then combines them "
    print " by determining the highest node and element number in the local "
    print " model and incrementing all numbering in the global model by as "
    print " much.  "
    print " ---- " 
    print " Usage:"
    print " C:\...>Combine.py [global filename] [local filename] (no extensions)] "
    print " "
    print " Output:"
    print " Global+Local.con, .nod, etc."
    print ' '

################

if __name__=="__main__":

    print "  ************************************************************* "
    print " * still working the kinks out..." 
    print "  ************************************************************* "

    # a helper for the keys...
    if "-help" in sys.argv or "-h" in sys.argv:
        HelpPrints()
        sys.exit()

    try:
        Globalfilename=sys.argv[1]
        Localfilename=sys.argv[2]
    except IndexError:
        HelpPrints()
        sys.exit()

    # do the merging #######
    merger = MPCFileJoiner(Globalfilename,Localfilename)
    # ####### 




