# quick hack to get a surface mesh from some file Gerd wrote for me
# hopefully i only have to do this once! 

import string, Vec3D

class surface:

    def __init__(self,quads,tris,nodes):

        self.SurfaceNodes={}
        self.Surf2Elem={} 
        self.SurfaceElements={}

        buff=quads.readline()
        while buff:
            line=string.split(buff,',')
            buff=quads.readline()
            qid=int(line[0])
            v0=int(line[1])
            if not self.Surf2Elem.has_key(v0): self.Surf2Elem[v0]=qid
            v1=int(line[2])
            if not self.Surf2Elem.has_key(v1): self.Surf2Elem[v1]=qid
            v2=int(line[3])
            if not self.Surf2Elem.has_key(v2): self.Surf2Elem[v2]=qid
            v3=int(line[4])
            if not self.Surf2Elem.has_key(v3): self.Surf2Elem[v3]=qid

            self.SurfaceElements[qid]=[[v0,v1,v2,v3],[]]

##        print self.SurfaceElements

        buff=tris.readline()
        while buff:
            line=string.split(buff,',')
            buff=tris.readline()
            qid=int(line[0])
            v0=int(line[1])
            if not self.Surf2Elem.has_key(v0): self.Surf2Elem[v0]=qid
            v1=int(line[2])
            if not self.Surf2Elem.has_key(v1): self.Surf2Elem[v1]=qid
            v2=int(line[3])
            if not self.Surf2Elem.has_key(v2): self.Surf2Elem[v2]=qid

            self.SurfaceElements[qid]=[[v0,v1,v2],[]]


        buff=nodes.readline()
        while buff:
            line=string.split(buff)
            buff=nodes.readline()
            nid=int(line[0])
            if self.Surf2Elem.has_key(nid):
                self.SurfaceNodes[nid]=Vec3D.Vec3D(float(line[1]),float(line[2]),float(line[3]))




    def WriteIt(self):
        '''
        go through the surface mesh and output to screen the appropriate text:

        VTX: id x y z
        FACET: id # A B C region0 region1 retained
        REGION: 0 0

        facets are never "retained", hence retained=0.
        '''

        vtxs=self.SurfaceNodes.keys()
        vtxs.sort()
        s=str()
        for vtx in vtxs:
            try:
                s+="VTX: "
                s+=str(vtx)+' '
                s+=str(self.SurfaceNodes[vtx].x())+' '
                s+=str(self.SurfaceNodes[vtx].y())+' '
                s+=str(self.SurfaceNodes[vtx].z())+"\n"
            except AttributeError:
                print vtx, self.SurfaceNodes[vtx]

        facets=self.SurfaceElements.keys()
        facets.sort()
        for facet in facets:
            s+="FACET: "+str(facet)+' '
            verts=self.SurfaceElements[facet]
            if len(verts[0]) == 3:
                a='3 '
                b=str(verts[0][0])+' '+str(verts[0][1])
                b+=' '+str(verts[0][2])
                if verts[1]:
                    a='6 '
                    b+=' '+str(verts[1][0])+' '+str(verts[1][1])
                    b+=' '+str(verts[1][2])
                c=" -1 0 0 \n"
                s+=a+b+c
            else: # i.e. if len(verts[0]) == 4:
                a='4 '
                b=str(verts[0][0])+' '+str(verts[0][1])
                b+=' '+str(verts[0][2])+' '+str(verts[0][3])
                if verts[1]:
                    a='8 '
                    b+=' '+str(verts[1][0])+' '+str(verts[1][1])
                    b+=' '+str(verts[1][2])+' '+str(verts[1][3])
                c=" -1 0 0 \n"
                s+=a+b+c
                
        s+="REGION: 0 0"

        print s

################################################

def WriteOut(fi,fo,surf):
    buff=fi.readline()
    while buff:
        line=string.split(buff)
        buff=fi.readline()
        if surf.SurfaceNodes.has_key(int(line[0])): fo.write(buff)

################################################

if __name__=="__main__":

    quads=open('SurfaceQuadrilateralVertices.csv','r')
    tris=open('SurfaceTriangleVertices.csv','r')
    nodes=open('Nodes','r')

    surf=surface(quads,tris,nodes)

    quads.close()
    tris.close()
    nodes.close()
    
##    surf.WriteIt()

##    disp=open('h:\\users\\jme32\\research\\Dissertation\\DDSimLIII\\SIPS3002_refined\\Plus_70GrainLumber\\LEI_Displacement\\Displacement.L00','r')
##    fo=open('h:\\users\\jme32\\research\\Dissertation\\DDSimLIII\\SIPS3002_refined\\Plus_70GrainLumber\\LEI_Displacement\\Displacement.L00.tab.0','w')
##    fo.write("Displacement 1 \n")

##    WriteOut(disp,fo,surf)

##    disp.close()
##    fo.close()

    sig=open('h:\\users\\jme32\\research\\Dissertation\\DDSimLIII\\SIPS3002_refined\\Plus_70GrainLumber\\LEI_Displacement\\NodalStress.L00','r')
    fo=open('h:\\users\\jme32\\research\\Dissertation\\DDSimLIII\\SIPS3002_refined\\Plus_70GrainLumber\\LEI_Displacement\\NodalStress.L00.tab.0','w')
    fo.write("Stress 2 \n")

    WriteOut(sig,fo,surf)

    sig.close()
    fo.close()
