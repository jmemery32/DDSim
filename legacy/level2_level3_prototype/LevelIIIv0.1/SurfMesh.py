
import sys, MeshTools, Vec3D, ColTensor, string, os

class SurfaceMeshWriter:
    '''
    a class to extract the surface mesh from the MeshTools object and
    write the surfae mesh to the screen.

    LIMITATIONS:
    Currently i don't support multiple regions... this is only a temporary way
    of doing this anyway.  Gerd is going to write the surface mesh of
    polycrystals.  ha-ha-ha... when, exactly, is Gerd going to follow through
    with anything?  Let alone supply the surface mesh?  it's almost august... 
    '''

    def __init__(self,model):
        self.model=model
        self.SurfaceNodes = {} # key = nid, values = coordinates (no mids)
        self.SurfaceElements = {} # key = seid, values = corner nids, CCW
                                  # points out.  
        self.__MakeCMeshObj() # method to fill the above

    def __MakeCMeshObj(self):
        # loop over all nodes to populate self.SurfaceElements
        node_list=self.model.GetNodeList()
        for nid in node_list:
            if self.model.IsSurfaceNode(nid) == 1:
                self.SurfaceNodes[nid]=self.model.GetNodeInfo(nid)[0]
                list_seids = self.model.GetAdjacentSurfElems(nid)
                for seid in list_seids:
                    if self.SurfaceElements.has_key(seid):
                        continue
                    else:
                        SurfaceNidList = self.model.GetSurfElemInfo(seid)
                        if len(SurfaceNidList) == 6: # triangular surface facet
                            self.SurfaceElements[seid] = ((SurfaceNidList[0], \
                                                           SurfaceNidList[1], \
                                                           SurfaceNidList[2]),\
                                                          (SurfaceNidList[3], \
                                                           SurfaceNidList[4], \
                                                           SurfaceNidList[5]))
                        elif len(SurfaceNidList) == 8: # quad surface facet
                            self.SurfaceElements[seid] = ((SurfaceNidList[0], \
                                                           SurfaceNidList[1], \
                                                           SurfaceNidList[2], \
                                                           SurfaceNidList[3]),\
                                                          (SurfaceNidList[4], \
                                                           SurfaceNidList[5], \
                                                           SurfaceNidList[6], \
                                                           SurfaceNidList[7]))
                        elif len(SurfaceNidList) == 3: # linear tri surf facet
                            self.SurfaceElements[seid] = ((SurfaceNidList[0], \
                                                           SurfaceNidList[1], \
                                                           SurfaceNidList[2]),\
                                                          (None))
                        else: # linear quadrilateral surface facet
                            self.SurfaceElements[seid] = ((SurfaceNidList[0], \
                                                           SurfaceNidList[1], \
                                                           SurfaceNidList[2], \
                                                           SurfaceNidList[3]),\
                                                          (None))

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
            s+="VTX: "
            s+=str(vtx)+' '
            s+=str(self.SurfaceNodes[vtx].x())+' '
            s+=str(self.SurfaceNodes[vtx].y())+' '
            s+=str(self.SurfaceNodes[vtx].z())+"\n"

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

def HelpPrints():
    print " " 
    print " This script reads a finite element model, in RDB format, and"
    print " returns the surface mesh in Bruce's format.  The surface mesh"
    print " will print to the screen, so use pipe out to write to file"
    print " ---- " 
    print " Arguments are:"
    print " ...>SurfMesh.py [input filename (no extension)] "
    print ' '

if __name__=="__main__":

    # a helper for the keys...
    if "-help" in sys.argv or "-h" in sys.argv:
        HelpPrints()
        sys.exit()

    try:
        filename=sys.argv[1]
    except IndexError:
        HelpPrints()
        sys.exit()

    #
    try:
        blah=open(filename+'.sig','r')
    except IOError:
        blah=open(filename+'.nod','r')
        buff=blah.readline()
        s=''
        while buff:
            line=string.splitfields(buff)
            s+=line[0]+" 0 0 0 0 0 0 \n"
            buff=blah.readline()
        s=s[0:-1]
        blah.close()
        sig=open(filename+'.sig','w')
        sig.write(s)
        sig.close()

    try:
        blah=open(filename+'.smp','r')
    except IOError:
        blah=open(filename+".ShapeMap",'r')
        blah.close()
        os.system("copy /Y "+filename+".ShapeMap "+filename+".smp")
    except IOError:
        print "you also need the ShapeMap file for this model:"
        print "  ", filename+".smp"
        print "  exiting..."
        sys.exit()

    # read and store the finite element model
    model = MeshTools.MeshTools(filename,'RDB')
    model.SetPointInsideTolerance(1.0e-7)

    # store object
    writer = SurfaceMeshWriter(model)

    # write surface mesh file to screen
    writer.WriteIt()
