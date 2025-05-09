import Vec3D

#  Heavily modified for DDSim V1.3.  JME 12/8/05

class SurfaceMesh:
    """
    this class stores the connectivity of a surface mesh and
    provides some tools to use with the surface mesh
    """

    def __init__(self):
        self.__elems = {} #values are element classes below
        self.__seid = -1 # surface element id start at 0 and keep udating,
                        # initialized to -1 to indicate no elements exist yet.
        self.__revcon = {} # "reverse connectivity" i.e. key = node; value =
                           # attached elements.
##        self.__norms = {} # a dictionary the contains
        self.__tris = {} # a dictionary of all triangular facets (not just
                         # surface facets)
        self.__quads = {} # a dictionary of all quadriateral facets
        self.__Nodes = {} # a dictionary of nodes (key) and position (values)
        self.maxX = 0.0 # to create a bounding box for the surface mesh.
        self.maxY = 0.0
        self.maxZ = 0.0
        self.minX = 0.0
        self.minY = 0.0
        self.minZ = 0.0

# ------ SurfaceMesh --------

    def AddNode(self,id,xyz):
        '''
        id = integer node id
        xyz = Vec3D object of nodal position
        '''
        self.__nodes[id] = xyz

# ------ SurfaceMesh --------

    def AddSurfEl(self,id,surf):
        '''
        id - integer element id
        surf - List that contains the nodes of the element; corner nodes
         preceeding midside nodes.
        '''

        # Linear or quadratic triangles
        if len(surf) == 6 or len(surf) == 3:
            self.AddTri(id,surf) # adds a triangular surface element and
                              # the reverse connetivity of the triangle.
        else:
            self.AddQuad(id,surf) # ditto for quads.

# ------ SurfaceMesh --------

    def AddTri(self,id,surf):
        '''
        method to add a triangle object to the __tris dictionary.
        '''
        SurfElem = triangle(surf[0:3],surf[3:])
        self.__tris[id]=SurfElem
        self.__BuildRevCon(id,SurfElem)

# ------ SurfaceMesh --------

    def AddQuad(self,surf):
        '''
        method to add a quadrilateral object to the __quads dictionary.
        '''
        SurfElem = quadrilateral(surf[0:4],surf[4:])
        self.__quads[id]=SurfElem
        self.__BuildRevCon(id,SurfElem)

# ------ SurfaceMesh --------

    def GiveBounds(self):
        '''
        return two Vec3D obects that describe the bounding box of the surface
        mesh.
        '''

        for node in self.__Nodes.values():
            if node.x() > self.maxX: self.maxX = node.x()
            if node.y() > self.maxY: self.maxY = node.y()
            if node.z() > self.maxZ: self.maxZ = node.z()
            if node.x() < self.minX: self.minX = node.x()
            if node.y() < self.minY: self.minY = node.y()
            if node.z() < self.minZ: self.minZ = node.z()

        return (Vec3D.Vec3D(self.maxX,self.maxY,self.maxZ),\
                Vec3D.Vec3D(self.minX,self.minY,self.minZ))

# ------ SurfaceMesh --------

    def GiveListofNodes(self):
        x = self.__Nodes.keys()
        return x.sort()

# ------ SurfaceMesh --------

    def GiveNodalPosition(self,nid):
        return self.__Nodes[nid]

# ------ SurfaceMesh --------

    def __BuildRevCon(self,seid,surfelem):
        '''
        builds revcon dictionay...
        ''' 

        # do CORNER nodes
        for i in xrange(surfelem.num_corn_nodes):
            elist=self.__revcon.get(surfelem.corners[i],[])
            elist.append(seid)
            self.__revcon[surfelem.corners[i]]=elist

        # do MIDSIDE nodes
        for i in xrange(surfelem.num_mid_nodes):
            elist=self.__revcon.get(surfelem.mid_sides[i],[])
            elist.append(seid)
            self.__revcon[surfelem.mid_sides[i]]=elist

# ------ SurfaceMesh --------

    def Issurf(self,nid):
        '''
        a method that will tell you if a particular node is in the surface
        mesh.
        '''

        is_surf=0 # 0 = is not on surface
                  # 1 = is on surface
        if self.__revcon.has_key(nid):
            is_surf = 1
        else:
            is_surf = 0
                
        return is_surf

# ------ SurfaceMesh --------

    def elgiveme(self,seid):
        '''
        return the surface element connectivity.
        ''' 
        if seid == 'all':
            corn_nodes=[]
            mid_nodes=[]
            for i in xrange(len(self.__elems)):
                corn_nodes+=[self.__elems[i].corners]
                mid_nodes+=[self.__elems[i].mid_sides]
            return corn_nodes,mid_nodes

        else:
            corn_nodes=self.__elems[seid].corners
            mid_nodes=self.__elems[seid].mid_sides
            return corn_nodes,mid_nodes

# ------ SurfaceMesh --------

    def elshowme(self,seid):
        if seid == 'all': 
            for seid in self.__elems.keys():
                corn_nodes=self.__elems[seid].corners
                mid_nodes=self.__elems[seid].mid_sides
                print ' Surface Element -', seid,'-'
                print ' corner nodes: ', corn_nodes
                print ' midside nodes: ', mid_nodes
        else:
            corn_nodes=self.__elems[seid].corners
            mid_nodes=self.__elems[seid].mid_sides
            print ' Surface Element -', seid,'-'
            print ' corner nodes: ', corn_nodes
            print ' midside nodes: ', mid_nodes

# ------ SurfaceMesh --------

    def revgiveme(self,nid):
        if nid == 'all': 
            for i in xrange(len(self.__revcon)):
                node=i
                elems=self.__revcon[i]
                return node,elems
        else:
            node=nid
            elems=self.__revcon[nid]
            return node,elems

# ------ SurfaceMesh --------

    def revshowme(self,nid):
        if nid == 'all': 
            for id in self.__revcon:
                node=id
                elems=self.__revcon[id]
                print ' Surface Node -', node,'-'
                print ' Reverse Connectivity: ', elems
        else:
            node=nid
            elems=self.__revcon[nid]
            print ' Surface Node -', node,'-'
            print ' Reverse Connectivity: ', elems

# ------ SurfaceMesh --------

    # Old code that won't work anymore, but might be useful anyway... 
##    def BuildSrfMsh(self):
##        self.__seid = 0
##        for obj in self.__tris:
##            if self.__tris[obj][0] == 1:
##                self.__elems[self.__seid]=self.__tris[obj][1]
##                self.__AddElem(self.__seid,self.__elems[self.__seid])
##                self.__seid+=1
##            elif self.__tris[obj][0] > 2:
##                print 'not valid mesh; see SurfaceMesh.py'
##            else:
##                continue
##
##        for obj in self.__quads:
##            if self.__quads[obj][0] == 1:
##                self.__elems[self.__seid]=self.__quads[obj][1]
##                self.__AddElem(self.__seid,self.__elems[self.__seid])
##                self.__seid+=1
##            elif self.__quads[obj][0] > 2:
##                print 'not valid mesh; see SurfaceMesh.py'
##            else:
##                continue
##
### ------ SurfaceMesh --------
##
##    def ReNum(self):
##        dumm={}
##        id=0
##        for sel in self.__elems.values():
##            dumm[id]=sel
##            self.__AddElem(id,sel)
##            id+=1
##
##        self.__seid=id
##        self.__elems=dumm

# ------ END of SurfaceMesh --------

class triangle:
    """ this class is a triagular surface element"""

#   triangle numbering
#
#                            
#                         2 o 
#                          / \  
#                         /   \ 
#                        /     \ 
#                       /       \
#                    5 o         o4
#                     /           \ 
#                    /             \ 
#                   /               \ 
#                  o--------o--------o 
#                   0       3        1

## note:  obviously this is bare bones right now and only really stores the
##        connectivity of the triangular surface element    

    def __init__(self,corners,mid_sides):
        self.corners=corners
        self.mid_sides=mid_sides
        self.num_corn_nodes=len(self.corners)
        self.num_mid_nodes=len(self.mid_sides)

    def compare(self,check):
        sum=0 # check each vert in check to see if in self.corners. If yes
              # then add 1 to sum.  At end check if sum = 3 indicating that
              # all three nodes are duplicate and thus the element is a
              # duplicate.
        if check.corners[0] in self.corners:
            sum+=1
        if check.corners[1] in self.corners:
            sum+=1
        if check.corners[2] in self.corners:
            sum+=1
        if sum == 3:
            return 1
        else:
            return 0

class quadrilateral:
    """ this class is a quadrilateral surface element"""

## note:  obviously this is bare bones right now and only really stores the
##        connectivity of the quadrilateral surface element        

    def __init__(self,corners,mid_sides):
        self.corners=corners
        self.mid_sides=mid_sides
        self.num_corn_nodes=len(self.corners)
        self.num_mid_nodes=len(self.mid_sides)

    def compare(self,check):
        sum=0 # check each vert in check to see if in self.corners. If yes
              # then add 1 to sum.  At end check if sum = 3
        if check.corners[0] in self.corners:
            sum+=1
        if check.corners[1] in self.corners:
            sum+=1
        if check.corners[2] in self.corners:
            sum+=1
        if check.corners[4] in self.corners:
            sum+=1
        if sum == 4:
            return 1
        else:
            return 0        