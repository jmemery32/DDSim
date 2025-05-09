
import MeshTools
import Vec3D
import ColTensor

# Exceptions

CompStress = "Craped out in CompareStress()"

CompDisp = "Craped out in CompareDisp()"

def IsPointOutsideTest():
    '''
    meshtools fails for poorly shaped elements.
    ''' 
    path='''h:\\users\\jme32\\for_wash\\'''
    filename='''Tray_Wall'''
    model1 = MeshTools.MeshTools(path+filename,'RDB')
    model1.SetPointInsideTolerance(0.00000000001)

    # qpnt = u + trans
    qpnt=Vec3D.Vec3D(-1.63307890498,-6.22368282107,-17.7654257505)
    print model1.IsPointOutsideMesh(qpnt)

    u = Vec3D.Vec3D(-0.00030890498332,0.0010671789323,0.000474249487647)
    trans = Vec3D.Vec3D(-1.63277,-6.22475,-17.7659)
    print "   0.9*u:", model1.IsPointOutsideMesh(0.9*u+trans)
    print "0.9999*u:", model1.IsPointOutsideMesh(0.9999*u+trans)
    print "   1.0*u:", model1.IsPointOutsideMesh(u+trans)
    print " 1.001*u:", model1.IsPointOutsideMesh(1.001*u+trans)
    print "   1.1*u:", model1.IsPointOutsideMesh(1.1*u+trans)

def CompareStress(stress1,stress2,zero):
    '''
    function to compare the individual componenets of two stress tensors.
    '''

    found = 0
    stuff=[() for i in xrange(6)]
    if abs(stress1.xx()-stress2.xx()) > zero:
        found = 1
        stuff[0] +=('stress xx not the same ',stress1.xx(),stress2.xx())
    if abs(stress1.yy()-stress2.yy()) > zero:
        found = 1
        stuff[1] +=('stress yy not the same',stress1.yy(),stress2.yy())
    if abs(stress1.zz()-stress2.zz()) > zero:
        found = 1
        stuff[2] +=('stress zz not the same',stress1.zz(),stress2.zz())
    if abs(stress1.xy()-stress2.xy()) > zero:
        found = 1
        stuff[3] +=('stress xy not the same',stress1.xy(),stress2.xy())
    if abs(stress1.yz()-stress2.yz()) > zero:
        found = 1
        stuff[4] +=('stress yz not the same',stress1.yz(),stress2.yz())
    if abs(stress1.zx()-stress2.zx()) > zero:
        found = 1
        stuff[5] +=('stress zx not the same',stress1.zx(),stress2.zx())

    if found:
        raise CompStress, stuff
    else:
        return 1

def CompareDisp(disp1,disp2,zero):
    '''
    function to compare the individual componenets of two stress tensors.
    '''

    found = 0
    stuff=[() for i in xrange(6)]
    if abs(disp.x()-disp2.x()) > zero:
        found = 1
        stuff[0] +=('disp x not the same ',disp1.x(),disp2.x())
    if abs(disp.y()-disp2.y()) > zero:
        found = 1
        stuff[1] +=('disp y not the same',disp1.y(),disp2.y())
    if abs(disp.z()-disp2.z()) > zero:
        found = 1
        stuff[2] +=('disp z not the same',disp1.z(),disp2.z())

    if found:
        raise CompDisp, stuff
    else:
        return 1

def CompareToItself(model1,model2,list1,list2,zero):
    '''
    function that loops over doid_list to compare each MeshTools function
    from the two instance to each other to look for potential inconsistencies
    '''

    #   1
    print "  TEST locations, stress and disp from GetNodeInfo:"
    for nid1 in list1:
        nid2=list2[list1.index(nid1)]
        if nid1 != nid2:
            print 'lists are not identical'
        (loc1,disp1,stress1) = model1.GetNodeInfo(nid1)
        (loc2,disp2,stress2) = model2.GetNodeInfo(nid2)
        if (loc1-loc2).Magnitude() > zero:
            raise 'locations are not the same'
        try: 
            samestress = CompareStress(stress1,stress2,zero)
        except KeyboardInterrupt:
                raise 'KeyboardInterrupt'
        except CompStress, msg:
            errormsg=''
            for ms in msg:
                errormsg += ms[0] + ' at node %s and %s' %(nid1,nid2)+"\n"
                errormsg += '   the stress are %1.6e and %1.6e' % (ms[1],ms[2])
                errormsg += "\n"
            print errormsg

        # compare displacements: 
        try: 
            samedisp = CompareDisp(disp1,disp2,zero)
        except KeyboardInterrupt:
                raise 'KeyboardInterrupt'
        except CompDisp, msg:
            errormsg=''
            for ms in msg:
                errormsg += ms[0] + ' at node %s and %s' %(nid1,nid2)+"\n"
                errormsg += '    the disp are %1.6e and %1.6e' % (ms[1],ms[2])
                errormsg += "\n"
            print errormsg
    print "\n"

    #   2
    print "  TEST GetPtStress vs. GetNodeInfo (model1):"
    for nid1 in list1:
        (loc1,disp1,stress1) = model1.GetNodeInfo(nid1)
        stress11 = model1.GetPtStress(loc1)
        try:
            samestress = CompareStress(stress11,stress1,zero)
        except KeyboardInterrupt:
                raise 'KeyboardInterrupt'
        except CompStress, msg:
            errormsg=''
            for ms in msg:
                errormsg += ms[0] + ' at node %s in model 1' %(nid1)+"\n"
                errormsg += '    the stress are %1.6e and %1.6e' % (ms[1],ms[2])
                errormsg += "\n"
            print errormsg
    print "\n"

    #   3
    print "  TEST GetPtStress vs. GetNodeInfo (model2):"
    for nid2 in list2:
        (loc2,disp2,stress2) = model2.GetNodeInfo(nid2)
        stress22 = model2.GetPtStress(loc2)
        try:
            samestress = CompareStress(stress22,stress2,zero)
        except KeyboardInterrupt:
                raise 'KeyboardInterrupt'
        except CompStress, msg:
            errormsg=''
            for ms in msg:
                errormsg += ms[0] + ' at node %s in model 2' %(nid2)+"\n"
                errormsg += '    the stress are %1.6e and %1.6e' % (ms[1],ms[2])
                errormsg += "\n"
            print errormsg
    print "\n"

    #   4
    print "  TEST GetPtStress at a Vec3D"
    for nid1 in list1:
        nid2=list2[list1.index(nid1)]
        if nid1 != nid2:
            print 'lists are not identical'
        (loc1,disp1,stress1) = model1.GetNodeInfo(nid1)
        stress11 = model1.GetPtStress(loc1)
        (loc2,disp2,stress2) = model2.GetNodeInfo(nid2)
        stress22 = model2.GetPtStress(loc2)
        try:
            samestress = CompareStress(stress11,stress22,zero)
        except KeyboardInterrupt:
                raise 'KeyboardInterrupt'
        except CompStress, msg:
            errormsg=''
            for ms in msg:
                errormsg += ms[0] + ' at node %s and %s' %(nid1,nid2)+"\n"
                errormsg += '    the stress are %1.6e and %1.6e' % (ms[1],ms[2])
                errormsg += "\n"
            print errormsg
    print "\n"

        #currently not tested:  
##        .GetPtDisp(pnt)
##        .GetElemInfo(eid)
##        .GetAdjacentElems(nid)
##        .GetSurfElemInfo(eid)
##        .GetAdjacentSurfElems(nid)
##        .IsPointOutsideMesh(pnt1)
##        .IsSurfaceNode(nid)
##        .SurfaceNormal(nid,5)
##        .GetAdjacentSurfEdgeLengths(nid)

def Func(model,qpnt):
    '''
    a play funciton to test how MeshTools performs if called from within a
    function.
    '''

    print ' Pnt stress:'
    print 'MeshTools',
    sig = model.GetPtStress(qpnt)
    sigG=[[sig.xx(), sig.xy(), sig.zx()],
          [sig.xy(), sig.yy(), sig.yz()],
          [sig.zx(), sig.yz(), sig.zz()]]
    print sig

def MakeLists(filename,conpath,surface,doid_range,model):
    '''
    a function to make lists of nodes from the model stored as a MeshTools
    object.
    '''

    xyz_data  = open(conpath+filename+'''.nod''', 'r')
    xyz=xyz_data.readlines()
    doid_list =[]
    node_list = []

    for i in range(len(xyz)):
        line=xyz[i].split()
        nid=int(line[0])
        x_val=float(line[1])
        y_val=float(line[2])
        z_val=float(line[3])
        node_list += [nid]
        if doid_range: # pick damage origins from this domain only
            if x_val >= doid_range[0][0] and x_val <=doid_range[1][0] and \
               y_val >= doid_range[0][1] and y_val <=doid_range[1][1] and \
               z_val >= doid_range[0][2] and z_val <=doid_range[1][2]:
               doid_list+=[nid] # damage origin id
        else:
            doid_list+=[nid]

    temp_list=[]

    if surface == 1:
        for i in doid_list:
            issurf = model.IsSurfaceNode(i)
            if issurf == 1:
                temp_list+=[i]
        doid_list=temp_list

    doid_list.sort()
    node_list.sort()

    xyz_data.close()

    return doid_list, node_list

def OrgTest(mmodel,pnt,eid,nid):
    '''
    Original test of MeshTools.  Compares to answers we got from FemModel for
    the simple cube model.

    pnt - a Vec3D object representing a point in space
    eid - an integer element i.d.
    nid - an integer node i.d.
    '''

    print ' Pnt stress:'
    print 'MeshTools',mmodel.GetPtStress(pnt)
    print 'FemModel (8.08523e-013 3.4 2.62242e-012 4.326e-013 1.00179e-014 4.7595e-012)'
    print ' '
    print ' Pnt Disp:'
    print 'MeshTools', mmodel.GetPtDisp(pnt)
    print 'FemModel None'
    print ' '
    print ' Node Info'
    print 'MeshTools', mmodel.GetNodeInfo(nid)
    print 'FemModel  (10 10 0, None, (-1.48313e-013 3.4 6.53739e-013 5.40886e-013 -4.63987e-013 1.45791e-012))'
    print ' '
    print ' Element Info' 
    print 'MeshTools', mmodel.GetElemInfo(eid)
    print 'FemModel [0, 6, 3, 8, 14, 24, 11, 25, 33, 15]'
    print ' '
    print ' Adjacent Elements'
    print 'MeshTools', mmodel.GetAdjacentElems(nid)
    print 'FemModel [12, 5]'
    print ' '
    print ' Surface Element info'
    print 'MeshTools', mmodel.GetSurfElemInfo(eid)
    print 'FemModel ([0, 3, 6], [11, 24, 14])'
    print ' '
    print ' Adjacent surface elements'
    print 'MeshTools', mmodel.GetAdjacentSurfElems(nid)
    print 'FemModel [2, 5]'
    print ' '
    print ' Is point outside mesh?' 
    print 'MeshTools', mmodel.IsPointOutsideMesh(pnt)
    print 'FemModel 1'
    print ' '
    print ' Is node on the surface?' 
    print 'MeshTools', mmodel.IsSurfaceNode(nid)
    print 'FemModel 1'
    print ' '
    print ' Surface Normal' 
    print 'MeshTools', mmodel.SurfaceNormal(nid)
    print 'FemModel (0.707107 0.707107 0, 0.35355339059327379)'
    print ' '
    print ' Adjacent surface edge lengths' 
    print 'MeshTools', mmodel.GetAdjacentSurfEdgeLengths(nid)
    print 'FemModel [10.0, 10.0]'
    print ' '

def SelfConsistent(model,list,zero):
    '''
    function to query at a location a couple times and check if MeshTools
    returns the same stress.
    '''

    print " TEST stress from GetNodeInfo:"
    for nid in list:
        (loc1,disp1,stress1) = model.GetNodeInfo(nid)
        (loc2,disp2,stress2) = model.GetNodeInfo(nid)

        try:
            samestress = CompareStress(stress1,stress2,zero)
        except KeyboardInterrupt:
                raise 'KeyboardInterrupt'
        except CompStress, msg:
            errormsg=''
            for ms in msg:
                errormsg += ms[0] + ' at node %s' %(nid)+"\n"
                errormsg += '    the stress are %1.6e and %1.6e' % (ms[1],ms[2])
                errormsg += "\n"
            print errormsg
    print "\n"

    print "  TEST GetPtStress:"
    for nid in list:
        (loc1,disp1,stress1) = model.GetNodeInfo(nid)
        (loc2,disp2,stress2) = model.GetNodeInfo(nid)
        stress11 = model.GetPtStress(loc1)
        stress22 = model.GetPtStress(loc2)
        try:
            samestress = CompareStress(stress11,stress22,zero)
        except KeyboardInterrupt:
                raise 'KeyboardInterrupt'
        except CompStress, msg:
            errormsg=''
            for ms in msg:
                errormsg += ms[0] + ' at node %s' %(nid)+"\n"
                errormsg += '    the stress are %1.6e and %1.6e' % (ms[1],ms[2])
                errormsg += "\n"
            print errormsg
    print "\n"
        

if __name__ == "__main__":
    '''
    Test MeshTools.pyd functionality.  Functions include:

    .GetPtStress(pnt)
    .GetPtDisp(pnt)
    .GetNodeInfo(nid)
    .GetElemInfo(eid)
    .GetAdjacentElems(nid)
    .GetSurfElemInfo(eid)
    .GetAdjacentSurfElems(nid)
    .IsPointOutsideMesh(pnt1)
    .IsSurfaceNode(nid)
    .SurfaceNormal(nid,5)
    .GetAdjacentSurfEdgeLengths(nid)
    .GetMaxDimension() --> max,min (Vec3D's)
    '''


    # Different models....
##    path='''c:\\john\\research\\SIPS\\franc3d\\coupon01\\'''
##    filename='''coupon01'''
##    path='''examples\\'''
##    filename='''example1'''
##    path='''Z:\\research\\URETI\\code\\python\\DDSimv1.4\\examples\\'''
##    path='''z:\\research\\SIPS\\franc3d\\coupon02\\c2\\Determ_016\\'''
##    filename='''coupon02X3'''
##    path='''h:\\users\\jme32\\research\\URETI\\franc3d\\Tray_wall\\'''
##    path='''h:\\users\\jme32\\for_wash\\'''
##    filename='''Tray_Wall'''
##    path='''z:\\Gerd_FEM\\cube_20x20_bend\\'''
##    filename="cube_20x20_bend"
##    path='''Z:\\Gerd_FEM\\WedgeElements\\LinearWedges\\TwoWedges\\'''
##    filename='''TwoWedgies'''
    path='''Z:\\research\\sips\\franc3d\\SIPS3002\\'''
    filename='''SIPS3002'''

    model1 = MeshTools.MeshTools(path+filename,'RDB')
    model1.SetPointInsideTolerance(0.00000000001)

    print model1.GetMaxDimension()
##    model2 = MeshTools.MeshTools(path+filename,'FEM')
##
##    node_list1,doid_list1=MakeLists(filename,path,1,None,model1)
##    node_list2,doid_list2=MakeLists(filename,path,1,None,model2)
##
##    zero=1.0e-3
##
##    print "CompareToItself"
##    # function that loops over doid_list to compare each MeshTools function
##    # from the two instance to each other to look for potential inconsistencies
##    CompareToItself(model1,model2,doid_list1,doid_list2,zero)
##
##    print "SeflConsistent"
##    # function to compare multiple stress queries at the same point to look for
##    # inconsistencies
##    SelfConsistent(model1,doid_list1,zero)

##    for i in range(1000):
##        OrgTest(model1,Vec3D.Vec3D(1,1,1),5,0)
    











