'''
Usage: DDSimFiles [options] file

Options:
        -v              : verbose mode
        -t              : type
        --binary        : create binary files
'''

import os, struct

verbose = 0
binary = 0

types = {
    "con" : "Elements",
    "edg" : "Edges",
    "nod" : "Nodes",
    "sig" : "Stress",
    "smp" : "ShapeMap"
}

surface_triangle = {}
surface_quadrilateral = {}

################################################################################

def usage(status, msg=''):
    "Error message and usage"
    print __doc__
    if msg:
        print '-- ERROR: %s' % msg
    sys.exit(status)

################################################################################

def write_nodes(fname):
    fin = open(fname, 'r')
    if binary: fmt = 'wb'
    else: fmt = 'w'
    fnod = open(os.path.dirname(fname)+'\Nodes.csv', fmt)
    buff = fin.readline().strip().split()
    while len(buff):
        if len(buff) != 4:
            usage(1, 'Invalid entry in nodes file found')
            return None
        NodeID = int(buff[0])
        x = float(buff[1])
        y = float(buff[2])
        z = float(buff[3])
        if binary: fnod.write(struct.pack('=i3d', NodeID, x, y, z))
        else: fnod.write('%d,%.15e,%.15e,%.15e\n' % (NodeID, x, y, z))

        buff = fin.readline().strip().split()

    fnod.close()
    fin.close()

################################################################################

def add_brick_surface_quadrilaterals(ElemID, v0, v1, v2, v3, v4, v5, v6, v7):
    global surface_quadrilateral
    key = []
    value = []

    value.append(ElemID + ':' + '0:' + v0 + ':' + v1 + ':' + v3 + ":" + v2)
    t = [v0, v1, v3, v2]
    t.sort()
    key.append(t[0] + ':' + t[1] + ':' + t[2] + ':' + t[3])

    value.append(ElemID + ':' + '1:' + v4 + ':' + v6 + ':' + v7 + ":" + v5)
    t = [v4, v6, v7, v5]
    t.sort()
    key.append(t[0] + ':' + t[1] + ':' + t[2] + ':' + t[3])

    value.append(ElemID + ':' + '2:' + v0 + ':' + v4 + ':' + v5 + ":" + v1)
    t = [v0, v4, v5, v1]
    t.sort()
    key.append(t[0] + ':' + t[1] + ':' + t[2] + ':' + t[3])

    value.append(ElemID + ':' + '3:' + v2 + ':' + v3 + ':' + v7 + ":" + v6)
    t = [v2, v3, v7, v6]
    t.sort()
    key.append(t[0] + ':' + t[1] + ':' + t[2] + ':' + t[3])

    value.append(ElemID + ':' + '4:' + v0 + ':' + v2 + ':' + v6 + ":" + v4)
    t = [v0, v2, v6, v4]
    t.sort()
    key.append(t[0] + ':' + t[1] + ':' + t[2] + ':' + t[3])

    value.append(ElemID + ':' + '5:' + v1 + ':' + v5 + ':' + v7 + ":" + v3)
    t = [v1, v5, v7, v3]
    t.sort()
    key.append(t[0] + ':' + t[1] + ':' + t[2] + ':' + t[3])


    for i in range(len(key)):
        if surface_quadrilateral.has_key(key[i]):
            del surface_quadrilateral[key[i]]
        else:
            surface_quadrilateral[key[i]] = value[i]

################################################################################

def add_pyramid_faces(ElemID, v0, v1, v2, v3, v4):
    global surface_triangle
    global surface_quadrilateral
    tri_key = []
    tri_value = []

    quad_key = []
    quad_value = []

    quad_value.append(ElemID + ':' + '0:' + v0 + ':' + v1 + ':' + v2 + ":" + v3)
    t = [v0, v1, v2, v3]
    t.sort()
    quad_key.append(t[0] + ':' + t[1] + ':' + t[2] + ':' + t[3])

    tri_value.append(ElemID + ':' + '1:' + v0 + ':' + v4 + ':' + v1)
    t = [v0, v4, v1]
    t.sort()
    tri_key.append(t[0] + ':' + t[1] + ':' + t[2])

    tri_value.append(ElemID + ':' + '2:' + v1 + ':' + v4 + ':' + v2)
    t = [v1, v4, v2]
    t.sort()
    tri_key.append(t[0] + ':' + t[1] + ':' + t[2])

    tri_value.append(ElemID + ':' + '3:' + v2 + ':' + v4 + ':' + v3)
    t = [v2, v4, v3]
    t.sort()
    tri_key.append(t[0] + ':' + t[1] + ':' + t[2])

    tri_value.append(ElemID + ':' + '4:' + v3 + ':' + v4 + ':' + v0)
    t = [v3, v4, v0]
    t.sort()
    tri_key.append(t[0] + ':' + t[1] + ':' + t[2])

    for i in range(len(quad_key)):
        if surface_quadrilateral.has_key(quad_key[i]):
            del surface_quadrilateral[quad_key[i]]
        else:
            surface_quadrilateral[quad_key[i]] = quad_value[i]

    for i in range(len(tri_key)):
        if surface_triangle.has_key(tri_key[i]):
            del surface_triangle[tri_key[i]]
        else:
            surface_triangle[tri_key[i]] = tri_value[i]

################################################################################

def add_tet_surface_triangles(ElemID, v0, v1, v2, v3):
    global surface_triangle
    key = []
    value = []

    value.append(ElemID + ':' + '0:' + v0 + ':' + v1 + ':' + v2)
    t = [v0, v1, v2]
    t.sort()
    key.append(t[0] + ':' + t[1] + ':' + t[2])

    value.append(ElemID + ':' + '1:' + v1 + ':' + v3 + ':' + v2)
    t = [v1, v3, v2]
    t.sort()
    key.append(t[0] + ':' + t[1] + ':' + t[2])

    value.append(ElemID + ':' + '2:' + v2 + ':' + v3 + ':' + v0)
    t = [v2, v3, v0]
    t.sort()
    key.append(t[0] + ':' + t[1] + ':' + t[2])

    value.append(ElemID + ':' + '3:' + v0 + ':' + v3 + ':' + v1)
    t = [v0, v3, v1]
    t.sort()
    key.append(t[0] + ':' + t[1] + ':' + t[2])

    for i in range(len(key)):
        if surface_triangle.has_key(key[i]):
            del surface_triangle[key[i]]
        else:
            surface_triangle[key[i]] = value[i]

################################################################################

def add_wedge_faces(ElemID, v0, v1, v2, v3, v4, v5):
    global surface_triangle
    global surface_quadrilateral
    tri_key = []
    tri_value = []

    quad_key = []
    quad_value = []

    quad_value.append(ElemID + ':' + '0:' + v0 + ':' + v3 + ':' + v4 + ":" + v1)
    t = [v0, v3, v4, v1]
    t.sort()
    quad_key.append(t[0] + ':' + t[1] + ':' + t[2] + ':' + t[3])

    quad_value.append(ElemID + ':' + '1:' + v2 + ':' + v5 + ':' + v3 + ":" + v0)
    t = [v2, v5, v3, v0]
    t.sort()
    quad_key.append(t[0] + ':' + t[1] + ':' + t[2] + ':' + t[3])

    quad_value.append(ElemID + ':' + '2:' + v1 + ':' + v4 + ':' + v5 + ":" + v2)
    t = [v1, v4, v5, v2]
    t.sort()
    quad_key.append(t[0] + ':' + t[1] + ':' + t[2] + ':' + t[3])

    tri_value.append(ElemID + ':' + '3:' + v0 + ':' + v1 + ':' + v2)
    t = [v0, v1, v2]
    t.sort()
    tri_key.append(t[0] + ':' + t[1] + ':' + t[2])

    tri_value.append(ElemID + ':' + '4:' + v3 + ':' + v5 + ':' + v4)
    t = [v3, v5, v4]
    t.sort()
    tri_key.append(t[0] + ':' + t[1] + ':' + t[2])

    for i in range(len(quad_key)):
        if surface_quadrilateral.has_key(quad_key[i]):
            del surface_quadrilateral[quad_key[i]]
        else:
            surface_quadrilateral[quad_key[i]] = quad_value[i]

    for i in range(len(tri_key)):
        if surface_triangle.has_key(tri_key[i]):
            del surface_triangle[tri_key[i]]
        else:
            surface_triangle[tri_key[i]] = tri_value[i]

################################################################################

def write_elements(fname):

    fin = open(fname, 'r')
    if binary: fmt = 'wb'
    else: fmt = 'w'

    fhex = open(os.path.dirname(fname)+'\Bricks.csv', fmt)
    fhexvert = open(os.path.dirname(fname)+'\BrickVertices.csv', fmt)

    fpyr = open(os.path.dirname(fname)+'\Pyramids.csv', fmt)
    fpyrvert = open(os.path.dirname(fname)+'\PyramidVertices.csv', fmt)

    ftet = open(os.path.dirname(fname)+'\Tetrahedra.csv', fmt)
    ftetvert = open(os.path.dirname(fname)+'\TetrahedronVertices.csv', fmt)

    fwed = open(os.path.dirname(fname)+'\Wedges.csv', fmt)
    fwedvert = open(os.path.dirname(fname)+'\WedgeVertices.csv', fmt)

    buff = fin.readline().strip().split()
    while len(buff):
        ElemID = buff[0]
        MatID = buff[2]
        num_vert = int(buff[3])
        if num_vert+4 != len(buff):
            usage(1, 'Invalid entry in elements file found')
            return None

        if num_vert == 4: # tet
            if binary:
                ftet.write(struct.pack('=iBi', int(ElemID), 4, int(MatID)))
            else:
                ftet.write('%s,4,%s\n' % (ElemID, MatID))

            v0 = buff[4]
            v1 = buff[5]
            v2 = buff[6]
            v3 = buff[7]

            if binary:
                ftetvert.write(struct.pack('=iBi', int(ElemID), 0, int(v0)))
                ftetvert.write(struct.pack('=iBi', int(ElemID), 1, int(v1)))
                ftetvert.write(struct.pack('=iBi', int(ElemID), 2, int(v2)))
                ftetvert.write(struct.pack('=iBi', int(ElemID), 3, int(v3)))
            else:
                ftetvert.write('%s,%s,%s\n' % (ElemID, 0, v0))
                ftetvert.write('%s,%s,%s\n' % (ElemID, 1, v1))
                ftetvert.write('%s,%s,%s\n' % (ElemID, 2, v2))
                ftetvert.write('%s,%s,%s\n' % (ElemID, 3, v3))

            add_tet_surface_triangles(ElemID, v0, v1, v2, v3)
            
        elif num_vert == 5: # pyramid 
            if binary:
                fpyr.write(struct.pack('=iBi', int(ElemID), 5, int(MatID)))
            else:
                fpyr.write('%s,5,%s\n' % (ElemID, MatID))

            v0 = buff[4]
            v1 = buff[5]
            v2 = buff[6]
            v3 = buff[7]
            v4 = buff[8]

            if binary:
                fpyrvert.write(struct.pack('=iBi', int(ElemID), 0, int(v0)))
                fpyrvert.write(struct.pack('=iBi', int(ElemID), 1, int(v1)))
                fpyrvert.write(struct.pack('=iBi', int(ElemID), 2, int(v2)))
                fpyrvert.write(struct.pack('=iBi', int(ElemID), 3, int(v3)))
                fpyrvert.write(struct.pack('=iBi', int(ElemID), 4, int(v4)))
            else:
                fpyrvert.write('%s,%s,%s\n' % (ElemID, 0, v0))
                fpyrvert.write('%s,%s,%s\n' % (ElemID, 1, v1))
                fpyrvert.write('%s,%s,%s\n' % (ElemID, 2, v2))
                fpyrvert.write('%s,%s,%s\n' % (ElemID, 3, v3))
                fpyrvert.write('%s,%s,%s\n' % (ElemID, 4, v4))

            add_pyramid_faces(ElemID, v0, v1, v2, v3, v4)

        elif num_vert == 6: # wedge 
            if binary:
                fwed.write(struct.pack('=iBi', int(ElemID), 6, int(MatID)))
            else:
                fwed.write('%s,6,%s\n' % (ElemID, MatID))

            v0 = buff[4]
            v1 = buff[5]
            v2 = buff[6]
            v3 = buff[7]
            v4 = buff[8]
            v5 = buff[9]

            if binary:
                fwedvert.write(struct.pack('=iBi', int(ElemID), 0, int(v0)))
                fwedvert.write(struct.pack('=iBi', int(ElemID), 1, int(v1)))
                fwedvert.write(struct.pack('=iBi', int(ElemID), 2, int(v2)))
                fwedvert.write(struct.pack('=iBi', int(ElemID), 3, int(v3)))
                fwedvert.write(struct.pack('=iBi', int(ElemID), 4, int(v4)))
                fwedvert.write(struct.pack('=iBi', int(ElemID), 5, int(v5)))
            else:
                fwedvert.write('%s,%s,%s\n' % (ElemID, 0, v0))
                fwedvert.write('%s,%s,%s\n' % (ElemID, 1, v1))
                fwedvert.write('%s,%s,%s\n' % (ElemID, 2, v2))
                fwedvert.write('%s,%s,%s\n' % (ElemID, 3, v3))
                fwedvert.write('%s,%s,%s\n' % (ElemID, 4, v4))
                fwedvert.write('%s,%s,%s\n' % (ElemID, 5, v5))

            add_wedge_faces(ElemID, v0, v1, v2, v3, v4, v5)

        elif num_vert == 8: # brick 
            if binary:
                fhex.write(struct.pack('=iBi', int(ElemID), 8, int(MatID)))
            else:
                fhex.write('%s,8,%s\n' % (ElemID, MatID))

            v0 = buff[4]
            v1 = buff[5]
            v2 = buff[6]
            v3 = buff[7]
            v4 = buff[8]
            v5 = buff[9]
            v6 = buff[10]
            v7 = buff[11]

            if binary:
                fhexvert.write(struct.pack('=iBi', int(ElemID), 0, int(v0)))
                fhexvert.write(struct.pack('=iBi', int(ElemID), 1, int(v1)))
                fhexvert.write(struct.pack('=iBi', int(ElemID), 2, int(v2)))
                fhexvert.write(struct.pack('=iBi', int(ElemID), 3, int(v3)))
                fhexvert.write(struct.pack('=iBi', int(ElemID), 4, int(v4)))
                fhexvert.write(struct.pack('=iBi', int(ElemID), 5, int(v5)))
                fhexvert.write(struct.pack('=iBi', int(ElemID), 6, int(v6)))
                fhexvert.write(struct.pack('=iBi', int(ElemID), 7, int(v7)))
            else:
                fhexvert.write('%s,%s,%s\n' % (ElemID, 0, v0))
                fhexvert.write('%s,%s,%s\n' % (ElemID, 1, v1))
                fhexvert.write('%s,%s,%s\n' % (ElemID, 2, v2))
                fhexvert.write('%s,%s,%s\n' % (ElemID, 3, v3))
                fhexvert.write('%s,%s,%s\n' % (ElemID, 4, v4))
                fhexvert.write('%s,%s,%s\n' % (ElemID, 5, v5))
                fhexvert.write('%s,%s,%s\n' % (ElemID, 6, v6))
                fhexvert.write('%s,%s,%s\n' % (ElemID, 7, v7))

            add_brick_surface_quadrilaterals(ElemID, v0, v1, v2, v3, v4, v5, v6, v7)

        buff = fin.readline().strip().split()

    fwedvert.close()
    fwed.close()
    ftetvert.close()
    ftet.close()
    fpyrvert.close()
    fpyr.close()
    fhexvert.close()
    fhex.close()

    fin.close()

################################################################################

def write_edges(fname):

    fin = open(fname, 'r')
    if binary: fmt = 'wb'
    else: fmt = 'w'
    fedg = open(os.path.dirname(fname)+'\Edges.csv', fmt)
    buff = fin.readline().strip().split()
    while len(buff):
        if len(buff) != 5:
            usage(1, 'Invalid entry in edges file found')
            return None
        EdgeID = buff[0]
        v0 = buff[3]
        v1 = buff[4]
        if binary: fedg.write(struct.pack('=3i', int(v0), int(v1), int(EdgeID)))
        else: fedg.write('%s,%s,%s\n' % (v0, v1, EdgeID))

        buff = fin.readline().strip().split()

    fedg.close()
    fin.close()

################################################################################

def write_displacement(fname):

    fin = open(fname, 'r')
    if binary: fmt = 'wb'
    else: fmt = 'w'
    fdsp = open(os.path.dirname(fname)+'\Displacement.csv', fmt)
    buff = fin.readline().strip().split()
    while len(buff):
        if len(buff) != 4:
            usage(1, 'Invalid entry in displacement file found')
            return None
        NodeID = buff[0]
        dx = float(buff[1])
        dy = float(buff[2])
        dz = float(buff[3])
        if binary: fdsp.write(struct.pack('=i3d', NodeID, dx, dy, dz))
        else:
            fdsp.write('%s,%.15e,%.15e,%.15e\n' % (NodeID, dx, dy, dz))

        buff = fin.readline().strip().split()

    fdsp.close()
    fin.close()

################################################################################

def write_nodal_stress(fname):

    fin = open(fname, 'r')
    if binary: fmt = 'wb'
    else: fmt = 'w'
    fsig = open(os.path.dirname(fname)+'\NodalStress.csv', fmt)
    buff = fin.readline().strip().split()
    while len(buff):
        if len(buff) != 7:
            usage(1, 'Invalid entry in nodal_stress file found')
            return None
        NodeID = int(buff[0])
        sxx = float(buff[1])
        syy = float(buff[2])
        szz = float(buff[3])
        sxy = float(buff[4])
        syz = float(buff[5])
        szx = float(buff[6])
        if binary:
            fsig.write(struct.pack('=i6d', NodeID, sxx, syy, szz, sxy, syz, szx))
        else:
            fsig.write('%d,%.15e,%.15e,%.15e,%.15e,%15e,%.15e\n' % \
                (NodeID, sxx, syy, szz, sxy, syz, szx))

        buff = fin.readline().strip().split()

    fsig.close()
    fin.close()

################################################################################

def write_surface_triangles(fname):

    if binary: fmt = 'wb'
    else: fmt = 'w'
    fp = open(os.path.dirname(fname)+'\SurfaceTriangles.csv', fmt)
    fp1 = open(os.path.dirname(fname)+'\SurfaceTriangleVertices.csv', fmt)

    global surface_triangle
    count = 0
    for k in surface_triangle.keys():
        s = surface_triangle[k].split(':')
        if binary:
            fp.write(struct.pack('=iBiB', count, 3, int(s[0]), int(s[1])))
            fp1.write(struct.pack('=iBi', count, 0, int(s[2])))
            fp1.write(struct.pack('=iBi', count, 1, int(s[3])))
            fp1.write(struct.pack('=iBi', count, 2, int(s[4])))
        else:
            fp.write('%d,%d,%s,%s\n' % (count, 3, s[0], s[1]))
            fp1.write('%d,%s,%s,%s\n' % (count, s[2], s[3], s[4]))
            #fp1.write('%d,1,%s\n' % (count, s[3]))
            #fp1.write('%d,2,%s\n' % (count, s[4]))

        count += 1

    fp1.close()
    fp.close()

################################################################################

def write_surface_quadrilaterals(fname):

    if binary: fmt = 'wb'
    else: fmt = 'w'
    fp = open(os.path.dirname(fname)+'\SurfaceQuadrilaterals.csv', fmt)
    fp1 = open(os.path.dirname(fname)+'\SurfaceQuadrilateralVertices.csv', fmt)

    global surface_quadrilateral
    count = len(surface_triangle) 
    for k in surface_quadrilateral.keys():
        s = surface_quadrilateral[k].split(':')
        if binary:
            fp.write(struct.pack('=iBiB', count, 4, int(s[0]), int(s[1])))
            fp1.write(struct.pack('=iBi', count, 0, int(s[2])))
            fp1.write(struct.pack('=iBi', count, 1, int(s[3])))
            fp1.write(struct.pack('=iBi', count, 2, int(s[4])))
            fp1.write(struct.pack('=iBi', count, 3, int(s[5])))
        else:
            fp.write('%d,%d,%s,%s\n' % (count, 4, s[0], s[1]))
            fp1.write('%d,%s,%s,%s,%s\n' % (count, s[2], s[3], s[4], s[5]))
#            fp1.write('%d,1,%s\n' % (count, s[3]))
#            fp1.write('%d,2,%s\n' % (count, s[4]))
#            fp1.write('%d,3,%s\n' % (count, s[5]))

        count += 1

    fp1.close()
    fp.close()

################################################################################

def main(argv):
    import getopt, os

    type = ''

    try:
        optlist, args = getopt.getopt(argv[1:], 't:hv?', ['binary'])
        if not optlist: usage(0)
    except getopt.error, e:
        usage(1, e)

    for o, a in optlist:
        if   o == '-v' :
            global verbose
            verbose = 1
        elif o == '-t' : type = a
        elif o == '--binary' :
            global binary
            binary = 1
        else: usage(0)

    if not types.has_key(type): usage(1, 'Invalid type \"' + type + '\" specified')

    if not len(args):
        usage(1, 'No filename specified!')
    fname = args[0]

    try:
        file = open(fname, 'r')
        file.close()
    except IOError, e:
        usage(1, e)

    if type == 'con':
        write_elements(fname)
        write_surface_triangles(fname)
        write_surface_quadrilaterals(fname)
    elif type == 'nod': write_nodes(fname)
    elif type == 'edg': write_edges(fname)
    elif type == 'sig': write_nodal_stress(fname)
    elif type == 'dsp': write_displacement(fname)

################################################################################

if __name__ == '__main__':
    import sys
    main(sys.argv)
