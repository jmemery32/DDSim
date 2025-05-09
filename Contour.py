#	Contour - routine for doing three dimensional contours
#		    on planer polygons.
#
#	This routine uses an algorithm originaly created by TJ Boone.
#
#	Contour(coords,values,cvalues,pgon_func,pgon_data)
#
#	Given the coordinates and contour values for 
#       the verticies of a planer patch, this routine computes the
#       coordinates of the verticies of the descrete contours which
#       pass through the patch and calls the pgon_func to draw
#       the polygons.
#
#	coords    - list of polygon vertex coordinates (x,y,z)
#	values    - values of the contoured parameter at the
#                   polygon verticies
#	cvalues   - values of the "limits" of the contour value
#                   ranges.  If there are n contour colors available
#                   then there should be n-1 of these values.
#                   Vertices with parameter values < cvalues[0]
#                   will be assigned a color index of 0.  Vertices
#                   with parameter values > cvalues[n-2] will be
#                   assigned a color index of n-1.
#	pgon_func - function to call to draw a polygon
#	pgon_data - data to pass to the functions
#
#	the pgon_func routines is called with arguments defined as
#       follows:
#
#       pgon_func(pgon_data,coords,c_index)
#
#	pgon_data -      data passed to Contour
#	coords -         coordinates of the vertices of the
#                        contour polygon
#	c_index -        contour color index

#  Find the index in the cvalues array of a vertex value 

def _GetIndex(value,cvalues):
    if value < cvalues[0]:
        return 0
    elif value >= cvalues[-1]:
        return len(cvalues)
    else:
        for i in xrange(len(cvalues)-1):
            if value >= cvalues[i] and value < cvalues[i+1]:
                return i+1


#  Add one mid point into the lists

def _AddOneMid(index,w1,w2,coord_1,coord_2,con_list):
    con_list[index].append((w1*coord_1[0] + w2*coord_2[0],
                            w1*coord_1[1] + w2*coord_2[1],
                            w1*coord_1[2] + w2*coord_2[2]))


#  Add information for points between vertex points

def _AddMidPts(index_1,index_2,value_1,value_2,coord_1,coord_2,
               cvalues,con_list):
    if index_1 < index_2:
        for i in xrange(index_1,index_2):
            w1 = (value_2-cvalues[i]) / (value_2-value_1)
            w2 = 1.0 - w1
            _AddOneMid(i,w1,w2,coord_1,coord_2,con_list)
            if w2 != 1.0:
                _AddOneMid(i+1,w1,w2,coord_1,coord_2,con_list)
    else:
        for i in xrange(index_1,index_2,-1):
            w1 = (value_2-cvalues[i-1]) / (value_2-value_1)
            w2 = 1.0 - w1
            if w1 != 1.0:
                _AddOneMid(i,w1,w2,coord_1,coord_2,con_list)
            _AddOneMid(i-1,w1,w2,coord_1,coord_2,con_list)


def Contour(coords,values,cvalues,pgon_func,pgon_data):

##    print
##    print 'coords: ',coords
##    print 'values: ',values
##    print 'cvalues: ',cvalues

    # first determine contour level for the corners

    c_index = []
    for value in values:
        c_index.append(_GetIndex(value,cvalues))
##    print 'c_index:', c_index

    # create the num_points list

    con_list = []
    for i in xrange(len(cvalues)+1):
        tmp_list = []
        con_list.append(tmp_list)

    # add the vertex nodes to the lists

    for i in xrange(len(coords)):
        j = (i+1) % len(coords)
        con_list[c_index[i]].append((coords[i][0],coords[i][1],coords[i][2]))

        # if necessary add a mid point
##        print con_list

        if c_index[i] != c_index[j]: 
            _AddMidPts(c_index[i],c_index[j],values[i],values[j],
                       coords[i],coords[j],cvalues,con_list)

    # now display all the necessary polygons
##    print 'con_list:', con_list

    for i in xrange(len(cvalues)+1):
        if len(con_list[i]) >= 3:
            pgon_func(pgon_data,con_list[i],i)


if __name__ == '__main__':
    def pgon_func(file,coords,level):
##        print coords,level
        # writes each surface patch to Sview.svw
        rgb=[[240, 0  , 0  ],
             [255, 120, 120],
             [255, 155, 155],
             [255, 190, 190],
             [255, 250, 250],
             [220, 220, 245],
             [184, 184, 255],
             [148, 148, 255],
             [112, 112, 255],
             [ 0 ,  0 , 240]]

        # loop over polygon coords and write to file
        text='p '+str(len(coords))+' '
        for i in coords:
            text+=str(i[0])+' '+str(i[1])+' '+str(i[2])+' '

        text+=str(rgb[level-1][0]/255.)+' '+str(rgb[level-1][1]/255.)+ \
               ' '+str(rgb[level-1][2]/255.)+' '
        file.write(text+"\n")

    # open the file
    file = open('''test_contour.svw''', 'w')

    # first make cvalues
    cvalues = range(0,1000,100)

    # make up some geometry (a square with two 6-noded surface triangles)
    coords=[[[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0]], \
            [[1, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0]], \
            [[0, 1, 0], [0, 0.5, 0], [0.5, 0.5, 0]], \
            [[0.5, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0]], \
            [[1, 0, 0], [1, 0.5, 0], [0.5, 0.5, 0]], \
            [[1, 1, 0], [0.5, 1, 0], [1, 0.5, 0]], \
            [[0.5, 1, 0], [0, 1, 0], [0.5, 0.5, 0]], \
            [[0.5, 0.5, 0], [1, 0.5, 0], [0.5, 1, 0]]]
    values=[[800,800,800], [800,900,800], [800,800,900], [800,900,800], \
            [800,800,900], [800,800,800], [800,800,900], [900,800,800]]

##    coords=[[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], \
##            [1, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0], \
##            [0, 1, 0], [0, 0.5, 0], [0.5, 0.5, 0], \
##            [0.5, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0]]
##    values=[900,900,900,   900,900,900,     900,900,900, 900,900,900]

##    coords=[[0,0,0],[1,0,0],[0,1,0]]
##    values=[800,700,700]

    for i in xrange(len(values)):
        coord=coords[i]
        value=values[i]
        Contour(coord,value,cvalues,pgon_func,file)
