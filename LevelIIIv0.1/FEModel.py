import Vec3D
import ElementClasses

class FEModel:
    '''
    class that reads Gerd's file format and stores the finite element model.

    I read the following files:
        MyFirstName.con - element connectivity
        MyFirstName.edg - edges of quadratic elements
        MyFirstName.nod - nodes
        MyFirstName.MPC - the multipoint constraint files
        MyFirstName.pr  - pressure boundary conditions
        MyFirstName.Materials - Materials file
        MyFirstName.MaterialModelMap - MaterialModelMap
        MyFirstName.ShapeMap - ShapeMap
    '''

    def __init__(self,MyFirstName):
        