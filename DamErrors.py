
# exceptions
class FitPolyError(Exception):
    '''
    for when the fit poly subroutine dies
    '''
    pass

class BuildPhiError(Exception):
    '''
    when __BuildPhi routine returns 1 or 3
    '''
    def __repr__(self):
        return '__BuildPhi routine returns 1 or 3'

class ListsIndexError(Exception):
    '''
    when Damage.Lists() has problems
    '''
    pass

class dAdNError(Exception):
    '''
    raised in __.GrowDam() when zero increment is returned from dadN.pyd
    '''
    pass 