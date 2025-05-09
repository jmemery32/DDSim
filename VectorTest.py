from JohnsVectorTools import *

if __name__ == "__main__":

    a=[1.0, 2.0, 3.0, 5.0]
    b=[-2, 3, 4, -2.0]
    e=0.5

    for i in xrange(10000000):

        c= Plus(a,b)
        d= Minus(a,b)
        f= Star(a,b)
        a= Divide(a,b)
        aa= ScalarMult(a,e)
        i= Divide(a,Star(a,b)) # This one eats memory
        j= Plus(a,Plus(a,b)) # So does this one
        k= Divide(Minus(a,b),Minus(a,b))
        l= Plus(a,Star(a,Minus(a,b)))