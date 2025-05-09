
class DamHistory:
    '''
    a container class to store data at a crack tip.  this could be a lot more
    robust, like could store the plastic zone size and pass from one damage
    element to the next... but i don't currently... 
    '''

    def __init__(self,c,K=None):
        self.c=[c] # length of crack at 
        if K: self.K=[K] # SIF corresponding to a
        else: self.K=[]

    def __getitem__(self,key):
        if key=='c' or key== 0: return self.c
        elif key=='K' or key ==1: return self.K

    def __setitem__(self,key,value):
        if key=='c' or key== 0: self.c=value
        elif key=='K' or key ==1: self.K=value

    def __repr__(self):
        s=" crack lengths (%i): \n" %(len(self.c))
        for i in range(len(self.c)):
            s+=" %1.6e, " %(self.c[i])
        s=s[0:-2]
        s+="\n \n Stress Intenstiy Factors (%i): \n" %(len(self.K))
        for i in range(len(self.K)):
            s+=" %1.6e, " %(self.K[i])
        s=s[0:-2]
        return s+"\n"

    def __len__(self):
        return len(self.c)

if __name__=="__main__":

    a=DamHistory(2.3e-3,.134)
    for i in range(10):
        a.c+=[a.c[i]+0.01]
        a.K+=[a.K[i]+0.1]
##        a[1]+=[a.K[i]+0.1]

    print a
##    print a['c'][0]
##    print a[0][0]
##    print a[1][0]


