import math


class OneVariable:
    '''
    store points of a function of one variable, y = f(x), and interpolate
    for function values at a point (a) not stored.
    '''

    def __init__(self,x,y=None):
        '''
        x = list of coorinate points
        y = list of function values.  Should match x so that:
                    f(x[i]) = y[i]
            (can be added one at a time with .AddFuncVal(i,y) or .Insert(x,y))
        '''
        self.x=x
        self.Data={} # use corresponding x index as key!
        if y:
            for i in range(len(x)):
                self.AddFuncVal(i,y[i])

        self.Ordered = False
        self.Interp='LOG' # set LOG as default

        self.__Order()

# ######################################################

    def AddFuncVal(self,i,y):
        '''
        add the function value, y, for coordinate x[i].
        '''

        assert (i < len(self.x))
        self.Data[i]=y
        self.Ordered = False

# ######################################################

    def Insert(self,x,y):
        '''
        insert one function value in the table.
        x = coordinate point
        y = function value at x
        '''
        self.x.append(x)
        self.Data[len(self.x)-1]=y
        self.Ordered = False

# ######################################################

    def SetInterpolateType(self,interp):
        '''
        arguements are string values either:
        LINEAR
        -or-
        LOG (Default)
        '''
        self.Interp=interp

# ######################################################

    def __Order(self):
        '''
        order the x values and rearrange the data. 
        '''
        for i in range(len(self.x)-1):
            for j in range(i,len(self.x)): 
                if self.x[i] > self.x[j]:
                    # fix x
                    tmp = self.x[i]
                    self.x[i]=self.x[j]
                    self.x[j]=tmp

                    # fix Data
                    tmp = self.Data[i]
                    self.Data[i]=self.Data[j]
                    self.Data[j]=tmp

        self.Ordered=True

# ######################################################

    def IsInDomain(self,a):
        if len(self.x) == 0: return False
        if not self.Ordered: self.__Order()
        return (a >= self.x[0] and a <= self.x[-1])

# ######################################################

    def Evaluate(self,a):
        '''
        return the approximate (interpolated) value, f(a).
        '''

        # self.IsInDomain calls self.__Order()...
        if self.IsInDomain(a): pass
        else:
            assert (self.IsInDomain(a)) 

        #binary search to find the row segment
        l = 0
        u = len(self.x)-1
        while (1):
            m = (l+u)/2 
            if (self.x[m] < a): l = m
            elif (self.x[m] > a): u = m
            else: break # self.x[m] = a
            if (l == u-1): break 

        if (l == u-1):
            rl = l
            ru = u
        else:
            if (l == u): m = u
            if (m == len(self.x)-1):
                rl = m-1
                ru = m
            else:
                rl = m
                ru = m+1

        # interpolate 

        if (self.Interp == 'LINEAR'):
            r = (a - self.x[rl])/(self.x[ru]-self.x[rl])
            val = self.Data[rl] + \
                     r * (self.Data[ru] - self.Data[rl]) 
        else:
            r = (math.log10(a) - math.log10(self.x[rl])) / \
                (math.log10(self.x[ru]) - math.log10(self.x[rl])) 
            val = math.log10(self.Data[rl]) + \
                     r * (math.log10(self.Data[ru]) - \
                                        math.log10(self.Data[rl]))
            val = pow(10.0,val)

        return val

# ######################################################

if __name__ == "__main__":

    # the function:
    myfun = lambda a: a**0.5
    print " The function is: myfun = lambda a: a**0.5"
    # insert this value... see below
    p = 4.0

    x=[0.1,0.25,0.5,1.0,2.0,3.0,6.0]
    y=map(myfun,x)

    fun=OneVariable(x,y)

    print 'x values: ', fun.x
    print 'Data', fun.Data
    print ' '

    for i in range(len(x)-1):
        a = x[i] + 0.5
        print " Linear: %f, %f, % f" % (a, fun.Evaluate(a), myfun(a))
        fun.SetInterpolateType('LOG')
        print " Log   : %f, %f, % f" % (a, fun.Evaluate(a), myfun(a))
        fun.SetInterpolateType('LINEAR')

    print ' '
    fun.Insert(p,myfun(p))
    print 'x values: ', fun.x
    print 'Data', fun.Data
    print ' '

    for i in range(len(x)-1):
        a = x[i] + 0.5
        print " Linear: %f, %f, % f" % (a, fun.Evaluate(a), myfun(a))
        fun.SetInterpolateType('LOG')
        print " Log   : %f, %f, % f" % (a, fun.Evaluate(a), myfun(a))
        fun.SetInterpolateType('LINEAR')

    print ' '
    print 'x values: ', fun.x
    print 'Data', fun.Data
    print ' '







