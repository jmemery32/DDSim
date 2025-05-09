import math


class InterpTable:

    def __init__(self,rowvals,colvals):
        '''
        rowVals, colvals = lists of the values of rows and columns
        '''
        self.RowVals=rowvals.sort()
        self.ColVals=colvals.sort()
        self.Ordered=False 
        self.Data={} # use tuble
        self.RInterp='LINEAR' # set Linear as default
        self.CInterp='LINEAR' # set Linear as default

# ######################################################

    def AppendRowsAndInsert(self,new_rowvals,new_fxy):  #not tested
        '''
        len(new_fxy) = len(new_rowvals)*len(self.ColVals) AND:

        new_fxy[0:len(self.ColVals)] = f(new_rowval[0],self.ColVals[0:]) and 
        new_fxy[len(self.ColVals):2*len(self.ColVals)] = \
                                       f(new_rowval[1],self.ColVals[0:]) and
        new_fxy[2*len(self.ColVals):3*len(self.ColVals)] = \
                                       f(new_rowval[2],self.ColVals[0:]) etc...
        '''

        lc=len(self.ColVals)
        for i in range(len(new_rowvals)):
            self.AppendRowAndInsert(new_rowvals[i],\
                        new_fxy[i*lc:(i+1)*lc]

# ######################################################

    def AppendRowAndInsert(self,new_rowval,new_fxy): #not tested
        '''
        new_rowval = float
        new_fxy = list, len(new_fxy) == len(self.ColVals)
        '''
        assert(len(new_fxy) == len(self.ColVals))
        for i in xrange(len(self.ColVals)):
            self.Data[(len(self.RowVals),i)] = new_fxy[i]

        self.AppendRow(new_rowval)

# ######################################################

    def AppendColsAndInsert(self,new_colvals,new_fxy): #not tested
        '''
        len(new_fxy) = len(new_colvals)*len(self.RowVals) AND:

        new_fxy[0:len(self.RowVals)] = f(new_colval[0],self.RowVals[0:]) and 
        new_fxy[len(self.RowVals):2*len(self.RowVals)] = \
                                       f(new_colval[1],self.RowVals[0:]) and
        new_fxy[2*len(self.RowVals):3*len(self.RowVals)] = \
                                       f(new_colval[2],self.RowVals[0:]) etc...
        '''

        lr=len(self.RowVals)
        for i in range(len(new_colvals)):
            self.AppendColAndInsert(new_colvals[i],\
                        new_fxy[i*lr:(i+1)*lr]

# ######################################################

    def AppendColAndInsert(self,new_colval,new_fxy): #not tested
        '''
        new_colval = float
        new_fxy = list, len(new_fxy) == len(self.RowVals)
        '''
        assert(len(new_fxy) == len(self.RowVals))
        for i in xrange(len(self.RowVals)):
            self.Data[(i,len(self.ColVals))] = new_fxy[i]

        self.AppendCol(new_colval)

# ######################################################

    def AppendRow(self,rowval): #not tested
        self.Ordered = False
        self.RowVals+=[rowval]

# ######################################################


    def AppendCol(self,colval): #not tested
        self.Ordered = False
        self.ColVals += colval

# ######################################################

    def AppendRows(self,rowvals): #not tested
        self.Order = False
        self.RowVals.append(rowvals)

# ######################################################

    def AppendCols(self,colvals): #not tested
        self.Order = False
        self.ColVals.append(colvals)

# ######################################################

    def Insert(self,row,col,fxy):
        '''
        insert one function value in the table.
        row = int
        col = int
        '''
        assert (row < len(self.RowVals) and col < len(self.ColVals))
        self.Data[(row,col)]=fxy

# ######################################################

    def SetInterpolateType(self,row_interp,col_interp):
        '''
        arguements are string values either:
        LINEAR
        -or-
        LOGARITHMIC (but any string value other than LINEAR is treated as Log!
        '''
        self.RInterp=row_interp
        self.CInterp=col_interp

# ######################################################

    def __Order(self):
        # sort rows first... 
        for i in range(len(self.RowVals)-1):
            for j in range(i,len(self.RowVals)): 
                if self.RowVals[i] > self.RowVals[j]:
                    # fix RowVals
                    tmp = self.RowVals[i]
                    self.RowVals[i]=self.RowVals[j]
                    self.RowVals[j]=tmp

                    # fix Data
                    for k in range(len(self.ColVals)):
                        tmp = self.Data[(i,k)]
                        self.Data[(i,k)]=self.Data[(j,k)]
                        self.Data[(j,k)]=tmp

        # now sort cols... 
        for i in range(len(self.ColVals)-1):
            for j in range(i,len(self.ColVals)): 
                if self.ColVals[i] > self.ColVals[j]:
                    # fix ColVals
                    tmp = self.ColVals[i]
                    self.ColVals[i]=self.ColVals[j]
                    self.ColVals[j]=tmp

                    # fix Data
                    for k in range(len(self.RowVals)):
                        tmp = self.Data[(k,i)]
                        self.Data[(k,i)]=self.Data[(k,j)]
                        self.Data[(k,j)]=tmp

        self.Ordered=True

# ######################################################

    def IsInRowDomain(self,row):
        if len(self.RowVals) == 0: return False
        if not self.Ordered: self.__Order()
        return (row >= self.RowVals[0] and row <= self.RowVals[-1])

# ######################################################

    def IsInColDomain(self,col):
        if len(self.ColVals) == 0: return False
        if not self.Ordered: self.__Order()
        return (col >= self.ColVals[0] and col <= self.ColVals[-1])

# ######################################################

    def Evaluate(self,row,col):
        assert (self.IsInRowDomain(row))
        assert (self.IsInColDomain(col))
        if not self.Ordered: self.__Order()

        #binary search to find the row segment
        l = 0
        u = len(self.RowVals)-1
        while (1):
            m = (l+u)/2 
            if (self.RowVals[m] < row): l = m
            elif (self.RowVals[m] > row): u = m
            else: break
            if (l == u-1): break 

        if (l == u-1):
            rl = l
            ru = u
        else:
            if (l == u): m = u
            if (m == len(self.RowVals)-1):
                rl = m-1
                ru = m
            else:
                rl = m
                ru = m+1

        #binary search to find the col segment
        l = 0
        u = len(self.ColVals)-1
        while (1):
            m = (l+u)/2 
            if (self.ColVals[m] < col): l = m
            elif (self.ColVals[m] > col): u = m
            else: break
            if (l == u-1): break 

        if (l == u-1):
            cl = l
            cu = u
        else:
            if (l == u): m = u
            if (m == len(self.ColVals)-1):
                cl = m-1
                cu = m
            else:
                cl = m
                cu = m+1

        # interpolate rows

        if (self.RInterp == 'LINEAR'):
            r = (row - self.RowVals[rl])/(self.RowVals[ru]-self.RowVals[rl])
            cv_low = self.Data[rl,cl] + \
                     r * (self.Data[ru,cl] - self.Data[rl,cl]) 
            cv_upr = self.Data[rl,cu] + \
                     r * (self.Data[ru,cu] - self.Data[rl,cu]) 
        else:
            r = (row - math.log10(self.RowVals[rl])) / \
                (math.log10(self.RowVals[ru]) - math.log10(self.RowVals[rl])) 
            cv_low = self.Data[rl,cl] + \
                     math.pow(10.0,r * (math.log10(self.Data[ru,cl]) - \
                                        math.log10(self.Data[rl,cl]))) 
            cv_upr = self.Data[rl,cu] + \
                     math.pow(10.0,r * (math.log10(self.Data[ru,cu]) - \
                                        math.log10(self.Data[rl,cu]))) 

        if (self.CInterp == 'LINEAR'):
            c = (col - self.ColVals[cl]) / (self.ColVals[cu]-self.ColVals[cl]) 
            val = cv_low + c * (cv_upr - cv_low) 
        else:
            c = (col - math.log10(self.ColVals[cl])) / \
                (math.log10(self.ColVals[cu]) - math.log10(self.ColVals[cl]))
            val = cv_low + \
                     math.pow(10.0,c*(math.log10(cv_upr) - math.log10(cv_low)))

        return val

# ######################################################

if __name__ == "__main__":

    row_vals=[]
    row_vals.append(1.0) 
    row_vals.append(2.0) 
    row_vals.append(3.0)
    row_vals.append(4.0) 
    row_vals.append(5.0) 

    col_vals=[]
    col_vals.append(10.0) 
    col_vals.append(20.0) 
    col_vals.append(30.0) 
    col_vals.append(40.0) 

    table = InterpTable(row_vals,col_vals)

    for i in range(len(row_vals)):
        for j in range(len(col_vals)):
            table.Insert(i,j,row_vals[i]*col_vals[j]) 

    print table.RowVals
    print ' [1.0, 2.0, 3.0, 4.0, 5.0]'
    print table.ColVals
    print ' [10.0, 20.0, 30.0, 40.0]'
    print table.Data
    print ' {(0, 1): 20.0, (1, 2): 60.0, (3, 2): 120.0, (0, 0): 10.0, '
    print ' (3, 3): 160.0, (3, 0): 40.0, (3, 1): 80.0, (2, 1): 60.0, '
    print ' (0, 2): 30.0, (2, 0): 30.0, (1, 3): 80.0, (2, 3): 120.0, '
    print ' (4, 3): 200.0, (2, 2): 90.0, (1, 0): 20.0, (4, 2): 150.0, '
    print ' (0, 3): 40.0, (4, 1): 100.0, (1, 1): 40.0, (4, 0): 50.0} '

    print "f(1.5,15.0): ", table.Evaluate(1.5,15.0)
    print ' f(1.5,15.0):  22.5'
    print "f(4.5,15.0): ", table.Evaluate(4.5,15.0)
    print ' f(4.5,15.0):  67.5'
    print "f(1.5,35.0): ", table.Evaluate(1.5,35.0)
    print ' f(1.5,35.0):  52.5'
    print "f(4.75,30.0): ", table.Evaluate(4.75,30.0)
    print ' f(4.75,30.0):  142.5'
    print "f(4.75,40.0): ", table.Evaluate(4.75,40.0)
    print ' f(4.75,40.0):  190.0'
    print "f(4.75,37.5): ", table.Evaluate(4.75,37.5)
    print ' f(4.75,37.5):  178.125'

##    new_rows=[10.0,1.5]
##    new_cols=[15.0,100.0,5.0]
##
##    table.AppendRows(new_rows)
##    table.AppendCols(new_cols)
##    
##    for i in range(len(row_vals)):
##        for j in range(len(col_vals)):
##            table.Insert(i,j,row_vals[i]*col_vals[j]) 
    
    
    
    
    
    








