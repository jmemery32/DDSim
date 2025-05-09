import string

def MPC(fin,fout,inc):
    '''
    change node id's of global nodes in mpc file
    '''

    # the column in the MPC file that tells how many mpc's are coming...
    # currently column #3... but file formats are always subject to change.
    numCol=3 

    buff=fin.readline()
    while buff:
        list=string.split(buff)
        buff=fin.readline()

        number=int(list[numCol])
        for i in range(number):
            list[numCol+1+(2*i)] = int(list[numCol+1+(2*i)]) + inc

        s=''
        for p in list: s+=str(p)+' '

        # strip off the last space
        s=s[0:-1]+"\n"

        fout.write(s)

if __name__=="__main__":

    nodeinc = 1932129
    elementinc = 1416088

    fin=open('constraints.mpc','r')
    fout=open('MPC.polycrystal','w')
    MPC(fin,fout,nodeinc)






