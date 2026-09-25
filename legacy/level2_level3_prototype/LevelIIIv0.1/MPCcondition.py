import string

if __name__=="__main__":

    b=open('MPC','r')

    buff=b.readline()
    weights=[]

    while buff:
        line=string.split(buff)
        buff=b.readline()
        for i in range(5,len(line),2): weights+=[float(line[i])]

    print abs(max(weights))/abs(min(weights))