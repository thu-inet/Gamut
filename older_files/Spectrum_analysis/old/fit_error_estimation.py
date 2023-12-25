import numpy.random as rd
import numpy as np
import matplotlib.pyplot as plt
def Fiterror(error):
    sum_error=0
    sum_real=0
    for j in range(1000):
        for i in range(-10,10):
            h=1+rd.randn()*error
            sigma=1+rd.randn()*error
            sum_real+=np.exp(-i**2/2)**2
            sum_error+=(h*np.exp(-i**2/2/sigma**2)-np.exp(-i**2/2))**2
    return 1-sum_error/sum_real

def Fiterror2(error):
    N=100
    c=0.29317
    n=10
    sigma=3
    b=c*n
    list_e=list(rd.randn(N)*sigma)
    list_ind=np.linspace(-b-c,b+c,2*n+1)
    list_e.extend(list_ind)
    list_e.sort()
    list_ind_i=[]
    ind_old=-1
    for i in range(-n,n+1):
        ind_i=list_e.index(list_ind[i+n])
        l=ind_i-ind_old-1
        ind_old=ind_i
        list_ind_i.append(l)
    list_ind_i_thoe=[np.exp(-i**2/2/sigma**2)*N*c*2 for i in np.linspace(-b,b,2*n+1)]
    print(sum(list_ind_i))
    s1=0
    s2=0
    for i in range(len(list_ind_i)):
        s1+=list_ind_i_thoe[i]**2
        s2+=(list_ind_i[i]-list_ind_i_thoe[i])**2
    # print(1-s2/s1)
    plt.plot(list_ind_i)
    plt.plot(list_ind_i_thoe)
    plt.show()
    return list_ind_i
print(Fiterror2(0.04))
# for i in np.linspace(0,0.5,51):
#     print(Fiterror(i))