import numpy as np
import matplotlib.pyplot as plt

def get_nboxes(N):
    N = int(N)
    m=1
    nt=N
    while nt!=0:
        nt = nt >> 1
        m = m << 1
    nboxes = 2*N - m<<1
    if m<nboxes:
        nboxes = m
    return nboxes, 
#    return "N = "+str(N)+"  nboxes = "+str(nboxes-1)
Nbs = []
Ns =np.linspace(10,10**9,100000)
for n in Ns:
    Nbs.append(get_nboxes(n))
Nbs=np.array(Nbs)
plt.figure(1)
plt.plot(np.log10(Ns),Nbs/Ns)
#plt.plot(np.log10(Ns),np.log10(Nbs))
#plt.plot(Ns,Ns) 
plt.show()
