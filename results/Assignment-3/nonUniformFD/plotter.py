import numpy as np
from matplotlib import pyplot as plt

infile=open("results.dat","r")
inlines=infile.readlines()

x=[]
d_exact=[]
d1=[]
d2=[]

for i in range(0,len(inlines)):
        inline=inlines[i].strip().split()
        x.append(float(inline[0]))
        d1.append(float(inline[1]))
        d2.append(float(inline[2]))
        d_exact.append(float(inline[3]))

plt.rcParams.update({'font.size': 12})

plt.figure(0)
plt.plot(x,d1,"*")
plt.plot(x,d2)
plt.plot(x,d_exact)
plt.grid()
plt.legend(["Non-uniform-grid-aware scheme","Standard central difference","Exact"])
plt.show()
