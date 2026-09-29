from matplotlib import pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

solverID=1 #CHANGE THIS TO 0 FOR FIRST ORDER UPWIND SCHEME, 1 FOR COMPACT SCHEME

fig,ax=plt.subplots()
plt.grid()
plt.xlabel("x")
plt.ylabel("u")
plt.ylim(-1.0,1.0)
plt.xlim(0.0,1.0)
line, = ax.plot([], [], lw=2)

infile=open("s_"+str(solverID)+"_"+"dump_"+str(0)+".dat")
infile.readline()
infile.readline()
inline=infile.readline().strip().split()
imx=int(inline[2][2:])

x=np.zeros(imx)
sol=np.zeros(imx)

def update(frame):
        infile=open("s_"+str(solverID)+"_"+"dump_"+str(frame)+".dat")
        infile.readline()
        infile.readline()
        inline=infile.readline().strip().split()
        infile.readline()
        infile.readline()

        for i in range(0,imx):
                inline=infile.readline().strip().split()
                x[i]=float(inline[0])
        
        for i in range(0,imx):
                inline=infile.readline().strip().split()
                sol[i]=float(inline[0])

        line.set_data(x,sol)
        infile.close()
        return line,

anim=FuncAnimation(
        fig,
        func=update,
        frames=100,
        interval=100,
        blit=False
        )

plt.show()

writer=PillowWriter(fps=20)
anim.save("animation"+str(solverID)+".gif",writer=writer)


#Comparison of dispersive error

deltaX=x[1]-x[0]
dt=1.0*deltaX
uExact=np.sin(np.pi*(x-1.0*dt*10000)/(10.0*deltaX));


plt.figure()
infile=open("s_"+str(solverID)+"_"+"dump_"+str(9999)+".dat")
infile.readline()
infile.readline()
inline=infile.readline().strip().split()
infile.readline()
infile.readline()

for i in range(0,imx):
        inline=infile.readline().strip().split()
        x[i]=float(inline[0])

for i in range(0,imx):
        inline=infile.readline().strip().split()
        sol[i]=float(inline[0])

infile.close()

plt.rcParams.update({'font.size': 12}) 
plt.plot(x,sol)
plt.plot(x,uExact,"*")
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.title("Solution error: 10K iters@CFL 1")
plt.legend(["Numerical","Analytical"])
plt.show()
