import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import copy
import time

L = 8
N = 384
radius = 0.5
diameter = radius*2
nparticles = []
dr = 0.05
area_fraction = (N*math.pi*(diameter**3))/(6*(L**3))
avg_no_density = N/L**3
states = int(4/dr)
Nrk = np.zeros(states)

class Particle:
    def __init__(self, x,y,z,radius): 
        self.position = [x,y,z]
        self.radius = radius
    
    
def init_particles():
    rem2d = N%int(L)
    base3d = N//(int(L)**2)
    rem3d = N%(int(L)**2)
    rem3d2 = rem3d//L
    for x in np.arange(0,int(L),diameter):
        for y in np.arange(0,int(L),diameter):
            for z in np.arange(0,base3d,diameter):
                nparticles.append(Particle((x+radius),(y+radius),(z+radius),radius))
    for j in np.arange(0, int(L),diameter):
        for k in np.arange(0,rem3d2,diameter):
                nparticles.append(Particle((j+radius),(k+radius), base3d+radius,radius))
    for l in np.arange(0,rem2d,diameter):
        nparticles.append(Particle((l+radius),radius,base3d+radius, radius))

def graph(x,y,z):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlim3d(0, L)
    ax.set_ylim3d(0,L)
    ax.set_zlim3d(0, L)
    ax.scatter(x, y, z ,depthshade=1,s =1500)
    # ax.set_edgecolors = ax.set_facecolors = lambda *args:None
    plt.show()

def init_plot():
    x = []
    y = []
    z = []
    for i in range(len(nparticles)):
        x.append(nparticles[i].position[0])
        y.append(nparticles[i].position[1])
        z.append(nparticles[i].position[2])
    return(x,y,z)

print(area_fraction)
init_particles()
print(len(nparticles))
# graph(init_plot()[0],init_plot()[1],init_plot()[2])

