import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import math
import copy
import time


N = 85
L = 100
radius = 0.5
diameter = radius*2
nparticles = []
plt.figure(figsize=(6,6))
dr = 0.2
area_fraction = N/L
avg_no_density = N/L
states = int((L/4)/dr)
Nrk = np.zeros(states)
repeat_N = 400

class Particle:
    def __init__(self, x, radius): 
        self.position = [x]
        self.radius = radius
    
    
def init_particles():
    if L < N:
        raise Exception
    for j in np.arange(0,N,diameter):
        nparticles.append(Particle((j+radius), radius))
    

def get_mc_move(xamp):
    dx = np.random.uniform(-0.5,1.5)
    xstep = xamp * (dx - 0.5)
    return [xstep]
    
"""Chooses a random particle and changes the particle's position if the move is allowed"""
def particle_bc(nparticles, move):              #nparticles is a list of particle objects
    rnd = np.random.randint(0,N)
    rnd_particle = nparticles[rnd]
    new_pos = [move[0]+ rnd_particle.position[0]]
    # rnd_particle.position = new_pos
    if new_pos[0] > L:
        new_pos[0] = new_pos[0] - L
    if new_pos[0] < 0:
        new_pos[0] = new_pos[0] + L
    temp_part =  copy.deepcopy(rnd_particle)   #Need to copy particle into temporary variable
    temp_part.position = new_pos
    return (temp_part, rnd_particle, rnd)

def overlaps(p1, p2):  # Checks if two particles overlap
    # intersects = math.hypot(p1.position[0]-p2.position[0], p1.position[1] - p2.position[1]) <= (p1.radius + p2.radius)
    intersects = abs(p1.position[0]-p2.position[0]) <= (p1.radius + p2.radius)
    if intersects:
        return True # Circle radii intersects
    else:
        return False

""" Change radius to diameter"""
def check_mc_move(temp_particle,org_particle, rnd): # Checks if the trial particle overlaps with any other particle considering boundary conditions
    x = None
    moves_accepted = 0
    for i in nparticles:
        if i is not org_particle:
            if overlaps(temp_particle, i) is True:
                x = True
                break

            if temp_particle.position[0] < diameter:
                temp_particle.position[0] += L
                if overlaps(temp_particle,i) is True:
                    x = True
                    break
                temp_particle.position[0] -= L

            if temp_particle.position[0] > L - diameter:
                temp_particle.position[0] -= L
                if overlaps(temp_particle,i) is True:
                    x = True
                    break
                temp_particle.position[0] += L

    if x is None:
        org_particle.position = temp_particle.position
        nparticles[rnd] = org_particle
        moves_accepted = 1
    return moves_accepted


def count_particles2(ref, particle,k):              # Checks if a particle is within the shell radius originated from ref
    if ref != particle:
        xdist = particle.position[0] - ref.position[0]
        if xdist < -L/2:
            xdist = xdist + L
        if xdist > L/2:
            xdist = xdist - L
        if xdist >= (k*dr) and xdist <= (k*dr+dr):
            return True
        else:
            return False

def count_shell_radii():
    for k in range(states):
        for i in nparticles:
            for j in nparticles:
                if count_particles2(i,j,k):
                    Nrk[k] += 1

def calc_gr():
    gr = []
    for i in range (1,len(Nrk)):
        gr.append(Nrk[i]/((avg_no_density)*(N*dr*repeat_N)))
    return gr


def init_plot():
    x = []
    for i in range(len(nparticles)): 
        x.append(nparticles[i].position[0])
    return(x)


def graph(x):
    # plt.figure(figsize=(8, 0.5))
    # scatter = ax.scatter(x, s =900)
    # plt.savefig('1D scatter.jpeg', dpi = 250)
    plt.figure()
    plt.axes(xlim=(0,L))
    # plt.figure(figsize=(8, 0.5))
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.scatter(x, len(x)*[1], s=130)
    plt.show()
    
def graph2(a,b):
    plt.plot(a,b)
    # plt.axhline(y=1,linestyle='dashed',color='black')
    # plt.title(f"Area Fraction - {afL}")
    plt.xlabel('Radius (r)')
    plt.ylabel('Pair Distribution Function g(r)')
    # plt.savefig('1D plot(4).jpeg', dpi = 250)
    plt.show()

def repeat():
    xamp = 1
    for j in range(repeat_N):
        moves_accepted = 0
        for i in range(1000):
            trial_part = particle_bc(nparticles,get_mc_move(xamp))  #returns trial particle with a new position and original particle
            moves_accepted += check_mc_move(trial_part[0],trial_part[1], trial_part[2])
        if moves_accepted/1000 < 0.25:
            # print(moves_accepted)
            xamp = xamp * 0.95
        if moves_accepted/1000 > 0.5:
            # print(moves_accepted)
            xamp = xamp * 1.05
        count_shell_radii()
        print(j)



if __name__ == '__main__':
    start_time = time.time()
    init_particles()
    xamp = 1
    for i in range (100):
        moves_accepted = 0
        for i in range(1000):
            trial_part = particle_bc(nparticles,get_mc_move(xamp))  #returns trial particle with a new position and original particle
            moves_accepted += check_mc_move(trial_part[0],trial_part[1], trial_part[2])  
        if moves_accepted/1000 < 0.25:
            xamp = xamp * 0.95
        if moves_accepted/1000 > 0.5:
            xamp = xamp * 1.05


    repeat()
    # for i in range(len(nparticles)):
    #     print(nparticles[i].position)
    executionTime = (time.time() - start_time)
    graph(init_plot())
    rx = []
    for i in range(1,states):
        rx.append(dr*i)
    graph2(rx,calc_gr())
    np.savetxt(f"/Users/ThomasMatheickal/Projects/Colloidal Fluids/Data/1D/1d plot-{area_fraction}.dat", [rx,calc_gr()])
    print("area fraction: " + str(area_fraction))
    print('Execution time in minutes: ' + str(executionTime/60)) 

