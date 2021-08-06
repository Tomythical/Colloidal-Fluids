import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import math
import copy
import time


N = 100
radius = 0.5
diameter = radius*2
nparticles = []
plt.figure(figsize=(6,6))
dr = 0.05
af = 0.66
afL = math.sqrt(N*math.pi*(diameter**2)/(4*af))
area_fraction = (N*math.pi*(diameter**2))/(4*(afL**2))
L = round(afL,3)
print(L)
avg_no_density = N/L**2
print(avg_no_density)
states = int(5/dr)
Nrk = np.zeros(states)

class Particle:
    def __init__(self, x,y, radius): 
        self.position = [x,y]
        self.radius = radius
    
    
def init_particles():
    rem = N%int(L)
    base = N//int(L)
    if L**2 < N:
        raise Exception
    for j in np.arange(0,int(L),diameter):
        for k in np.arange(0,base,diameter):
            nparticles.append(Particle((j+radius),(k+radius), radius))
    for l in np.arange(0,int(rem),diameter):
        nparticles.append(Particle((l+radius),(base+radius), radius))


def get_mc_move(xamp,yamp):
    dx = np.random.uniform(-0.5,1.5)
    dy = np.random.uniform(-0.5,1.5)
    xstep = xamp * (dx - 0.5)
    ystep = yamp * (dy - 0.5)

    return [xstep,ystep]
    
"""Chooses a random particle and changes the particle's position if the move is allowed"""
def particle_bc(nparticles, move):              #nparticles is a list of particle objects
    rnd = np.random.randint(0,N)
    rnd_particle = nparticles[rnd]
    new_pos = [move[0]+ rnd_particle.position[0],move[1]+ rnd_particle.position[1]]  
    # rnd_particle.position = new_pos
    if new_pos[0] > L:
        new_pos[0] = new_pos[0] - L
    if new_pos[0] < 0:
        new_pos[0] = new_pos[0] + L
    if new_pos[1] > L:
        new_pos[1] = new_pos[1] - L
    if new_pos[1] < 0:
        new_pos[1] = new_pos[1] + L
    temp_part =  copy.deepcopy(rnd_particle)   #Need to copy particle into temporary variable
    temp_part.position = new_pos
    return (temp_part, rnd_particle, rnd)

def overlaps(p1, p2):  # Checks if two particles overlap
    intersects = math.hypot(p1.position[0]-p2.position[0], p1.position[1] - p2.position[1]) <= (p1.radius + p2.radius)
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
            if i.position[0] > temp_particle.position[0] - diameter and i.position[0] < temp_particle.position[0] + diameter and i.position[1] > temp_particle.position[1] - diameter and i.position[1] < temp_particle.position[1] + diameter:
                if overlaps(temp_particle,i) is True:
                    x = True
                    break
            if i.position[0] < diameter or i.position[0] > L-diameter or i.position[1] < diameter or i.position[1] > L-diameter:
                if temp_particle.position[0] < diameter and temp_particle.position[1] > diameter and temp_particle.position[1] < L -  diameter: #for x<0.5
                    temp_particle.position[0] += L
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[0] -= L
                if temp_particle.position[0] > L - diameter and temp_particle.position[1] > diameter and temp_particle.position[1] < L -  diameter: #for x>L-0.5
                    temp_particle.position[0] -= L
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[0] += L
                if temp_particle.position[1] < diameter and temp_particle.position[0] > diameter and temp_particle.position[0] < L -  diameter: #for y < 0.5
                    temp_particle.position[1] += L
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[1] -= L
                if temp_particle.position[1] > L - diameter and temp_particle.position[0] > diameter and temp_particle.position[0] < L -  diameter: #for y> L- 0.5
                    temp_particle.position[1] -= L
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[1] += L

                if temp_particle.position[0] < diameter and temp_particle.position[1] < diameter: #for x<0.5 and y < 0.5 i.e Bottom-left
                    temp_particle.position[0] += L                                            # Check overlap with top-right
                    temp_particle.position[1] += L
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[0] -= L
                    temp_particle.position[1] -= L

                    temp_particle.position[1] += L                                            # Check overlap with top-left
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[1] -= L

                    temp_particle.position[0] += L                                            # Check overlap with bottom-right
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[0] -= L

                if temp_particle.position[0] < diameter and temp_particle.position[1] > L - diameter: #for x<0.5 and y > L -0.5 i.e Top-left
                    temp_particle.position[0] += L                                                # Check overlap with bottom-right
                    temp_particle.position[1] -= L
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[0] -= L
                    temp_particle.position[1] += L

                    temp_particle.position[0] += L                                            # Check overlap with top-right
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[0] -= L

                    temp_particle.position[1] -= L                                            # Check overlap with bottom-left
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[1] += L

                if temp_particle.position[0] > L - diameter and temp_particle.position[1] > L - diameter: #for x>L-0.5 and y > L - 0.5 i.e. Top-Right
                    temp_particle.position[0] -= L                                                     # Check overlap with bottom-left
                    temp_particle.position[1] -= L
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[0] += L
                    temp_particle.position[1] += L
                
                    temp_particle.position[0] -= L                                            # Check overlap with top-left
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[0] += L

                    temp_particle.position[1] -= L                                            # Check overlap with bottom-right
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[1] += L

                if temp_particle.position[0] > L - diameter and temp_particle.position[1] < diameter: #for x> L - 0.5 and y < 0.5 i.e. Bottom-right
                    temp_particle.position[0] -= L                                                  # Check overlap with top-left
                    temp_particle.position[1] += L
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[0] += L
                    temp_particle.position[1] -= L
                    
                    temp_particle.position[1] += L                                            # Check overlap with top-right
                    if overlaps(temp_particle, i) is True:
                        x = True
                        break
                    temp_particle.position[1] -= L

                    temp_particle.position[0] -= L                                            # Check overlap with bottom-left
                    if overlaps(temp_particle, i) is True:
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
        ydist = particle.position[1] - ref.position[1]
        if xdist < -L/2:
            xdist = xdist + L
        if xdist > L/2:
            xdist = xdist - L
        if ydist < -L/2:
            ydist = ydist + L
        if ydist > L/2:
            ydist = ydist - L
        rdist = math.sqrt(xdist**2 + ydist**2)
        if rdist >= (k*dr) and rdist <= (k*dr+dr):
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
        gr.append((Nrk[i]/50000)/(2*math.pi*i*dr*dr*avg_no_density))
    return gr


def init_plot():
    x = []
    y = []
    for i in range(len(nparticles)):
        x.append(nparticles[i].position[0])
        y.append(nparticles[i].position[1])
    return(x,y)


def graph(x,y):
    ax = plt.axes(xlim=(0,L),ylim=(0,L))
    scatter = ax.scatter(x, y, s =900)
    plt.show()
    
def graph2(a,b):
    plt.plot(a,b)
    plt.axhline(y=1,linestyle='dashed',color='black')
    # plt.title(f"Area Fraction - {afL}")
    plt.xlabel('Radius (r)')
    plt.ylabel('Pair Distribution Function g(r)')
    plt.show()

def repeat():
    xamp = 1
    yamp = 1
    for j in range(500):
        moves_accepted = 0
        for i in range(1000):
            trial_part = particle_bc(nparticles,get_mc_move(xamp,yamp))  #returns trial particle with a new position and original particle
            moves_accepted += check_mc_move(trial_part[0],trial_part[1], trial_part[2])
        if moves_accepted/1000 < 0.25:
            # print(moves_accepted)
            xamp = xamp * 0.95
            yamp = yamp * 0.95
        if moves_accepted/1000 > 0.5:
            # print(moves_accepted)
            xamp = xamp * 1.05
            yamp = yamp * 1.05
        count_shell_radii()
        print(j)



if __name__ == '__main__':
    start_time = time.time()
    init_particles()
    xamp = 1
    yamp = 1
    for i in range (100):
        moves_accepted = 0
        for i in range(1000):
            trial_part = particle_bc(nparticles,get_mc_move(xamp,yamp))  #returns trial particle with a new position and original particle
            moves_accepted += check_mc_move(trial_part[0],trial_part[1], trial_part[2])  
        if moves_accepted/1000 < 0.25:
            xamp = xamp * 0.95
            yamp = yamp * 0.95
        if moves_accepted/1000 > 0.5:
            xamp = xamp * 1.05
            yamp = yamp * 1.05


    repeat()
    
    executionTime = (time.time() - start_time)
    graph(init_plot()[0],init_plot()[1])
    rx = []
    for i in range(1,states):
        rx.append(dr*i)
    graph2(rx,calc_gr())
    np.savetxt(f"/Users/ThomasMatheickal/Projects/Colloidal Fluids/Data/2D/2D plot-{area_fraction}.dat", [rx,calc_gr()])
    print('Execution time in minutes: ' + str(executionTime/60)) 

