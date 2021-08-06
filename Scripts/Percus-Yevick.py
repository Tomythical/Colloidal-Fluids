import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import math

L = 8
N = 300
radius = 0.5
diameter = radius*2
dr = 0.05
area_fraction = (N*math.pi*(diameter**3))/(6*(L**3))
avg_no_density = N/L**3
states = int(1/dr)

p1 = ((1+2*area_fraction)**2)/((1-area_fraction)**4)
p2 = (-(1+area_fraction/2)**2)/((1-area_fraction)**4)

sf = []
qr = []

def calc_dcf(qd):
    ncq1 = p1*((math.sin(qd)-qd*math.cos(qd))/qd**3)
    ncq2 = -6*area_fraction*p2*((qd**2*math.cos(qd)-2*qd*math.sin(qd)-2*math.cos(qd)+2)/qd**4)  
    ncq3 = -area_fraction*p1/2*((qd**4*math.cos(qd)-4*qd**3*math.sin(qd)-12*qd**2*math.cos(qd)+24*qd*math.sin(qd)+24*math.cos(qd)-24)/qd**6)
    ncq = -24*area_fraction*(ncq1+ncq2+ncq3)
    sq = 1/(1-ncq)
    return sq

for i in np.arange(dr,30,0.01):
    qr.append(i/2)
    sf.append(calc_dcf(i))

print(area_fraction)
plt.plot(qr,sf)
plt.xlabel('QR')
plt.ylabel('Structure Factor(Q)')
plt.savefig('Percus-Yevick.jpeg', dpi =250)
plt.show()