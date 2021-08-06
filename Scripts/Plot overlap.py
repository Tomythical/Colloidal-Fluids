import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.signal import find_peaks
import time

# data = np.loadtxt('Data/2D/2D plot-0.68.dat')
# data1 = np.loadtxt('Data/2D/2D plot-{0.69}.dat')
# data2 = np.loadtxt('Data/2D/2D plot-{0.70}.dat')
# data3 = np.loadtxt('Data/2D/2D plot-{0.71}.dat')
# x = data[0]
# y = data[1]/20
# x1 = data1[0]
# y1 = data1[1]
# x2 = data2[0]
# y2 = data2[1]
# x3 = data3[0]
# y3 = data3[1]
# # x4 = data4[0]
# # y4 = data4[1]
# for i in range(len(y1)):
#     y1[i] = y1[i] + 0.5
# for i in range(len(y2)):
#     y2[i] = y2[i] + 1
# for i in range(len(y3)):
#     y3[i] = y3[i] + 1.5
# # for i in range(len(y4)):
# #     y4[i] = y4[i] + 1.5
# plt.plot(x,y)
# plt.plot(x1,y1)
# plt.plot(x2,y2)
# plt.plot(x3,y3)
# # plt.plot(x4,y4)
# plt.xlabel('Radius (r)')
# plt.ylabel('Pair Distribution Function g(r)')
# plt.legend(["Area Fraction = 0.68", "Area Fraction = 0.69", "Area Fraction = 0.70", "Area Fraction = 0.71"], loc='upper right')
# # plt.savefig("3D Mixed Plot(3).jpeg", dpi = 250)
# plt.show()
def d3():
    N3d = [296,320,344,360,372,384,400]
    y3d = []
    x3d = []
    for i in N3d:
        x3d.append((i*math.pi)/(6*(8**3)))
    af = 0.66
    for i in range(len(N3d)):
        data = np.loadtxt(f"Data/3D/3D plot-{N3d[i]}.dat")
        N2d = 100
        diameter = 1
        afL = math.sqrt(N2d*math.pi*(diameter**2)/(4*af))
        L = round(afL,3)
        print(L)
        avg_no_density = N2d/(L**2)
        print(avg_no_density)
        density1d = 0.75
        density3d = N3d[i]/(8**3)
        xvals = data[0]
        yvals = data[1]

        troughs,_ = find_peaks(-yvals,height=[-20,-0.001])

        # plt.plot(xvals,yvals)
        # plt.plot(xvals[troughs],yvals[troughs],'o')
        # plt.title("Area Fraction - 0.71")
        # plt.savefig('Data/2D/2D_Mixed_Plot.png', dpi =250)
        # plt.show()
        for i in range(len(xvals)):
            yvals[i] = (yvals[i])* (xvals[i]**2)
        # print(xvals[troughs[0]])
        # print(yvals[troughs[0]])
        # ytrough = troughs[0]

        cn = np.trapz(yvals[0:32],x=xvals[0:32])
        final = cn*4*math.pi*density3d
        y3d.append(final)

    print(y3d)
    plt.plot(x3d,y3d,marker= 'x')
    plt.xlabel('Packing Fraction')
    plt.ylabel('Coordination Number')
    # plt.errorbar(x3d,y3d,yerr=0.1)
    plt.show()

def d2():
    n2d = [0.66,0.67,0.68,0.69,0.7,0.71,0.72]
    y2d = []
    for i in n2d:
        data = np.loadtxt(f"Data/2D/2D plot-{i}.dat")
        N2d = 100
        diameter = 1
        afL = math.sqrt(N2d*math.pi*(diameter**2)/(4*i))
        L = round(afL,3)
        avg_no_density = N2d/(L**2)
        xvals = data[0]
        yvals = data[1]

        troughs,_ = find_peaks(-yvals,height=[-20,-0.001])

        plt.plot(xvals,yvals)
        plt.plot(xvals[troughs],yvals[troughs],'o')
        # plt.title("Area Fraction - 0.71")
        # plt.savefig('Data/2D/2D_Mixed_Plot.png', dpi =250)
        plt.show()
        for i in range(len(xvals)):
            yvals[i] = (yvals[i])* (xvals[i]**2)
        # print(xvals[troughs[0]])
        # print(yvals[troughs[0]])
        # ytrough = troughs[0]
        # print(troughs)
        cn = np.trapz(yvals[0:27],x=xvals[0:27])
        final = cn*4*math.pi*avg_no_density
        y2d.append(final)

    print(y2d)
    plt.plot(n2d,y2d,marker= 'x')
    plt.xlabel('Packing Fraction')
    plt.ylabel('Coordination Number')
    plt.show()

def d1():
    n2d = [0.75,0.8,0.85]
    y2d = []
    for i in n2d:
        data = np.loadtxt(f"Data/1D/1d plot-{i}.dat")
        avg_no_density = i
        xvals = data[0]
        yvals = data[1]

        troughs,_ = find_peaks(-yvals,height=[-20,-0.001])

        # plt.plot(xvals,yvals)
        # plt.plot(xvals[troughs],yvals[troughs],'o')
        # # plt.title("Area Fraction - 0.71")
        # # plt.savefig('Data/2D/2D_Mixed_Plot.png', dpi =250)
        # plt.show()
        for i in range(len(xvals)):
            yvals[i] = (yvals[i])* (xvals[i]**2)
        # print(xvals[troughs[0]])
        # print(yvals[troughs[0]])
        # ytrough = troughs[0]
        print(troughs)
        cn = np.trapz(yvals[0:7],x=xvals[0:7])
        final = cn*4*math.pi*avg_no_density
        y2d.append(final)

    print(y2d)
    plt.plot(n2d,y2d)
    plt.xlabel('Packing Fraction')
    plt.ylabel('Coordination Number')
    plt.show()

d1()