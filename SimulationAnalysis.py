from matplotlib import pyplot as plt
import numpy as np
from SolarSystem import SolarSystem
from Particle import Particle
from os import path

"""
This file creates the graphs from the npy data, to access one graph plotting simply change the doc strings.
"""


fig = plt.figure(figsize=(8,8)) #was 8,8 for solar system part and 8,6 for other graphs

#GET RID OF THIS FOR THE ENERGY SYSTEM PART
ax = fig.add_subplot(111, projection='3d')
from mpl_toolkits.mplot3d import Axes3D

Data = np.load("/Users/annabellecarver/Desktop/Common/86400_10000stepECWJ.npy", allow_pickle=True)
Data_T = Data.T
#print(Data_T[1])


for i in range(len(Data[0][1].the_planets)):
    x_i = [q.the_planets[i].position[0] for q in Data_T[1]]
    y_i = [q.the_planets[i].position[1] for q in Data_T[1]]
    z_i = [q.the_planets[i].position[2] for q in Data_T[1]]
    planet_name = Data_T[1][0].the_planets[i].name
    ax.plot(x_i, y_i, z_i, label='{0} Trajectory'.format(planet_name))
ax.set_xlabel('x-position (m)', fontfamily="Monospace", fontsize=10)
ax.set_ylabel('y-position (m)', fontfamily="Monospace", fontsize=10)
ax.set_zlabel('z-position (m)', fontfamily="Monospace", fontsize=10)
plt.legend(loc = 'upper left')
plt.show()



"""


#this is for testing individual pieces of code before putting them together!
Data = np.load("/Users/annabellecarver/Desktop/Common/1e44_EC_AM.npy", allow_pickle=True)

x= []
y= []

for i in range(len(Data)):
    y.append(Data[i][1])
    x.append(Data[i][0])

plt.plot(x, y, label = "Euler")
#plt.title("Total Kinetic Energy")
plt.xlabel("Time (s)")
plt.ylabel("Total Angular Momentum (kgm^2/s)")
plt.show()
axes = plt.gca()
axes.set_xlim([0,2.25e43])


"""

#THIS IS THE GRAPH MAKING FOR ENERGY PLOTS
#if path.exists("smallstepEuler1.npy"):
#    print("The file exists.")

#THIS IS FOR KINETIC ENERGY GRAPHS
euler = np.load("/Users/annabellecarver/Desktop/Common/1e5_4total_E_KE.npy", allow_pickle=True)
euler_cromer = np.load("/Users/annabellecarver/Desktop/Common/1e5_4total_EC_KE.npy", allow_pickle=True)
Verlet = np.load("/Users/annabellecarver/Desktop/Common/1e5_4total_V_KE.npy", allow_pickle=True)




x1=[]
y1=[]

x2=[]
y2=[]

x3=[]
y3=[]

for i in range(len(euler)):
    y1.append(euler[i][1])
    x1.append(euler[i][0])

for i in range(len(euler_cromer)):
    y2.append(euler_cromer[i][1])
    x2.append(euler_cromer[i][0])

for i in range(len(Verlet)):
    y3.append(Verlet[i][1])
    x3.append(Verlet[i][0])


plt.plot(x1, y1, label = "Euler")
plt.plot(x2, y2, label = "Euler-Cromer")
plt.plot(x3, y3, label = "Verlet")
#plt.title("Total Kinetic Energy")
plt.xlabel("Time (s)")
plt.ylabel("Total Kinetic Energy (J)")
plt.legend(loc = 'upper right')
axes = plt.gca()
axes.set_xlim([0,10000000])
plt.show()






#THIS IS FOR POTENTIAL ENERGY GRAPHS
Euler = np.load("/Users/annabellecarver/Desktop/Common/1e5_4total_E_PE.npy", allow_pickle=True)
EulerCromer = np.load("/Users/annabellecarver/Desktop/Common/1e5_4total_EC_PE.npy", allow_pickle=True)
Verlet = np.load("/Users/annabellecarver/Desktop/Common/1e5_4total_V_PE.npy", allow_pickle=True)


x1=[]
y1=[]

x2=[]
y2=[]

x3=[]
y3=[]

for i in range(len(Euler)):
    y1.append(Euler[i][1])
    x1.append(Euler[i][0])

for i in range(len(EulerCromer)):
    y2.append(EulerCromer[i][1])
    x2.append(EulerCromer[i][0])

for i in range(len(Verlet)):
    y3.append(Verlet[i][1])
    x3.append(Verlet[i][0])


plt.plot(x1, y1, label = "Euler")
plt.plot(x2, y2, label = "Euler-Cromer")
plt.plot(x3, y3, label = "Verlet")
#plt.title("Total Gravitational potential Energy")
plt.xlabel("Time (s)")
plt.ylabel("Total Gravitational Potential Energy (J)")
plt.legend(loc = 'upper right')
axes = plt.gca()
axes.set_xlim([0,10000000])
plt.show()







#THIS IS FOR TOTAL ENERGY GRAPHS
Euler = np.load("/Users/annabellecarver/Desktop/Common/1e5_4total_E_ENG.npy", allow_pickle=True)
EulerCromer = np.load("/Users/annabellecarver/Desktop/Common/1e5_4total_EC_ENG.npy", allow_pickle=True)
Verlet = np.load("/Users/annabellecarver/Desktop/Common/1e5_4total_V_ENG.npy", allow_pickle=True)


x1=[]
y1=[]

x2=[]
y2=[]

x3=[]
y3=[]

for i in range(len(Euler)):
    y1.append(Euler[i][1])
    x1.append(Euler[i][0])

for i in range(len(EulerCromer)):
    y2.append(EulerCromer[i][1])
    x2.append(EulerCromer[i][0])

for i in range(len(Verlet)):
    y3.append(Verlet[i][1])
    x3.append(Verlet[i][0])


plt.plot(x1, y1, label = "Euler")
plt.plot(x2, y2, label = "Euler-Cromer")
plt.plot(x3, y3, label = "Verlet")
#plt.title("Total Energy")
plt.xlabel("Time (s)")
plt.ylabel("Total Energy (J)")
plt.legend(loc = 'upper right')
axes = plt.gca()
axes.set_xlim([0,10000000])
plt.show()



#THIS IS FOR TOTAL MOMENTUM GRAPHS
Euler = np.load("/Users/annabellecarver/Desktop/Common/1e3_4total_E_MO.npy", allow_pickle=True)
EulerCromer = np.load("/Users/annabellecarver/Desktop/Common/1e3_4total_EC_MO.npy", allow_pickle=True)
Verlet = np.load("/Users/annabellecarver/Desktop/Common/1e3_4total_V_MO.npy", allow_pickle=True)


x1=[]
y1=[]

x2=[]
y2=[]

x3=[]
y3=[]

for i in range(len(Euler)):
    y1.append(Euler[i][1])
    x1.append(Euler[i][0])

for i in range(len(EulerCromer)):
    y2.append(EulerCromer[i][1])
    x2.append(EulerCromer[i][0])

for i in range(len(Verlet)):
    y3.append(Verlet[i][1])
    x3.append(Verlet[i][0])


plt.plot(x1, y1, label = "Euler")
plt.plot(x2, y2, label = "Euler-Cromer")
plt.plot(x3, y3, label = "Verlet")
#plt.title("Total Momentum")
plt.xlabel("Time (s)")
plt.ylabel("Total Momentum (kgm/s)")
plt.legend(loc = 'upper right')
axes = plt.gca()
axes.set_xlim([0,10000000])
plt.show()
