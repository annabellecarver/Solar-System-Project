import numpy as np
from Particle import Particle
from SolarSystem import SolarSystem
import copy
import scipy
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
from poliastro import constants
from spiceypy import sxform, mxvg
from matplotlib import pyplot as plt
from os import path
UPPER_G = 6.67408E-11


"""
This file is the testing file where the gravitational acceleration on the surface of the planets is tested and maximum and minimum of a graph can be determined.
"""

mass_massless = 1
position_massless = [3389600, 0, 0] #just above the radius of the earth or mars
velocity_massless = [0, 0, 0]
massless = Particle(position_massless, velocity_massless, np.array([0,0,0]), name='Massless', mass = mass_massless, algorithm = "EULER")
#we should see that this particle is not affected by the masses of the solar system 

mass_earthtest = (constants.GM_earth / UPPER_G).value
position_earthtest = [0, 0, 0]
velocity_earthtest = [0, 0, 0]
earthtest = Particle(position_earthtest, velocity_earthtest, np.array([0,0,0]), name='Earthtest', mass = mass_earthtest, algorithm = "EULER")

mass_mars = (constants.GM_mars / UPPER_G).value
position_mars = [0, 0, 0]
velocity_mars = [0, 0, 0]
mars = Particle(position_mars, velocity_mars, np.array([0,0,0]), name='Mars', mass = mass_mars, algorithm = "EULER")

the_planets = [mars, massless]
#test for gravitational acceleration on surface of a planet,


acceleration = SolarSystem(the_planets).update_gravitational_acceleration()
print(acceleration)






#this part is for plotting 3D plots to test the maximum and minimum

Data = np.load("/Users/annabellecarver/Desktop/Common/MasslessTest3D.npy", allow_pickle=True)
Data_T = Data.T

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
from mpl_toolkits.mplot3d import Axes3D

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

maximum = max(x_i)
print(maximum)
minimum = min(x_i)
print(minimum)





#this section is for plotting 2D graphs in order to get maximum and minimums.
"""


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

maximum = max(x)
print(maximum)
minimum = min(x)
print(minimum)

"""