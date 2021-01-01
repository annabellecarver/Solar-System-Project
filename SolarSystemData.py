
import numpy as np
from Particle import Particle
from SolarSystem import SolarSystem
import copy
import scipy
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
from poliastro import constants
from spiceypy import sxform, mxvg
UPPER_G = 6.67408E-11

"""
Npy files storing all critical data for any programmed simulation is made in this file. Here, data from NASA’s HORIZONS interface is ran through the solar system’s methods to predict each particle’s journey through space and its associated energy[4]. By changing one key word, the time step and range of iteration, six different data files can be made for a certain algorithm over the chosen time frame.
"""


t = Time("2020-12-1 17:00:00.0", scale="tdb")
#to change the approximation used change the method of each planet to EULER, EULERCROMER OR VERLET



#for the sun
mass_sun = (constants.GM_sun / UPPER_G).value
pos_sun, vel_sun = get_body_barycentric_posvel("sun", t, ephemeris="jpl")
sunarray = [
    pos_sun.xyz[0].to("m").value,
    pos_sun.xyz[1].to("m").value,
    pos_sun.xyz[2].to("m").value,
    vel_sun.xyz[0].to("m/s").value,
    vel_sun.xyz[1].to("m/s").value,
    vel_sun.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
sunarrayecl = mxvg(trans, sunarray, 6, 6)
position_sun = [sunarrayecl[0], sunarrayecl[1], sunarrayecl[2]]
velocity_sun = [sunarrayecl[3], sunarrayecl[4], sunarrayecl[5]]
sun = Particle(position_sun, velocity_sun, np.array([0,0,0]), name='Sun', mass = mass_sun, algorithm = "EULER")
#print(mass_sun)
#print(velocity_sun)

#for mercury
mass_mercury = (constants.GM_mercury / UPPER_G).value
pos_mercury, vel_mercury = get_body_barycentric_posvel("mercury", t, ephemeris="jpl")
mercuryarray = [
    pos_mercury.xyz[0].to("m").value,
    pos_mercury.xyz[1].to("m").value,
    pos_mercury.xyz[2].to("m").value,
    vel_mercury.xyz[0].to("m/s").value,
    vel_mercury.xyz[1].to("m/s").value,
    vel_mercury.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
mercuryarrayecl = mxvg(trans, mercuryarray, 6, 6)
position_mercury = [mercuryarrayecl[0], mercuryarrayecl[1], mercuryarrayecl[2]]
velocity_mercury = [mercuryarrayecl[3], mercuryarrayecl[4], mercuryarrayecl[5]]
mercury = Particle(position_mercury, velocity_mercury, np.array([0,0,0]), name='Mercury', mass = mass_mercury, algorithm = "EULER")
#print(velocity_mercury)

#for venus
mass_venus = (constants.GM_venus / UPPER_G).value
pos_venus, vel_venus = get_body_barycentric_posvel("venus", t, ephemeris="jpl")
venusarray = [
    pos_venus.xyz[0].to("m").value,
    pos_venus.xyz[1].to("m").value,
    pos_venus.xyz[2].to("m").value,
    vel_venus.xyz[0].to("m/s").value,
    vel_venus.xyz[1].to("m/s").value,
    vel_venus.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
venusarrayecl = mxvg(trans, venusarray, 6, 6)
position_venus = [venusarrayecl[0], venusarrayecl[1], venusarrayecl[2]]
velocity_venus = [venusarrayecl[3], venusarrayecl[4], venusarrayecl[5]]
venus = Particle(position_venus, velocity_venus, np.array([0,0,0]), name='Venus', mass = mass_venus, algorithm = "EULER")
#print(velocity_venus)


#for the earth
mass_earth = (constants.GM_earth / UPPER_G).value
pos_earth, vel_earth = get_body_barycentric_posvel("earth", t, ephemeris="jpl")
eartharray = [
    pos_earth.xyz[0].to("m").value,
    pos_earth.xyz[1].to("m").value,
    pos_earth.xyz[2].to("m").value,
    vel_earth.xyz[0].to("m/s").value,
    vel_earth.xyz[1].to("m/s").value,
    vel_earth.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
eartharrayecl = mxvg(trans, eartharray, 6, 6)
position_earth = [eartharrayecl[0], eartharrayecl[1], eartharrayecl[2]]
velocity_earth = [eartharrayecl[3], eartharrayecl[4], eartharrayecl[5]]
earth = Particle(position_earth, velocity_earth, np.array([0,0,0]), name='Earth', mass = mass_earth, algorithm = "EULER")
#print(velocity_earth)

#for the earth's moon
mass_moon = (constants.GM_moon / UPPER_G).value
pos_moon, vel_moon = get_body_barycentric_posvel("moon", t, ephemeris="jpl")
moonarray = [
    pos_moon.xyz[0].to("m").value,
    pos_moon.xyz[1].to("m").value,
    pos_moon.xyz[2].to("m").value,
    vel_moon.xyz[0].to("m/s").value,
    vel_moon.xyz[1].to("m/s").value,
    vel_moon.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
moonarrayecl = mxvg(trans, moonarray, 6, 6)
position_moon = [moonarrayecl[0], moonarrayecl[1], moonarrayecl[2]]
velocity_moon = [moonarrayecl[3], moonarrayecl[4], moonarrayecl[5]]
moon = Particle(position_moon, velocity_moon, np.array([0,0,0]), name='Moon', mass = mass_moon, algorithm = "EULER")
#print(velocity_moon)

#for mars
mass_mars = (constants.GM_mars / UPPER_G).value
pos_mars, vel_mars = get_body_barycentric_posvel("mars", t, ephemeris="jpl")
marsarray = [
    pos_mars.xyz[0].to("m").value,
    pos_mars.xyz[1].to("m").value,
    pos_mars.xyz[2].to("m").value,
    vel_mars.xyz[0].to("m/s").value,
    vel_mars.xyz[1].to("m/s").value,
    vel_mars.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
marsarrayecl = mxvg(trans, marsarray, 6, 6)
position_mars = [marsarrayecl[0], marsarrayecl[1], marsarrayecl[2]]
velocity_mars = [marsarrayecl[3], marsarrayecl[4], marsarrayecl[5]]
mars = Particle(position_mars, velocity_mars, np.array([0,0,0]), name='Mars', mass = mass_mars, algorithm = "EULER")
#print(velocity_mars)

#for jupiter
mass_jupiter = (constants.GM_jupiter / UPPER_G).value
pos_jupiter, vel_jupiter = get_body_barycentric_posvel("jupiter", t, ephemeris="jpl")
jupiterarray = [
    pos_jupiter.xyz[0].to("m").value,
    pos_jupiter.xyz[1].to("m").value,
    pos_jupiter.xyz[2].to("m").value,
    vel_jupiter.xyz[0].to("m/s").value,
    vel_jupiter.xyz[1].to("m/s").value,
    vel_jupiter.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
jupiterarrayecl = mxvg(trans, jupiterarray, 6, 6)
position_jupiter = [jupiterarrayecl[0], jupiterarrayecl[1], jupiterarrayecl[2]]
velocity_jupiter = [jupiterarrayecl[3], jupiterarrayecl[4], jupiterarrayecl[5]]
jupiter = Particle(position_jupiter, velocity_jupiter, np.array([0,0,0]), name='Jupiter', mass = mass_jupiter, algorithm = "EULER")
#print(velocity_jupiter)

#for saturn
mass_saturn = (constants.GM_saturn / UPPER_G).value
pos_saturn, vel_saturn = get_body_barycentric_posvel("saturn", t, ephemeris="jpl")
saturnarray = [
    pos_saturn.xyz[0].to("m").value,
    pos_saturn.xyz[1].to("m").value,
    pos_saturn.xyz[2].to("m").value,
    vel_saturn.xyz[0].to("m/s").value,
    vel_saturn.xyz[1].to("m/s").value,
    vel_saturn.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
saturnarrayecl = mxvg(trans, saturnarray, 6, 6)
position_saturn = [saturnarrayecl[0], saturnarrayecl[1], saturnarrayecl[2]]
velocity_saturn = [saturnarrayecl[3], saturnarrayecl[4], saturnarrayecl[5]]
saturn = Particle(position_saturn, velocity_saturn, np.array([0,0,0]), name='Saturn', mass = mass_saturn, algorithm = "EULER")
#print(velocity_saturn)

#for uranus
mass_uranus = (constants.GM_uranus / UPPER_G).value
pos_uranus, vel_uranus = get_body_barycentric_posvel("uranus", t, ephemeris="jpl")
uranusarray = [
    pos_uranus.xyz[0].to("m").value,
    pos_uranus.xyz[1].to("m").value,
    pos_uranus.xyz[2].to("m").value,
    vel_uranus.xyz[0].to("m/s").value,
    vel_uranus.xyz[1].to("m/s").value,
    vel_uranus.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
uranusarrayecl = mxvg(trans, uranusarray, 6, 6)
position_uranus = [uranusarrayecl[0], uranusarrayecl[1], uranusarrayecl[2]]
velocity_uranus = [uranusarrayecl[3], uranusarrayecl[4], uranusarrayecl[5]]
uranus = Particle(position_uranus, velocity_uranus, np.array([0,0,0]), name='Uranus',mass = mass_uranus, algorithm = "EULER")
#print(velocity_uranus)

#for neptune
mass_neptune = (constants.GM_neptune / UPPER_G).value
pos_neptune, vel_neptune = get_body_barycentric_posvel("neptune", t, ephemeris="jpl")
neptunearray = [
    pos_neptune.xyz[0].to("m").value,
    pos_neptune.xyz[1].to("m").value,
    pos_neptune.xyz[2].to("m").value,
    vel_neptune.xyz[0].to("m/s").value,
    vel_neptune.xyz[1].to("m/s").value,
    vel_neptune.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
neptunearrayecl = mxvg(trans, neptunearray, 6, 6)
position_neptune = [neptunearrayecl[0], neptunearrayecl[1], neptunearrayecl[2]]
velocity_neptune = [neptunearrayecl[3], neptunearrayecl[4], neptunearrayecl[5]]
neptune = Particle(position_neptune, velocity_neptune, np.array([0,0,0]), name='Neptune', mass = mass_neptune, algorithm = "EULER")
#print(velocity_neptune)

#for pluto
mass_pluto = (constants.GM_pluto / UPPER_G).value
pos_pluto, vel_pluto = get_body_barycentric_posvel("pluto", t, ephemeris="jpl")
plutoarray = [
    pos_pluto.xyz[0].to("m").value,
    pos_pluto.xyz[1].to("m").value,
    pos_pluto.xyz[2].to("m").value,
    vel_pluto.xyz[0].to("m/s").value,
    vel_pluto.xyz[1].to("m/s").value,
    vel_pluto.xyz[2].to("m/s").value,
]
trans = sxform("J2000", "ECLIPJ2000", t.jd)
plutoarrayecl = mxvg(trans, plutoarray, 6, 6)
position_pluto = [plutoarrayecl[0], plutoarrayecl[1], plutoarrayecl[2]]
velocity_pluto = [plutoarrayecl[3], plutoarrayecl[4], plutoarrayecl[5]]
pluto = Particle(position_pluto, velocity_pluto, name='Pluto', mass = mass_pluto, algorithm = "EULER")
#print(velocity_pluto)

# note that the planets I can use are: sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune, pluto. to change the planets involved simply copy and paste more in!
the_planets = [sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune, pluto]
#to change the approximation used change the method of each planet to EULER, EULERCROMER OR VERLET


Data = []
delta_t = 10000
time = 0
#To get different time steps change delta_t to smaller or bigger and to change the range do that in each graph formation.
#To get different graph outputs then simply make the others into docstrings.



#THIS IS THE CODE FOR DOING THE SOLAR SYSTEM PLOTS
shortcut = SolarSystem(the_planets)
for n in range(0, 10000):
    shortcut.position_velocity(delta_t)

    if (n % 100==0):
        item = [shortcut.time, copy.deepcopy(shortcut)]
        Data.append(item)

np.save("/Users/annabellecarver/Desktop/Common/Unitmass", Data, allow_pickle=True)



#THIS IS FOR FINDING THE ANGULAR MOMENTUM GRAPHS
for i in range(0, 100000):
    AM = SolarSystem(the_planets).total_angular_momentum(delta_t)
    time += delta_t
    if i % 100==0:
        item = [time, copy.deepcopy(AM)]
        Data.append(item)
np.save("/Users/annabellecarver/Desktop/Common/1e44_EC_AM", Data, allow_pickle=True)


#THIS IS FOR FINDING THE KE GRAPHS
for i in range(0, 10000):
    KE = SolarSystem(the_planets).total_kinetic_energy(delta_t)
    time += delta_t
    if i % 100==0:
        item = [time, copy.deepcopy(KE)]
        Data.append(item)
np.save("/Users/annabellecarver/Desktop/Common/1e3_4total_V_KE", Data, allow_pickle=True)



#THIS IS FOR FINDING THE PE GRAPHS
for i in range(0, 10000):
    PE = SolarSystem(the_planets).total_gravitational_potential(delta_t)
    time += delta_t
    if i % 100==0:
        item = [time, copy.deepcopy(PE)]
        Data.append(item)
np.save("/Users/annabellecarver/Desktop/Common/1e3_4total_V_PE", Data, allow_pickle=True)


#THIS IS FOR FINDING THE TOTAL ENERGY GRAPHS
for i in range(0, 10000):
    ENG = SolarSystem(the_planets).total_energy(delta_t)
    time += delta_t
    if i % 100==0:
        item = [time, copy.deepcopy(ENG)]
        Data.append(item)
np.save("/Users/annabellecarver/Desktop/Common/1e3_4total_V_ENG", Data, allow_pickle=True)



#THIS IS FOR FINDING THE TOTAL MOMENTUM GRAPHS
for i in range(0, 100000):
    MO = SolarSystem(the_planets).total_momentum(delta_t)
    time += delta_t
    if i % 100==0:
        item = [time, copy.deepcopy(MO)]
        Data.append(item)
np.save("/Users/annabellecarver/Desktop/Common/1e5_5total_EC_MO", Data, allow_pickle=True)


