import math
import numpy as np
from scipy.spatial import distance


class Particle:
    """
    Defines a generic body for the solar system simulation. 
    In addition to containing all base variables required for manipulating the system’s planets.
    This class also calculates kinetic energy and momentum of a body at a given position

    Parameters
    ----------
    position: array
        The position of the particle
    velocity: array
        The velocity of the particle
    acceleration: array
        The acceleration of the particle
    name: string
        The name of the particle        
    mass: float
        The mass of the particle
    algorithm: float
        which approximation to use on the particle.
    
    Returns
    -------
    kinetic energy: float
        The kinetic energy of the particle
    momentum: array
        The momentum of the particle
    angular momentum: array
        The angular momentum of the particle
    """
    def __init__(self, position=np.array([0, 0, 0], dtype=float), velocity=np.array([0, 0, 0], dtype=float), acceleration=np.array([0, -10, 0], dtype=float), name='Ball', mass=1.0, UPPER_G = 6.67408E-11, algorithm= "EULER"):
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.acceleration = np.array(acceleration, dtype=float)
        self.name = name
        self.mass = mass
        self.algorithm = algorithm
        self.UPPER_G = UPPER_G

    def __result__(self):
        return "Particle: {0}, Mass: {1:.3e}, Position: {2}, Velocity: {3}, Acceleration: {4}".format(
            self.name, self.mass, self.position, self.velocity, self.acceleration)
    #changed the method name to result as this makes more sense to me.

    def kinetic_energy(self):
        kinetic_e = 0.5*self.mass*np.dot(self.velocity,self.velocity)
        return kinetic_e
    #I have used kinetic_e as the variable to differentiate from the method name.

    def momentum(self):
        m = (self.mass * self.velocity)
        return m
    #I have used m as the variable to differentiate from the method name.

    def angular_momentum(self):
        L = self.mass * self.velocity * self.position
        return L



    

