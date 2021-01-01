import numpy as np
import math
import copy
from Particle import Particle

class SolarSystem:
    """
    Utilising the methods from the particle file, the solar system class coalesces all the given n- bodies to simulate the solar system.
    the file then calculates the following positions and velocities as well as the total kinetic energy, total gravitational potential energy, total energy, total momentum and total angular momentum of the system as a whole.

    Parameters
    ----------
    the_planets: list
        A list of defined planets which are Particles from the Particle class
    time: float
        The time of the system

    Returns
    -------
    updated gravitational acceleration: array
        The gravitational acceleration at that point.
    update: array
        updates the positions which regard to which algorithm chosen.
    position_velocity: array
        Uses the update and updated gravitational acceleration to give new position velocities by first updating the acceleration and then using this to get the new positions and velocities.
    total angular momentum: float
        Sum of all the planets angular momentums.
    total momentum: float
        Sum of all the planets momentums.
    Gravitational potential energy: float
        Calculates the gravitational potential energy of all the planets.
    total kinetic energy: float
        Sum of all the planets kinetic energy.
    total energy: float
        Sum of all the planets kinetic and potential energy.

    """


    def __init__(self, planets):
        self.the_planets = np.array(planets)
        self.time=0
        self.UPPER_G = 6.67408E-11


    def __result__(self):
        return '{} at t={}'.format(self.the_planets,self.time)

    
    def update_gravitational_acceleration(self):
        UPPER_G = 6.67408E-11
        for m in self.the_planets:
            m.acceleration = np.array([0, 0, 0], dtype=float)
            for M in self.the_planets:
                if m==M:
                    continue
                radius_vector = M.position - m.position
                radius_distance = np.linalg.norm(radius_vector)
                g_mass_radius = UPPER_G * M.mass * radius_vector
                m.acceleration += g_mass_radius/(radius_distance**3)
        return m.acceleration
                

    def update(self, delta_t):
        for planet in self.the_planets:
            print(planet.algorithm)
            if planet.algorithm != "EULERCROMER":
                old_position = planet.position
                old_velocity = planet.velocity
                old_acceleration = planet.acceleration
                new_position = planet.position + (planet.velocity * delta_t)
                planet.position = new_position
                new_velocity = planet.velocity + (planet.acceleration * delta_t)
                planet.velocity = new_velocity

                if planet.algorithm == "VERLET":
                    self.update_gravitational_acceleration()
                    new_position = old_position + (old_velocity * delta_t) + (0.5 * old_acceleration * ((delta_t)**2))
                    planet.position = new_position
                    new_velocity = old_velocity + (0.5 * (planet.acceleration + old_acceleration) * delta_t)
                    planet.velocity = new_velocity

            elif planet.algorithm == "EULERCROMER":
                new_velocity = planet.velocity + (planet.acceleration * delta_t)
                planet.velocity = new_velocity
                new_position = planet.position + (planet.velocity * delta_t)
                planet.position = new_position
            else:
                print("No algorithm chosen.")

    def position_velocity(self, delta_t):
        self.update_gravitational_acceleration()
        self.update(delta_t)
        self.time += delta_t
    
    def total_momentum(self, delta_t):
        mom = 0 #initial momentum
        for p in self.the_planets:
            self.position_velocity(delta_t)
            mom += np.linalg.norm(p.momentum())
        return mom
    
    def total_gravitational_potential(self, delta_t):
        UPPER_G = 6.67408E-11
        for x in self.the_planets:
            x.potential = 0
            for y in self.the_planets:
                if x==y:
                    continue
                radius_vector = y.position - x.position
                radius_distance = np.linalg.norm(radius_vector)
                x.potential += (UPPER_G*y.mass*x.mass)/(radius_distance)
        potential_energy = 0
        for p in self.the_planets:
            self.position_velocity(delta_t)
            potential_energy += p.potential
        return potential_energy
  
    def total_kinetic_energy(self, delta_t):
        kinetic_en = 0
        for p in self.the_planets:
            self.position_velocity(delta_t)
            kinetic_en += p.kinetic_energy()
        return kinetic_en


    def total_energy(self, delta_t):
        energy = self.total_gravitational_potential(delta_t) + self.total_kinetic_energy(delta_t)
        return energy

    def total_angular_momentum(self, delta_t):
        total_a_m = 0
        for p in self.the_planets:
            self.position_velocity(delta_t)
            total_a_m += np.linalg.norm(p.angular_momentum())
        return total_a_m


        
    



