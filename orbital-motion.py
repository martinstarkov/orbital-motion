import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

class Vec2D(object):

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def __add__(self, vec2):
        return Vec2D(self.x + vec2.x, self.y + vec2.y)

    def __sub__(self, vec2):
        return Vec2D(self.x - vec2.x, self.y - vec2.y)

    def __mul__(self, vec2):
        return self.x * vec2.x + self.y * vec2.y

    def magnitude(self):
        return math.sqrt(self.x ** 2 + self.y ** 2)

    def multiply(self, constant):
        return Vec2D(self.x * constant, self.y * constant)

    def unit(self):
        if self.__abs__() != 0:
            return Vec2D(self.x / self.__abs__(), self.y / self.__abs__())
        else:
            return Vec2D(self.x, self.y)

    def __abs__(self):
        return np.sqrt(self.x ** 2 + self.y ** 2)

    def __eq__(self, vec2):
        return self.x == vec2.x and self.y == vec2.y

    def __str__(self):
        return '(%g, %g)' % (self.x, self.y)

    def __ne__(self, other):
        return not self.__eq__(other)

    def isZero(self):
        return self.x == 0 and self.y == 0

class CelestialBody(object):

    def __init__(self, name, mass, position, velocity, planet_color, size):
        self.name = name
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.patch = plt.Circle((position.get_x(), position.get_y()), size, color=planet_color, animated=False)

    def get_name(self):
        return self.name

    def get_mass(self):
        return self.mass

    def set_next_accel(self, accel):
        self.next_accel = accel

    def set_next_vel(self, vel):
        self.next_vel = vel

    def set_next_pos(self, pos):
        self.next_pos = pos

    def get_next_accel(self):
        return self.next_accel

    def get_next_vel(self):
        return self.next_vel

    def get_next_pos(self):
        return self.next_pos

    def get_patch(self):
        return self.patch

    def update_patch_pos(self, new_pos):
        self.patch.center = (new_pos.get_x(), new_pos.get_y())

    def get_orbital_radius(self, external_body):
        return self.position - external_body.get_position()

    def get_position(self):
        return self.position

    def set_position(self, new_pos):
        self.position = new_pos

    def get_velocity(self):
        return self.velocity

    def set_velocity(self, new_vel):
        self.velocity = new_vel

class Simulation(object):

    def __init__(self, iterations, step, G, scale_factor, animation_step, limits):
        self.limits = limits
        self.animation_step = animation_step
        self.scale_factor = scale_factor
        self.iterations = iterations
        self.step = step
        self.G = G
        self.bodies = []
        self.min_kin = 1e500

    def define_bodies(self):
        # # eventually read from file, but for now hardcoded
        self.center = "Mars"
        self.bodies.append(CelestialBody(self.center, 6.4185e23 * self.scale_factor ** 9, Vec2D(0, 0), Vec2D(0, 0), 'r', 3389500 / self.scale_factor * 5))
        self.bodies.append(CelestialBody("Phobos", 1.06e16 * self.scale_factor ** 9, Vec2D(9.3773e6, 0), Vec2D(0, 0), 'b', 11267 / self.scale_factor * 150))
        self.bodies.append(CelestialBody("Deimos", 1.80e15 * self.scale_factor ** 9, Vec2D(23.463e6, 0), Vec2D(0, 0), 'g', 6200 / self.scale_factor * 150))
        
        # self.center = "Sun"
        # self.bodies.append(CelestialBody(self.center, 10000 * 1e7 * self.scale_factor, Vec2D(0, 0), Vec2D(0, 0), '#FDB813', 10))
        # self.bodies.append(CelestialBody("Mercury", 1 * 1e7, Vec2D(15, 0), Vec2D(0, 0), '#68696d', 1))
        # self.bodies.append(CelestialBody("Venus", 100 * 1e7, Vec2D(20, 0), Vec2D(0, 0), '#8B7D82', 2))
        # self.bodies.append(CelestialBody("Earth", 123 * 1e7, Vec2D(30, 0), Vec2D(0, 0), '#9fc164', 3))
        # self.bodies.append(CelestialBody("Mars", 1 * 1e7, Vec2D(45, 0), Vec2D(0, 0), '#c1440e', 2.5))
        # self.bodies.append(CelestialBody("Jupiter", 500 * 1e7, Vec2D(80, 0), Vec2D(0, 0), '#c99039', 7))
        # self.bodies.append(CelestialBody("Neptune", 200 * 1e7, Vec2D(105, 0), Vec2D(0, 0), '#3f54ba', 4))
        # self.bodies.append(CelestialBody("Saturn", 400 * 1e7, Vec2D(125, 0), Vec2D(0, 0), '#D8B763', 6))
        # self.bodies.append(CelestialBody("Uranus", 150 * 1e7, Vec2D(140, 0), Vec2D(0, 0), '#85addb', 4))

        pass

    def integrate(self, current, factor): # numerical intergration using Euler-Cromer
        return current + factor.multiply(self.step)

    def find_body(self, name): # find body by name from bodies array, very useful
        for body in self.bodies:
            if body.get_name() == name:
                return body
        return False

    def find_accel(self, body, external_body):
        orbital_radius = body.get_orbital_radius(external_body)
        abs_radius = abs(orbital_radius) ** 3
        acceleration = orbital_radius.multiply(-self.G * external_body.get_mass() / abs_radius) # divide force by mass
        return acceleration

    def find_vel(self, body, external_body):
        radius = body.get_orbital_radius(external_body)
        vel = Vec2D(0, math.sqrt(self.G * external_body.get_mass() / radius.magnitude())) # velocity equation
        return vel

    def net_property(self, body, property_type): # find net acceleration or velocity of a body
        net = Vec2D(0, 0) # initial
        for external_body in self.bodies: # cycle through all other bodies
            if body.get_name() != external_body.get_name(): # don't calculate accel/vel caused due to itself
                if property_type == "accel": # net acceleration
                    net += self.find_accel(body, external_body) 
                elif property_type == "vel": # net velocity
                    net += self.find_vel(body, external_body)
        return net

    def kinetic_energy(self, body): # simple kinetic energy calculation for a body
        return 0.5 * body.get_mass() * abs(body.get_velocity()) ** 2

    def update_body_properties(self): # update and set all states during one iteration of the animation / calculation

        total_kinetic_energy = 0

        for body in self.bodies: # find all vel / accels

            pos = body.get_position()
            vel = Vec2D(0, 0)

            # set initial velocity
            if body.get_velocity().isZero(): 
                # if requested that velocity only takes into account central object: 
                # if body.get_name() != self.center:
                #     vel = self.find_vel(body, self.find_body(self.center))
                # otherwise, take into account all objects: 
                vel = self.net_property(body, "vel")
            else:
                vel = body.get_velocity()

            accel = self.net_property(body, "accel") # net acceleration

            next_vel = self.integrate(vel, accel) # velocity numerical integration

            # keep central body still (optional, can be commented out)
            # if body.get_name() == self.center: 
            #     next_vel = Vec2D(0, 0)

            next_pos = self.integrate(pos, next_vel) # position numerical integration

            # store future position / velocity in body object for update loop
            body.set_next_pos(next_pos)
            body.set_next_vel(next_vel)
            body.set_next_accel(accel)
        
        patches = []

        for body in self.bodies: # update all vel / accels in one loop cycle

            # update body object
            body.update_patch_pos(body.get_next_pos())
            body.set_velocity(body.get_next_vel())
            body.set_position(body.get_next_pos())

            # tell animation to animate patch objects
            patches.append(body.get_patch())

            # calculate total kinetic energy
            total_kinetic_energy += self.kinetic_energy(body)

        # calculate miniumum and maximum kinetic energies over many cycles
        if total_kinetic_energy < self.min_kin:
            self.min_kin = total_kinetic_energy

        #print("Minimum kinetic energy of system: " + str(self.min_kin))

        print("Total kinetic energy of system: " + str(total_kinetic_energy))

        return patches # iterable patch object for animate function, refreshed every iteration

    def iterate(self):
        for i in range(self.iterations):
            self.update_body_properties()

    def animate(self, i):
        return self.update_body_properties()

    def display(self):

        # create matplotlib elements needed for graphing
        self.fig = plt.figure()
        self.ax = plt.axes()

        # add patches to axes
        for body in self.bodies:
            self.ax.add_patch(body.get_patch())

        # set axes scaling and limits
        self.ax.axis('scaled')
        self.ax.set_xlim(-self.limits / 2, self.limits / 2)
        self.ax.set_ylim(-self.limits / 2, self.limits / 2)

        # animator object
        anim = FuncAnimation(self.fig, self.animate, self.iterations, repeat=True, interval=self.animation_step, blit=True)

        plt.show()

    def run(self):

        # define initial conditions for bodies
        self.define_bodies()

        # iterate and animate objects
        self.display()

        # iterate objects
        #self.iterate()

# limits: 23.463e6 * 3 for Mars, 300 for Solar System
# scales: 5e0 for Mars, 1e3 for Solar System

#  iterations, step, G, scale_factor, animation_step, limits
sim = Simulation(1000, 0.1, 6.67408e-11, 5e0, 1, 23.463e6 * 3)
sim.run()