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

    def __init__(self):
        self.bodies = []
        self.min_kin = 1e500
        self.scale_factor = 0 #scale_factor
        self.counter = 0
        self.spacing = 500

    def read_file(self, read_bodies):
        try:
            # collect strings from database file
            import os
            basepath = os.path.dirname(os.path.realpath(__file__))
            path = ""
            for entry in os.listdir(basepath):
                if os.path.isfile(os.path.join(basepath, entry)):
                    if entry.endswith(".txt"):
                        path = basepath + "\\" + entry
            file = open(path, "r")
        except:
            print("Could not find database text file in current directory: " + str(basepath))
            raise SystemExit

        lines = file.readlines()
        file.close()

        # format and turn columns from database into celestial body objects
        for index, line in enumerate(lines):
            if len(line.strip()) > 0 and line.strip()[0] != "#": # ignore #s and empty lines
                column = line.split(",")
                if index == 1: # read constructor data from first line
                    self.iterations = int(column[0].strip())
                    self.step = float(column[1].strip())
                    self.G = float(column[2].strip()) * 1e-11
                    self.scale_factor = float(column[3].strip())
                    self.animation_step = int(column[4].strip())
                    self.limits = int(column[5].strip())
                elif read_bodies: # determine properties of body
                    name = column[0].strip().capitalize()
                    mass = float(column[1].strip()) * np.sqrt(self.scale_factor)
                    position = Vec2D(float(column[2].strip()), float(column[3].strip()))
                    velocity = Vec2D(0, 0) # initalized within update_body_properties function, easier to calculate there
                    color = "#" + column[4].strip()
                    size = float(column[5].strip())
                    self.bodies.append(CelestialBody(name, mass, position, velocity, color, size))


    def define_bodies(self):

        self.read_file(True)
        # determine the center (for use cases where the center doesn't move)
        self.center = "Mars"

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
        return 0.5 * body.get_mass() / np.sqrt(self.scale_factor) * abs(body.get_velocity()) ** 2

    def update_body_properties(self): # update and set all states during one iteration of the animation / calculation

        total_kinetic_energy = 0

        for body in self.bodies: # find all vel / accels
            pos = body.get_position()
            vel = Vec2D(0, 0)
            # set initial velocity
            if body.get_velocity().isZero():
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

        if self.counter % self.spacing == 0: #used for delaying the printing of total kinetic energy
            print("Total kinetic energy of system: " + str(total_kinetic_energy))

        self.counter += 1 # loop counter must increment every iteration

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

sim = Simulation()
sim.run()
