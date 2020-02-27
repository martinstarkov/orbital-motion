import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

from datetime import datetime
from matplotlib.animation import FuncAnimation
from random import randrange

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
        return math.sqrt(self.x**2 + self.y**2)

    def __eq__(self, vec2):
        return self.x == vec2.x and self.y == vec2.y

    def __str__(self):
        return '(%g, %g)' % (self.x, self.y)

    def __ne__(self, other):
        return not self.__eq__(other)

class CelestialBody(object):

    def __init__(self, name, mass, position, velocity):
        self.name = name
        self.mass = mass
        self.position = position
        self.velocity = velocity

    def get_name(self):
        return self.name

    def get_mass(self):
        return self.mass

    def get_orbital_radius(self, external_body):
        return external_body.get_position() - self.position

    def get_position(self):
        return self.position

    def set_position(self, new_pos):
        self.position = new_pos

    def get_velocity(self):
        return self.velocity

    def set_velocity(self, new_vel):
        self.velocity = new_vel
    
    def set_patch(self, new_patch):
        self.patch = new_patch

    def get_patch(self):
        return self.patch

class Simulation(object):

    def __init__(self, iterations, step, G):
        self.iterations = iterations
        self.step = step
        self.G = G
        self.bodies = []

    def integrate(self, current, factor):
        return current + factor.multiply(self.step)

    def define_bodies(self):
        # eventually read from file, but for now hardcoded
        self.bodies.append(CelestialBody("Mars", 6.4185e23, Vec2D(0, 0), Vec2D(0, 0)))
        self.bodies.append(CelestialBody("Phobos", 1.06e16, Vec2D(9.3773e6, 0), Vec2D(0, 0)))
        # self.bodies.append(CelestialBody("Deimos", 1.8e15, 23.463e6, Vec2D(23.463e6, 0), Vec2D(0, math.sqrt(self.G * 6.4185e23 / 23.463e6))))

    def find_body(self, name):
        for body in self.bodies:
            if body.get_name() == name:
                return body
        return False

    def find_acceleration(self, body, external_body):
        radius = body.get_orbital_radius(external_body) # orbital radius
        force = radius.unit().multiply(self.G * body.get_mass() * external_body.get_mass() / abs(radius) ** 2) # gravitation equation
        acceleration = force.multiply(1 / body.get_mass()) # divide force by mass
        return acceleration

    def find_velocity(self, body, external_body):
        radius = Vec2D(0, 0) - body.get_orbital_radius(external_body) # we want the opposite direction for velocity
        velocity = Vec2D(0, math.sqrt(self.G * external_body.get_mass() / radius.magnitude())) # velocity equation
        return velocity

    def net_property(self, body, motion_property):
        net = Vec2D(0, 0) # initial
        for external_body in self.bodies:
            if external_body.get_name() != body.get_name(): # don't calculate property caused due to itself
                if motion_property == "velocity":
                    net += self.find_velocity(body, external_body) # add net velocity to total
                elif motion_property == "acceleration":
                    net += self.find_acceleration(body, external_body) # add net acceleration to total
        return net


    def update_body_properties(self, iteration):
        for body in self.bodies:
            pos = body.get_position()
            vel = self.net_property(body, "velocity")
            accel = self.net_property(body, "acceleration")
            next_vel = self.integrate(vel, accel)
            next_pos = self.integrate(pos, next_vel)
            body.set_velocity(next_vel)
            body.set_position(next_pos) 

        
    def iterate(self):
        for i in range(self.iterations):
            self.update_body_properties(i)

    def update(self, frame):
        self.x_data.append(datetime.now())
        self.y_data.append(randrange(0, 100))
        self.line.set_data(self.x_data, self.y_data)
        self.figure.gca().relim()
        self.figure.gca().autoscale_view()
        return self.line,

    def run(self):

        # define initial conditions for bodies
        self.define_bodies()

        self.x_data, self.y_data = [], []

        self.figure = plt.figure()
        self.line, = plt.plot_date(self.x_data, self.y_data, '-')

        animation = FuncAnimation(self.figure, self.update, interval=200)

        plt.show()

        # run the amount of animations desired
        self.iterate()


sim = Simulation(100000, 0.1, 6.67408e-11)
sim.run()

