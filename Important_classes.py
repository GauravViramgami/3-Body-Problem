from Global_constants import *
import matplotlib.pyplot as plt

class Object2D:
    def __init__ (self, name, mass, position_initial, velocity_initial):
        self.name = name
        self.mass = mass
        self.position_initial = position_initial
        self.velocity_initial = velocity_initial
        self.position = position_initial
        self.velocity = velocity_initial
        self.acceleration = 0
        
    def euclideanDistance (self, x, y):
        r = (((self.position[0] - x)**2) + ((self.position[1] - y)**2))**(1/2)
        return r
    
    def updatePosition (self, x_new, y_new):
        self.position[0] = x_new
        self.position[1] = y_new
        
    def updateVelocity (self, vx_new, vy_new):
        self.velocity[0] = vx_new
        self.velocity[1] = vy_new

class TwoBodySystem:
    def __init__ (self, reference_object2D, orbiting_object2D):
        self.reference_object2D = reference_object2D
        self.orbiting_object2D = orbiting_object2D
        self.initialize_objects()
        self.euler_trajectory = [[], []]
        self.verlet_trajectory = [[], []]

    def initialize_objects(self):
        self.reference_object2D.position = self.reference_object2D.position_initial
        self.orbiting_object2D.position = self.orbiting_object2D.position_initial
        self.reference_object2D.velocity = self.reference_object2D.velocity_initial
        self.orbiting_object2D.velocity = self.orbiting_object2D.velocity_initial
        self.updateAcceleration()

    def updateAcceleration(self):
        self.orbiting_object2D.acceleration = [-1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * self.orbiting_object2D.position[0]) / ((self.orbiting_object2D.euclideanDistance(0, 0))**3), -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * self.orbiting_object2D.position[1]) / ((self.orbiting_object2D.euclideanDistance(0, 0))**3)]

    def euler_cromer_method(self, stepsize, num_iterations):
        self.initialize_objects()
        self.euler_trajectory = [[self.orbiting_object2D.position_initial[0]], [self.orbiting_object2D.position_initial[1]]]
        
        x_prev, y_prev = self.orbiting_object2D.position_initial
        vx_prev, vy_prev = self.orbiting_object2D.velocity_initial

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            r = self.orbiting_object2D.euclideanDistance(0, 0)

            self.orbiting_object2D.updateVelocity(vx_prev - ((CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * x_prev)/(r**3)) * stepsize, vy_prev - ((CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * y_prev)/(r**3)) * stepsize) 

            self.orbiting_object2D.updatePosition(x_prev + (self.orbiting_object2D.velocity[0] * stepsize), y_prev + (self.orbiting_object2D.velocity[1] * stepsize)) 
            # self.orbiting_object2D.updatePosition(x_prev + (vx_prev * stepsize), y_prev + (vy_prev * stepsize)) # It is working here as non-energy preserving euler method
            self.euler_trajectory[0].append(self.orbiting_object2D.position[0])
            self.euler_trajectory[1].append(self.orbiting_object2D.position[1])

            x_prev, y_prev = self.orbiting_object2D.position
            vx_prev, vy_prev = self.orbiting_object2D.velocity

    def verlet_method (self, stepsize, num_iterations):
        self.initialize_objects()
        self.verlet_trajectory = [[self.orbiting_object2D.position_initial[0]], [self.orbiting_object2D.position_initial[1]]]

        x_prev, y_prev = self.orbiting_object2D.position_initial
        vx_prev, vy_prev = self.orbiting_object2D.velocity_initial
        ax_prev, ay_prev = self.orbiting_object2D.acceleration

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            r = self.orbiting_object2D.euclideanDistance(0, 0)

            self.orbiting_object2D.updatePosition(x_prev + (vx_prev * stepsize) + ((ax_prev * (stepsize**2))/2), y_prev + (vy_prev * stepsize) + ((ay_prev * (stepsize**2))/2))
            self.updateAcceleration()
            self.orbiting_object2D.updateVelocity(vx_prev + (((ax_prev + self.orbiting_object2D.acceleration[0]) * stepsize)/2), vy_prev + (((ay_prev + self.orbiting_object2D.acceleration[1]) * stepsize)/2))

            self.verlet_trajectory[0].append(self.orbiting_object2D.position[0])
            self.verlet_trajectory[1].append(self.orbiting_object2D.position[1])

            x_prev, y_prev = self.orbiting_object2D.position
            vx_prev, vy_prev = self.orbiting_object2D.velocity
            ax_prev, ay_prev = self.orbiting_object2D.acceleration

    def plot_euler_trajectory(self):
        plt.plot(self.euler_trajectory[0],self.euler_trajectory[1])
        plt.show()

    def plot_verlet_trajectory(self):
        plt.plot(self.verlet_trajectory[0],self.verlet_trajectory[1])
        plt.show()