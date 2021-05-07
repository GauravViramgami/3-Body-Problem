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
        
    def euclideanDistance (self, from_position):
        r = (((self.position[0] - from_position[0])**2) + ((self.position[1] - from_position[1])**2))**(1/2)
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
        acceleration_reference_orbiting_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (self.orbiting_object2D.position[0] - self.reference_object2D.position[0])) / ((self.orbiting_object2D.euclideanDistance(self.reference_object2D.position))**3)
        acceleration_reference_orbiting_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (self.orbiting_object2D.position[1] - self.reference_object2D.position[1])) / ((self.orbiting_object2D.euclideanDistance(self.reference_object2D.position))**3)
        self.orbiting_object2D.acceleration = [acceleration_reference_orbiting_x, acceleration_reference_orbiting_y]

    def euler_cromer_method(self, stepsize, num_iterations):
        self.initialize_objects()
        self.euler_trajectory = [[self.orbiting_object2D.position_initial[0]], [self.orbiting_object2D.position_initial[1]]]
        
        x_prev, y_prev = self.orbiting_object2D.position_initial
        vx_prev, vy_prev = self.orbiting_object2D.velocity_initial

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            r = self.orbiting_object2D.euclideanDistance(self.reference_object2D.position)

            self.orbiting_object2D.updateVelocity(vx_prev + (self.orbiting_object2D.acceleration[0] * stepsize), vy_prev + (self.orbiting_object2D.acceleration[1] * stepsize)) 

            self.orbiting_object2D.updatePosition(x_prev + (self.orbiting_object2D.velocity[0] * stepsize), y_prev + (self.orbiting_object2D.velocity[1] * stepsize)) 

            # It is working here as non-energy preserving euler method
            # self.orbiting_object2D.updatePosition(x_prev + (vx_prev * stepsize), y_prev + (vy_prev * stepsize)) 

            self.updateAcceleration()
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
            r = self.orbiting_object2D.euclideanDistance(self.reference_object2D.position)

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

class ThreeBodySystem:
    def __init__ (self, reference_object2D, orbiting_object2D_1, orbiting_object2D_2):
        self.reference_object2D = reference_object2D
        self.orbiting_object2D_1 = orbiting_object2D_1
        self.orbiting_object2D_2 = orbiting_object2D_2
        self.initialize_objects()
        self.euler_trajectory = [[[], []], [[], []]]
        self.verlet_trajectory = [[[], []], [[], []]]

    def initialize_objects(self):
        self.reference_object2D.position = self.reference_object2D.position_initial
        self.orbiting_object2D_1.position = self.orbiting_object2D_1.position_initial
        self.orbiting_object2D_2.position = self.orbiting_object2D_2.position_initial
        self.reference_object2D.velocity = self.reference_object2D.velocity_initial
        self.orbiting_object2D_1.velocity = self.orbiting_object2D_1.velocity_initial
        self.orbiting_object2D_2.velocity = self.orbiting_object2D_2.velocity_initial
        self.updateAcceleration()

    def updateAcceleration(self):
        acceleration_reference_orbiting1_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (self.orbiting_object2D_1.position[0] - self.reference_object2D.position[0])) / ((self.orbiting_object2D_1.euclideanDistance(self.reference_object2D.position))**3)
        acceleration_orbiting2_orbiting1_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_2.name] * (self.orbiting_object2D_1.position[0] - self.orbiting_object2D_2.position[0])) / ((self.orbiting_object2D_1.euclideanDistance(self.orbiting_object2D_2.position))**3)
        acceleration_reference_orbiting1_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (self.orbiting_object2D_1.position[1] - self.reference_object2D.position[1])) / ((self.orbiting_object2D_1.euclideanDistance(self.reference_object2D.position))**3)
        acceleration_orbiting2_orbiting1_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_2.name] * (self.orbiting_object2D_1.position[1] - self.orbiting_object2D_2.position[1])) / ((self.orbiting_object2D_1.euclideanDistance(self.orbiting_object2D_2.position))**3)

        acceleration_reference_orbiting2_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (self.orbiting_object2D_2.position[0] - self.reference_object2D.position[0])) / ((self.orbiting_object2D_2.euclideanDistance(self.reference_object2D.position))**3)
        acceleration_orbiting1_orbiting2_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_1.name] * (self.orbiting_object2D_2.position[0] - self.orbiting_object2D_1.position[0])) / ((self.orbiting_object2D_2.euclideanDistance(self.orbiting_object2D_1.position))**3)
        acceleration_reference_orbiting2_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (self.orbiting_object2D_2.position[1] - self.reference_object2D.position[1])) / ((self.orbiting_object2D_2.euclideanDistance(self.reference_object2D.position))**3)
        acceleration_orbiting1_orbiting2_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_1.name] * (self.orbiting_object2D_2.position[1] - self.orbiting_object2D_1.position[1])) / ((self.orbiting_object2D_2.euclideanDistance(self.orbiting_object2D_1.position))**3)
        
        self.orbiting_object2D_1.acceleration = [acceleration_reference_orbiting1_x + acceleration_orbiting2_orbiting1_x, acceleration_reference_orbiting1_y + acceleration_orbiting2_orbiting1_y]
        self.orbiting_object2D_2.acceleration = [acceleration_reference_orbiting2_x + acceleration_orbiting1_orbiting2_x, acceleration_reference_orbiting2_y + acceleration_orbiting1_orbiting2_y]


    def euler_cromer_method(self, stepsize, num_iterations):
        self.initialize_objects()
        self.euler_trajectory = [[[self.orbiting_object2D_1.position_initial[0]], [self.orbiting_object2D_1.position_initial[1]]], [[self.orbiting_object2D_2.position_initial[0]], [self.orbiting_object2D_2.position_initial[1]]]]
        
        orbiting1_x_prev, orbiting1_y_prev = self.orbiting_object2D_1.position_initial
        orbiting2_x_prev, orbiting2_y_prev = self.orbiting_object2D_2.position_initial
        orbiting1_vx_prev, orbiting1_vy_prev = self.orbiting_object2D_1.velocity_initial
        orbiting2_vx_prev, orbiting2_vy_prev = self.orbiting_object2D_2.velocity_initial

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            r_reference_orbiting1 = self.orbiting_object2D_1.euclideanDistance(self.reference_object2D.position)
            r_reference_orbiting2 = self.orbiting_object2D_2.euclideanDistance(self.reference_object2D.position)
            r_orbiting1_orbiting2 = self.orbiting_object2D_1.euclideanDistance(self.orbiting_object2D_2.position)

            self.orbiting_object2D_1.updateVelocity(orbiting1_vx_prev + (self.orbiting_object2D_1.acceleration[0] * stepsize), orbiting1_vy_prev + (self.orbiting_object2D_1.acceleration[1] * stepsize)) 
            self.orbiting_object2D_2.updateVelocity(orbiting2_vx_prev + (self.orbiting_object2D_2.acceleration[0] * stepsize), orbiting2_vy_prev + (self.orbiting_object2D_2.acceleration[1] * stepsize)) 

            self.orbiting_object2D_1.updatePosition(orbiting1_x_prev + (self.orbiting_object2D_1.velocity[0] * stepsize), orbiting1_y_prev + (self.orbiting_object2D_1.velocity[1] * stepsize))
            self.orbiting_object2D_2.updatePosition(orbiting2_x_prev + (self.orbiting_object2D_2.velocity[0] * stepsize), orbiting2_y_prev + (self.orbiting_object2D_2.velocity[1] * stepsize)) 

            # It is working here as non-energy preserving euler method
            # self.orbiting_object2D_1.updatePosition(orbiting1_x_prev + (orbiting1_vx_prev * stepsize), orbiting1_y_prev + (orbiting1_vy_prev * stepsize))
            # self.orbiting_object2D_2.updatePosition(orbiting2_x_prev + (orbiting2_vx_prev * stepsize), orbiting2_y_prev + (orbiting2_vy_prev * stepsize)) 
            
            self.updateAcceleration()

            self.euler_trajectory[0][0].append(self.orbiting_object2D_1.position[0])
            self.euler_trajectory[0][1].append(self.orbiting_object2D_1.position[1])
            self.euler_trajectory[1][0].append(self.orbiting_object2D_2.position[0])
            self.euler_trajectory[1][1].append(self.orbiting_object2D_2.position[1])

            orbiting1_x_prev, orbiting1_y_prev = self.orbiting_object2D_1.position
            orbiting2_x_prev, orbiting2_y_prev = self.orbiting_object2D_2.position
            orbiting1_vx_prev, orbiting1_vy_prev = self.orbiting_object2D_1.velocity
            orbiting2_vx_prev, orbiting2_vy_prev = self.orbiting_object2D_2.velocity

    def verlet_method (self, stepsize, num_iterations):
        self.initialize_objects()
        self.verlet_trajectory = [[[self.orbiting_object2D_1.position_initial[0]], [self.orbiting_object2D_1.position_initial[1]]], [[self.orbiting_object2D_2.position_initial[0]], [self.orbiting_object2D_2.position_initial[1]]]]

        orbiting1_x_prev, orbiting1_y_prev = self.orbiting_object2D_1.position_initial
        orbiting2_x_prev, orbiting2_y_prev = self.orbiting_object2D_2.position_initial
        orbiting1_vx_prev, orbiting1_vy_prev = self.orbiting_object2D_1.velocity_initial
        orbiting2_vx_prev, orbiting2_vy_prev = self.orbiting_object2D_2.velocity_initial
        orbiting1_ax_prev, orbiting1_ay_prev = self.orbiting_object2D_1.acceleration
        orbiting2_ax_prev, orbiting2_ay_prev = self.orbiting_object2D_2.acceleration

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            r_reference_orbiting1 = self.orbiting_object2D_1.euclideanDistance(self.reference_object2D.position)
            r_reference_orbiting2 = self.orbiting_object2D_2.euclideanDistance(self.reference_object2D.position)
            r_orbiting1_orbiting2 = self.orbiting_object2D_1.euclideanDistance(self.orbiting_object2D_2.position)

            self.orbiting_object2D_1.updatePosition(orbiting1_x_prev + (orbiting1_vx_prev * stepsize) + ((orbiting1_ax_prev * (stepsize**2))/2), orbiting1_y_prev + (orbiting1_vy_prev * stepsize) + ((orbiting1_ay_prev * (stepsize**2))/2))
            self.orbiting_object2D_2.updatePosition(orbiting2_x_prev + (orbiting2_vx_prev * stepsize) + ((orbiting2_ax_prev * (stepsize**2))/2), orbiting2_y_prev + (orbiting2_vy_prev * stepsize) + ((orbiting2_ay_prev * (stepsize**2))/2))
            self.updateAcceleration()
            self.orbiting_object2D_1.updateVelocity(orbiting1_vx_prev + (((orbiting1_ax_prev + self.orbiting_object2D_1.acceleration[0]) * stepsize)/2), orbiting1_vy_prev + (((orbiting1_ay_prev + self.orbiting_object2D_1.acceleration[1]) * stepsize)/2))
            self.orbiting_object2D_2.updateVelocity(orbiting2_vx_prev + (((orbiting2_ax_prev + self.orbiting_object2D_2.acceleration[0]) * stepsize)/2), orbiting2_vy_prev + (((orbiting2_ay_prev + self.orbiting_object2D_2.acceleration[1]) * stepsize)/2))

            self.verlet_trajectory[0][0].append(self.orbiting_object2D_1.position[0])
            self.verlet_trajectory[0][1].append(self.orbiting_object2D_1.position[1])
            self.verlet_trajectory[1][0].append(self.orbiting_object2D_2.position[0])
            self.verlet_trajectory[1][1].append(self.orbiting_object2D_2.position[1])

            orbiting1_x_prev, orbiting1_y_prev = self.orbiting_object2D_1.position
            orbiting2_x_prev, orbiting2_y_prev = self.orbiting_object2D_2.position
            orbiting1_vx_prev, orbiting1_vy_prev = self.orbiting_object2D_1.velocity
            orbiting2_vx_prev, orbiting2_vy_prev = self.orbiting_object2D_2.velocity
            orbiting1_ax_prev, orbiting1_ay_prev = self.orbiting_object2D_1.acceleration
            orbiting2_ax_prev, orbiting2_ay_prev = self.orbiting_object2D_2.acceleration

    def plot_euler_trajectory(self):
        plt.plot(self.euler_trajectory[0][0], self.euler_trajectory[0][1], color = 'r', label = self.orbiting_object2D_1.name)
        plt.plot(self.euler_trajectory[1][0], self.euler_trajectory[1][1], color = 'g', label = self.orbiting_object2D_2.name)
        plt.show()

    def plot_verlet_trajectory(self):
        plt.plot(self.verlet_trajectory[0][0], self.verlet_trajectory[0][1], color = 'r', label = self.orbiting_object2D_1.name)
        plt.plot(self.verlet_trajectory[1][0], self.verlet_trajectory[1][1], color = 'g', label = self.orbiting_object2D_2.name)
        plt.show()