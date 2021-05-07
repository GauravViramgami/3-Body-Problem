from Global_constants import *
import matplotlib.pyplot as plt
import pygame

class Object2D:
    def __init__ (self, name, mass, position_initial, velocity_initial):
        self.name = name
        self.mass = mass
        self.position_initial = position_initial
        self.velocity_initial = velocity_initial
        self.position = position_initial
        self.velocity = velocity_initial
        self.acceleration = [0, 0]
        
    def euclideanDistance (self, from_position):
        r = (((self.position[0] - from_position[0])**2) + ((self.position[1] - from_position[1])**2))**(1/2)
        return r
    
    def updatePosition (self, x_new, y_new):
        self.position[0] = x_new
        self.position[1] = y_new
        
    def updateVelocity (self, vx_new, vy_new):
        self.velocity[0] = vx_new
        self.velocity[1] = vy_new
    
    def updateAcceleration (self, ax_new, ay_new):
        self.acceleration[0] = ax_new
        self.acceleration[1] = ay_new
    

class TwoBodySystem:
    def __init__ (self, reference_object2D, orbiting_object2D):
        self.reference_object2D = reference_object2D
        self.orbiting_object2D = orbiting_object2D
        self.initialize_objects()
        self.euler_trajectory = [[], []]
        self.verlet_trajectory = [[], []]
        self.ronald_ruth_3rdorder_trajectory = [[], []]
        self.ronald_ruth_4thorder_trajectory = [[], []]

    def initialize_objects(self):
        self.reference_object2D.position = self.reference_object2D.position_initial
        self.orbiting_object2D.position = self.orbiting_object2D.position_initial
        self.reference_object2D.velocity = self.reference_object2D.velocity_initial
        self.orbiting_object2D.velocity = self.orbiting_object2D.velocity_initial
        self.updateAcceleration()

    def updateAcceleration(self):
        acceleration_reference_orbiting_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (self.orbiting_object2D.position[0] - self.reference_object2D.position[0])) / ((self.orbiting_object2D.euclideanDistance(self.reference_object2D.position))**3)
        acceleration_reference_orbiting_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (self.orbiting_object2D.position[1] - self.reference_object2D.position[1])) / ((self.orbiting_object2D.euclideanDistance(self.reference_object2D.position))**3)
        self.orbiting_object2D.updateAcceleration(acceleration_reference_orbiting_x, acceleration_reference_orbiting_y)

    def euler_cromer_method(self, stepsize, num_iterations):
        self.initialize_objects()
        self.euler_trajectory = [[self.orbiting_object2D.position_initial[0]], [self.orbiting_object2D.position_initial[1]]]
        
        x_prev, y_prev = self.orbiting_object2D.position_initial
        vx_prev, vy_prev = self.orbiting_object2D.velocity_initial

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            
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
            
            self.orbiting_object2D.updatePosition(x_prev + (vx_prev * stepsize) + ((ax_prev * (stepsize**2))/2), y_prev + (vy_prev * stepsize) + ((ay_prev * (stepsize**2))/2))
            self.updateAcceleration()
            self.orbiting_object2D.updateVelocity(vx_prev + (((ax_prev + self.orbiting_object2D.acceleration[0]) * stepsize)/2), vy_prev + (((ay_prev + self.orbiting_object2D.acceleration[1]) * stepsize)/2))

            self.verlet_trajectory[0].append(self.orbiting_object2D.position[0])
            self.verlet_trajectory[1].append(self.orbiting_object2D.position[1])

            x_prev, y_prev = self.orbiting_object2D.position
            vx_prev, vy_prev = self.orbiting_object2D.velocity
            ax_prev, ay_prev = self.orbiting_object2D.acceleration

    def ronald_ruth_3rdorder_method (self, stepsize, num_iterations):
        self.initialize_objects()
        self.ronald_ruth_3rdorder_trajectory = [[self.orbiting_object2D.position_initial[0]], [self.orbiting_object2D.position_initial[1]]]

        coefficients = [[1, -2/3, 2/3], [-1/24, 3/4, 7/24]]

        x_prev, y_prev = self.orbiting_object2D.position_initial
        vx_prev, vy_prev = self.orbiting_object2D.velocity_initial
        ax_prev, ay_prev = self.orbiting_object2D.acceleration

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            
            x_array = [x_prev for i in range(len(coefficients[0]) + 1)]
            y_array = [y_prev for i in range(len(coefficients[0]) + 1)]
            vx_array = [vx_prev for i in range(len(coefficients[0]) + 1)]
            vy_array = [vy_prev for i in range(len(coefficients[0]) + 1)]
            ax_array = [ax_prev for i in range(len(coefficients[0]) + 1)]
            ay_array = [ay_prev for i in range(len(coefficients[0]) + 1)]

            for i in range(len(coefficients[0])):
                vx_array[i+1] = vx_array[i] + (coefficients[1][i] * ax_array[i] * stepsize)
                vy_array[i+1] = vy_array[i] + (coefficients[1][i] * ay_array[i] * stepsize)
                x_array[i+1] = x_array[i] + (coefficients[0][i] * vx_array[i+1] * stepsize)
                y_array[i+1] = y_array[i] + (coefficients[0][i] * vy_array[i+1] * stepsize)
                ax_array[i+1] = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (x_array[i+1] - self.reference_object2D.position[0])) / ((((x_array[i+1] - self.reference_object2D.position[0])**2) + ((y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))
                ay_array[i+1] = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (y_array[i+1] - self.reference_object2D.position[1])) / ((((x_array[i+1] - self.reference_object2D.position[0])**2) + ((y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))

            self.orbiting_object2D.updatePosition(x_array[len(coefficients[0])], y_array[len(coefficients[0])])
            self.orbiting_object2D.updateVelocity(vx_array[len(coefficients[0])], vy_array[len(coefficients[0])])
            self.orbiting_object2D.updateAcceleration(ax_array[len(coefficients[0])], ay_array[len(coefficients[0])])

            self.ronald_ruth_3rdorder_trajectory[0].append(self.orbiting_object2D.position[0])
            self.ronald_ruth_3rdorder_trajectory[1].append(self.orbiting_object2D.position[1])

            x_prev, y_prev = self.orbiting_object2D.position
            vx_prev, vy_prev = self.orbiting_object2D.velocity
            ax_prev, ay_prev = self.orbiting_object2D.acceleration

    def ronald_ruth_4thorder_method (self, stepsize, num_iterations):
        self.initialize_objects()
        self.ronald_ruth_4thorder_trajectory = [[self.orbiting_object2D.position_initial[0]], [self.orbiting_object2D.position_initial[1]]]

        coefficients = [[1/(2*(2-(2**(1/3)))), (1-(2**(1/3)))/(2*(2-(2**(1/3)))), (1-(2**(1/3)))/(2*(2-(2**(1/3)))), 1/(2*(2-(2**(1/3))))], [1/(2*(2-(2**(1/3)))), -(2**(1/3))/(2*(2-(2**(1/3)))), 1/(2*(2-(2**(1/3)))), 0]]

        x_prev, y_prev = self.orbiting_object2D.position_initial
        vx_prev, vy_prev = self.orbiting_object2D.velocity_initial
        ax_prev, ay_prev = self.orbiting_object2D.acceleration

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            
            x_array = [x_prev for i in range(len(coefficients[0]) + 1)]
            y_array = [y_prev for i in range(len(coefficients[0]) + 1)]
            vx_array = [vx_prev for i in range(len(coefficients[0]) + 1)]
            vy_array = [vy_prev for i in range(len(coefficients[0]) + 1)]
            ax_array = [ax_prev for i in range(len(coefficients[0]) + 1)]
            ay_array = [ay_prev for i in range(len(coefficients[0]) + 1)]

            for i in range(len(coefficients[0])):
                vx_array[i+1] = vx_array[i] + (coefficients[1][i] * ax_array[i] * stepsize)
                vy_array[i+1] = vy_array[i] + (coefficients[1][i] * ay_array[i] * stepsize)
                x_array[i+1] = x_array[i] + (coefficients[0][i] * vx_array[i+1] * stepsize)
                y_array[i+1] = y_array[i] + (coefficients[0][i] * vy_array[i+1] * stepsize)
                ax_array[i+1] = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (x_array[i+1] - self.reference_object2D.position[0])) / ((((x_array[i+1] - self.reference_object2D.position[0])**2) + ((y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))
                ay_array[i+1] = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (y_array[i+1] - self.reference_object2D.position[1])) / ((((x_array[i+1] - self.reference_object2D.position[0])**2) + ((y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))

            self.orbiting_object2D.updatePosition(x_array[len(coefficients[0])], y_array[len(coefficients[0])])
            self.orbiting_object2D.updateVelocity(vx_array[len(coefficients[0])], vy_array[len(coefficients[0])])
            self.orbiting_object2D.updateAcceleration(ax_array[len(coefficients[0])], ay_array[len(coefficients[0])])

            self.ronald_ruth_4thorder_trajectory[0].append(self.orbiting_object2D.position[0])
            self.ronald_ruth_4thorder_trajectory[1].append(self.orbiting_object2D.position[1])

            x_prev, y_prev = self.orbiting_object2D.position
            vx_prev, vy_prev = self.orbiting_object2D.velocity
            ax_prev, ay_prev = self.orbiting_object2D.acceleration

    def plot_euler_trajectory(self):
        plt.plot(self.euler_trajectory[0], self.euler_trajectory[1], label = self.orbiting_object2D.name)
        plt.legend(loc = 'best')
        plt.show()

    def plot_verlet_trajectory(self):
        plt.plot(self.verlet_trajectory[0], self.verlet_trajectory[1], label = self.orbiting_object2D.name)
        plt.legend(loc = 'best')
        plt.show()

    def plot_ronald_ruth_3rdorder_trajectory(self):
        plt.plot(self.ronald_ruth_3rdorder_trajectory[0], self.ronald_ruth_3rdorder_trajectory[1], label = self.orbiting_object2D.name)
        plt.legend(loc = 'best')
        plt.show()

    def plot_ronald_ruth_3rdorder_trajectory(self):
        plt.plot(self.ronald_ruth_3rdorder_trajectory[0], self.ronald_ruth_3rdorder_trajectory[1], label = self.orbiting_object2D.name)
        plt.legend(loc = 'best')
        plt.show()

    def visualize_trajectory(self):
        pygame.init()
        win = pygame.display.set_mode((720,720))
        radius=7.5
        time_stamp=0
        length=len(self.euler_trajectory[0])
        run=True
        pygame.draw.circle(win,(255,255,0),(360,360),radius)
        while run:
            pygame.time.delay(50)

            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    run = False
            pygame.draw.circle(win, (0, 200, 255), (self.euler_trajectory[0][time_stamp]*300/DISTANCE["SUN_"+self.orbiting_object2D.name] + 360 , self.euler_trajectory[1][time_stamp]*300/DISTANCE["SUN_"+self.orbiting_object2D.name] + 360), radius/3)

            pygame.display.update()
            if time_stamp<length-2:
                time_stamp+=1

        pygame.quit()


class ThreeBodySystem:
    def __init__ (self, reference_object2D, orbiting_object2D_1, orbiting_object2D_2):
        self.reference_object2D = reference_object2D
        self.orbiting_object2D_1 = orbiting_object2D_1
        self.orbiting_object2D_2 = orbiting_object2D_2
        self.initialize_objects()
        self.euler_trajectory = [[[], []], [[], []]]
        self.verlet_trajectory = [[[], []], [[], []]]
        self.ronald_ruth_3rdorder_trajectory = [[[], []], [[], []]]
        self.ronald_ruth_4thorder_trajectory = [[[], []], [[], []]]

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
        
        self.orbiting_object2D_1.updateAcceleration(acceleration_reference_orbiting1_x + acceleration_orbiting2_orbiting1_x, acceleration_reference_orbiting1_y + acceleration_orbiting2_orbiting1_y)
        self.orbiting_object2D_2.updateAcceleration(acceleration_reference_orbiting2_x + acceleration_orbiting1_orbiting2_x, acceleration_reference_orbiting2_y + acceleration_orbiting1_orbiting2_y)


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

    def ronald_ruth_3rdorder_method (self, stepsize, num_iterations):
        self.initialize_objects()
        self.ronald_ruth_3rdorder_trajectory = [[[self.orbiting_object2D_1.position_initial[0]], [self.orbiting_object2D_1.position_initial[1]]], [[self.orbiting_object2D_2.position_initial[0]], [self.orbiting_object2D_2.position_initial[1]]]]

        coefficients = [[1, -2/3, 2/3], [-1/24, 3/4, 7/24]]

        orbiting1_x_prev, orbiting1_y_prev = self.orbiting_object2D_1.position_initial
        orbiting2_x_prev, orbiting2_y_prev = self.orbiting_object2D_2.position_initial
        orbiting1_vx_prev, orbiting1_vy_prev = self.orbiting_object2D_1.velocity_initial
        orbiting2_vx_prev, orbiting2_vy_prev = self.orbiting_object2D_2.velocity_initial
        orbiting1_ax_prev, orbiting1_ay_prev = self.orbiting_object2D_1.acceleration
        orbiting2_ax_prev, orbiting2_ay_prev = self.orbiting_object2D_2.acceleration

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            
            orbiting1_x_array = [orbiting1_x_prev for i in range(len(coefficients[0]) + 1)]
            orbiting1_y_array = [orbiting1_y_prev for i in range(len(coefficients[0]) + 1)]
            orbiting1_vx_array = [orbiting1_vx_prev for i in range(len(coefficients[0]) + 1)]
            orbiting1_vy_array = [orbiting1_vy_prev for i in range(len(coefficients[0]) + 1)]
            orbiting1_ax_array = [orbiting1_ax_prev for i in range(len(coefficients[0]) + 1)]
            orbiting1_ay_array = [orbiting1_ay_prev for i in range(len(coefficients[0]) + 1)]

            orbiting2_x_array = [orbiting2_x_prev for i in range(len(coefficients[0]) + 1)]
            orbiting2_y_array = [orbiting2_y_prev for i in range(len(coefficients[0]) + 1)]
            orbiting2_vx_array = [orbiting2_vx_prev for i in range(len(coefficients[0]) + 1)]
            orbiting2_vy_array = [orbiting2_vy_prev for i in range(len(coefficients[0]) + 1)]
            orbiting2_ax_array = [orbiting2_ax_prev for i in range(len(coefficients[0]) + 1)]
            orbiting2_ay_array = [orbiting2_ay_prev for i in range(len(coefficients[0]) + 1)]

            for i in range(len(coefficients[0])):
                orbiting1_vx_array[i+1] = orbiting1_vx_array[i] + (coefficients[1][i] * orbiting1_ax_array[i] * stepsize)
                orbiting1_vy_array[i+1] = orbiting1_vy_array[i] + (coefficients[1][i] * orbiting1_ay_array[i] * stepsize)
                orbiting2_vx_array[i+1] = orbiting2_vx_array[i] + (coefficients[1][i] * orbiting2_ax_array[i] * stepsize)
                orbiting2_vy_array[i+1] = orbiting2_vy_array[i] + (coefficients[1][i] * orbiting2_ay_array[i] * stepsize)

                orbiting1_x_array[i+1] = orbiting1_x_array[i] + (coefficients[0][i] * orbiting1_vx_array[i+1] * stepsize)
                orbiting1_y_array[i+1] = orbiting1_y_array[i] + (coefficients[0][i] * orbiting1_vy_array[i+1] * stepsize)
                orbiting2_x_array[i+1] = orbiting2_x_array[i] + (coefficients[0][i] * orbiting2_vx_array[i+1] * stepsize)
                orbiting2_y_array[i+1] = orbiting2_y_array[i] + (coefficients[0][i] * orbiting2_vy_array[i+1] * stepsize)

                acc_reference_orbiting1_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting1_x_array[i+1] - self.reference_object2D.position[0])) / ((((orbiting1_x_array[i+1] - self.reference_object2D.position[0])**2) + ((orbiting1_y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))
                acc_orbiting2_orbiting1_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_2.name] * (orbiting1_x_array[i+1] - orbiting2_x_array[i+1])) / ((((orbiting1_x_array[i+1] - orbiting2_x_array[i+1])**2) + ((orbiting1_y_array[i+1] - orbiting2_y_array[i+1])**2))**(3/2))
                acc_reference_orbiting1_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting1_y_array[i+1] - self.reference_object2D.position[1])) / ((((orbiting1_x_array[i+1] - self.reference_object2D.position[0])**2) + ((orbiting1_y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))
                acc_orbiting2_orbiting1_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_2.name] * (orbiting1_y_array[i+1] - orbiting2_2_array[i+1])) / ((((orbiting1_x_array[i+1] - orbiting2_x_array[i+1])**2) + ((orbiting1_y_array[i+1] - orbiting2_y_array[i+1])**2))**(3/2))

                acc_reference_orbiting2_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting2_x_array[i+1] - self.reference_object2D.position[0])) / ((((orbiting2_x_array[i+1] - self.reference_object2D.position[0])**2) + ((orbiting2_y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))
                acc_orbiting1_orbiting2_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_1.name] * (orbiting2_x_array[i+1] - orbiting1_x_array[i+1])) / ((((orbiting1_x_array[i+1] - orbiting2_x_array[i+1])**2) + ((orbiting1_y_array[i+1] - orbiting2_y_array[i+1])**2))**(3/2))
                acc_reference_orbiting2_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting2_y_array[i+1] - self.reference_object2D.position[1])) / ((((orbiting2_x_array[i+1] - self.reference_object2D.position[0])**2) + ((orbiting2_y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))
                acc_orbiting1_orbiting2_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_1.name] * (orbiting2_y_array[i+1] - orbiting1_y_array[i+1])) / ((((orbiting1_x_array[i+1] - orbiting2_x_array[i+1])**2) + ((orbiting1_y_array[i+1] - orbiting2_y_array[i+1])**2))**(3/2))

                orbiting1_ax_array[i+1] = acc_reference_orbiting1_x + acc_orbiting2_orbiting1_x 
                orbiting1_ay_array[i+1] = acc_reference_orbiting1_y + acc_orbiting2_orbiting1_y
                orbiting2_ax_array[i+1] = acc_reference_orbiting2_x + acc_orbiting1_orbiting2_x
                orbiting2_ay_array[i+1] = acc_reference_orbiting1_y + acc_orbiting2_orbiting1_y

            self.orbiting1_object2D.updatePosition(orbiting1_x_array[len(coefficients[0])], orbiting1_y_array[len(coefficients[0])])
            self.orbiting1_object2D.updateVelocity(orbiting1_vx_array[len(coefficients[0])], orbiting1_vy_array[len(coefficients[0])])
            self.orbiting1_object2D.updateAcceleration(orbiting1_ax_array[len(coefficients[0])], orbiting1_ay_array[len(coefficients[0])])

            self.orbiting2_object2D.updatePosition(orbiting2_x_array[len(coefficients[0])], orbiting2_y_array[len(coefficients[0])])
            self.orbiting2_object2D.updateVelocity(orbiting2_vx_array[len(coefficients[0])], orbiting2_vy_array[len(coefficients[0])])
            self.orbiting2_object2D.updateAcceleration(orbiting2_ax_array[len(coefficients[0])], orbiting2_ay_array[len(coefficients[0])])

            self.ronald_ruth_3rdorder_trajectory[0][0].append(self.orbiting_object2D_1.position[0])
            self.ronald_ruth_3rdorder_trajectory[0][1].append(self.orbiting_object2D_1.position[1])
            self.ronald_ruth_3rdorder_trajectory[1][0].append(self.orbiting_object2D_2.position[0])
            self.ronald_ruth_3rdorder_trajectory[1][1].append(self.orbiting_object2D_2.position[1])

            orbiting1_x_prev, orbiting1_y_prev = self.orbiting_object2D_1.position
            orbiting2_x_prev, orbiting2_y_prev = self.orbiting_object2D_2.position
            orbiting1_vx_prev, orbiting1_vy_prev = self.orbiting_object2D_1.velocity
            orbiting2_vx_prev, orbiting2_vy_prev = self.orbiting_object2D_2.velocity
            orbiting1_ax_prev, orbiting1_ay_prev = self.orbiting_object2D_1.acceleration
            orbiting2_ax_prev, orbiting2_ay_prev = self.orbiting_object2D_2.acceleration
    

    def ronald_ruth_4thorder_method (self, stepsize, num_iterations):
        self.initialize_objects()
        self.ronald_ruth_4thorder_trajectory = [[[self.orbiting_object2D_1.position_initial[0]], [self.orbiting_object2D_1.position_initial[1]]], [[self.orbiting_object2D_2.position_initial[0]], [self.orbiting_object2D_2.position_initial[1]]]]

        coefficients = [[1/(2*(2-(2**(1/3)))), (1-(2**(1/3)))/(2*(2-(2**(1/3)))), (1-(2**(1/3)))/(2*(2-(2**(1/3)))), 1/(2*(2-(2**(1/3))))], [1/(2*(2-(2**(1/3)))), -(2**(1/3))/(2*(2-(2**(1/3)))), 1/(2*(2-(2**(1/3)))), 0]]

        orbiting1_x_prev, orbiting1_y_prev = self.orbiting_object2D_1.position_initial
        orbiting2_x_prev, orbiting2_y_prev = self.orbiting_object2D_2.position_initial
        orbiting1_vx_prev, orbiting1_vy_prev = self.orbiting_object2D_1.velocity_initial
        orbiting2_vx_prev, orbiting2_vy_prev = self.orbiting_object2D_2.velocity_initial
        orbiting1_ax_prev, orbiting1_ay_prev = self.orbiting_object2D_1.acceleration
        orbiting2_ax_prev, orbiting2_ay_prev = self.orbiting_object2D_2.acceleration

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            
            orbiting1_x_array = [orbiting1_x_prev for i in range(len(coefficients[0]) + 1)]
            orbiting1_y_array = [orbiting1_y_prev for i in range(len(coefficients[0]) + 1)]
            orbiting1_vx_array = [orbiting1_vx_prev for i in range(len(coefficients[0]) + 1)]
            orbiting1_vy_array = [orbiting1_vy_prev for i in range(len(coefficients[0]) + 1)]
            orbiting1_ax_array = [orbiting1_ax_prev for i in range(len(coefficients[0]) + 1)]
            orbiting1_ay_array = [orbiting1_ay_prev for i in range(len(coefficients[0]) + 1)]

            orbiting2_x_array = [orbiting2_x_prev for i in range(len(coefficients[0]) + 1)]
            orbiting2_y_array = [orbiting2_y_prev for i in range(len(coefficients[0]) + 1)]
            orbiting2_vx_array = [orbiting2_vx_prev for i in range(len(coefficients[0]) + 1)]
            orbiting2_vy_array = [orbiting2_vy_prev for i in range(len(coefficients[0]) + 1)]
            orbiting2_ax_array = [orbiting2_ax_prev for i in range(len(coefficients[0]) + 1)]
            orbiting2_ay_array = [orbiting2_ay_prev for i in range(len(coefficients[0]) + 1)]

            for i in range(len(coefficients[0])):
                orbiting1_vx_array[i+1] = orbiting1_vx_array[i] + (coefficients[1][i] * orbiting1_ax_array[i] * stepsize)
                orbiting1_vy_array[i+1] = orbiting1_vy_array[i] + (coefficients[1][i] * orbiting1_ay_array[i] * stepsize)
                orbiting2_vx_array[i+1] = orbiting2_vx_array[i] + (coefficients[1][i] * orbiting2_ax_array[i] * stepsize)
                orbiting2_vy_array[i+1] = orbiting2_vy_array[i] + (coefficients[1][i] * orbiting2_ay_array[i] * stepsize)

                orbiting1_x_array[i+1] = orbiting1_x_array[i] + (coefficients[0][i] * orbiting1_vx_array[i+1] * stepsize)
                orbiting1_y_array[i+1] = orbiting1_y_array[i] + (coefficients[0][i] * orbiting1_vy_array[i+1] * stepsize)
                orbiting2_x_array[i+1] = orbiting2_x_array[i] + (coefficients[0][i] * orbiting2_vx_array[i+1] * stepsize)
                orbiting2_y_array[i+1] = orbiting2_y_array[i] + (coefficients[0][i] * orbiting2_vy_array[i+1] * stepsize)

                acc_reference_orbiting1_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting1_x_array[i+1] - self.reference_object2D.position[0])) / ((((orbiting1_x_array[i+1] - self.reference_object2D.position[0])**2) + ((orbiting1_y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))
                acc_orbiting2_orbiting1_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_2.name] * (orbiting1_x_array[i+1] - orbiting2_x_array[i+1])) / ((((orbiting1_x_array[i+1] - orbiting2_x_array[i+1])**2) + ((orbiting1_y_array[i+1] - orbiting2_y_array[i+1])**2))**(3/2))
                acc_reference_orbiting1_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting1_y_array[i+1] - self.reference_object2D.position[1])) / ((((orbiting1_x_array[i+1] - self.reference_object2D.position[0])**2) + ((orbiting1_y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))
                acc_orbiting2_orbiting1_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_2.name] * (orbiting1_y_array[i+1] - orbiting2_2_array[i+1])) / ((((orbiting1_x_array[i+1] - orbiting2_x_array[i+1])**2) + ((orbiting1_y_array[i+1] - orbiting2_y_array[i+1])**2))**(3/2))

                acc_reference_orbiting2_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting2_x_array[i+1] - self.reference_object2D.position[0])) / ((((orbiting2_x_array[i+1] - self.reference_object2D.position[0])**2) + ((orbiting2_y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))
                acc_orbiting1_orbiting2_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_1.name] * (orbiting2_x_array[i+1] - orbiting1_x_array[i+1])) / ((((orbiting1_x_array[i+1] - orbiting2_x_array[i+1])**2) + ((orbiting1_y_array[i+1] - orbiting2_y_array[i+1])**2))**(3/2))
                acc_reference_orbiting2_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting2_y_array[i+1] - self.reference_object2D.position[1])) / ((((orbiting2_x_array[i+1] - self.reference_object2D.position[0])**2) + ((orbiting2_y_array[i+1] - self.reference_object2D.position[1])**2))**(3/2))
                acc_orbiting1_orbiting2_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_object2D_1.name] * (orbiting2_y_array[i+1] - orbiting1_y_array[i+1])) / ((((orbiting1_x_array[i+1] - orbiting2_x_array[i+1])**2) + ((orbiting1_y_array[i+1] - orbiting2_y_array[i+1])**2))**(3/2))

                orbiting1_ax_array[i+1] = acc_reference_orbiting1_x + acc_orbiting2_orbiting1_x 
                orbiting1_ay_array[i+1] = acc_reference_orbiting1_y + acc_orbiting2_orbiting1_y
                orbiting2_ax_array[i+1] = acc_reference_orbiting2_x + acc_orbiting1_orbiting2_x
                orbiting2_ay_array[i+1] = acc_reference_orbiting1_y + acc_orbiting2_orbiting1_y

            self.orbiting1_object2D.updatePosition(orbiting1_x_array[len(coefficients[0])], orbiting1_y_array[len(coefficients[0])])
            self.orbiting1_object2D.updateVelocity(orbiting1_vx_array[len(coefficients[0])], orbiting1_vy_array[len(coefficients[0])])
            self.orbiting1_object2D.updateAcceleration(orbiting1_ax_array[len(coefficients[0])], orbiting1_ay_array[len(coefficients[0])])

            self.orbiting2_object2D.updatePosition(orbiting2_x_array[len(coefficients[0])], orbiting2_y_array[len(coefficients[0])])
            self.orbiting2_object2D.updateVelocity(orbiting2_vx_array[len(coefficients[0])], orbiting2_vy_array[len(coefficients[0])])
            self.orbiting2_object2D.updateAcceleration(orbiting2_ax_array[len(coefficients[0])], orbiting2_ay_array[len(coefficients[0])])

            self.ronald_ruth_4thorder_trajectory[0][0].append(self.orbiting_object2D_1.position[0])
            self.ronald_ruth_4thorder_trajectory[0][1].append(self.orbiting_object2D_1.position[1])
            self.ronald_ruth_4thorder_trajectory[1][0].append(self.orbiting_object2D_2.position[0])
            self.ronald_ruth_4thorder_trajectory[1][1].append(self.orbiting_object2D_2.position[1])

            orbiting1_x_prev, orbiting1_y_prev = self.orbiting_object2D_1.position
            orbiting2_x_prev, orbiting2_y_prev = self.orbiting_object2D_2.position
            orbiting1_vx_prev, orbiting1_vy_prev = self.orbiting_object2D_1.velocity
            orbiting2_vx_prev, orbiting2_vy_prev = self.orbiting_object2D_2.velocity
            orbiting1_ax_prev, orbiting1_ay_prev = self.orbiting_object2D_1.acceleration
            orbiting2_ax_prev, orbiting2_ay_prev = self.orbiting_object2D_2.acceleration

    def plot_euler_trajectory(self):
        plt.plot(self.euler_trajectory[0][0], self.euler_trajectory[0][1], color = 'r', label = self.orbiting_object2D_1.name)
        plt.plot(self.euler_trajectory[1][0], self.euler_trajectory[1][1], color = 'g', label = self.orbiting_object2D_2.name)
        plt.legend(loc = 'best')
        plt.show()

    def plot_verlet_trajectory(self):
        plt.plot(self.verlet_trajectory[0][0], self.verlet_trajectory[0][1], color = 'r', label = self.orbiting_object2D_1.name)
        plt.plot(self.verlet_trajectory[1][0], self.verlet_trajectory[1][1], color = 'g', label = self.orbiting_object2D_2.name)
        plt.legend(loc = 'best')
        plt.show()

    def plot_ronald_ruth_3rd_trajectory(self):
        plt.plot(self.ronald_ruth_3rd_trajectory[0][0], self.ronald_ruth_3rd_trajectory[0][1], color = 'r', label = self.orbiting_object2D_1.name)
        plt.plot(self.ronald_ruth_3rd_trajectory[1][0], self.ronald_ruth_3rd_trajectory[1][1], color = 'g', label = self.orbiting_object2D_2.name)
        plt.legend(loc = 'best')
        plt.show()

    def plot_ronald_ruth_4th_trajectory(self):
        plt.plot(self.ronald_ruth_4th_trajectory[0][0], self.ronald_ruth_4th_trajectory[0][1], color = 'r', label = self.orbiting_object2D_1.name)
        plt.plot(self.ronald_ruth_4th_trajectory[1][0], self.ronald_ruth_4th_trajectory[1][1], color = 'g', label = self.orbiting_object2D_2.name)
        plt.legend(loc = 'best')
        plt.show()

    def visualize_trajectory(self):
        pygame.init()
        win = pygame.display.set_mode((720,720))
        radius=7.5
        time_stamp=0
        length=len(self.euler_trajectory[0][0])
        run=True
        pygame.draw.circle(win,(255,255,0),(360,360),radius)

        scaling=max(DISTANCE["SUN_"+self.orbiting_object2D_1.name],DISTANCE["SUN_"+self.orbiting_object2D_2.name])

        while run:
            pygame.time.delay(20)

            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    run = False
            pygame.draw.circle(win, (255, 0, 255), (self.euler_trajectory[0][0][time_stamp]*300/scaling + 360 , self.euler_trajectory[0][1][time_stamp]*300/scaling + 360), radius/3)
            pygame.draw.circle(win, (0, 200, 255), (self.euler_trajectory[1][0][time_stamp]*300/scaling + 360 , self.euler_trajectory[1][1][time_stamp]*300/scaling + 360), radius/3)

            pygame.display.update()
            if time_stamp<length-2:
                time_stamp+=1

        pygame.quit()

class NBodySystem:
    def __init__ (self, reference_object2D, orbiting_objects2D):
        self.reference_object2D = reference_object2D
        self.orbiting_objects2D = orbiting_objects2D
        self.initialize_objects()
        self.euler_trajectory = [[[], []] for i in range(len(orbiting_objects2D))]
        self.verlet_trajectory = [[[], []] for i in range(len(orbiting_objects2D))]
        self.ronald_ruth_3rdorder_trajectory = [[[], []] for i in range(len(orbiting_objects2D))]
        self.ronald_ruth_4thorder_trajectory = [[[], []] for i in range(len(orbiting_objects2D))]

    def initialize_objects(self):
        self.reference_object2D.position = self.reference_object2D.position_initial
        self.reference_object2D.velocity = self.reference_object2D.velocity_initial
        for i in range(len(self.orbiting_objects2D)):
            self.orbiting_objects2D[i].position = self.orbiting_objects2D[i].position_initial
            self.orbiting_objects2D[i].velocity = self.orbiting_objects2D[i].velocity_initial
        self.updateAcceleration()

    def updateAcceleration(self):
        for i in range(len(self.orbiting_objects2D)):
            acceleration_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (self.orbiting_objects2D[i].position[0] - self.reference_object2D.position[0])) / ((self.orbiting_objects2D[i].euclideanDistance(self.reference_object2D.position))**3)
            acceleration_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (self.orbiting_objects2D[i].position[1] - self.reference_object2D.position[1])) / ((self.orbiting_objects2D[i].euclideanDistance(self.reference_object2D.position))**3)

            for j in range(len(self.orbiting_objects2D)):
                if (j != i):
                    acceleration_x += -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_objects2D[j].name] * (self.orbiting_objects2D[i].position[0] - self.orbiting_objects2D[j].position[0])) / ((self.orbiting_objects2D[i].euclideanDistance(self.orbiting_objects2D[j].position))**3)
                    acceleration_y += -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_objects2D[j].name] * (self.orbiting_objects2D[i].position[1] - self.orbiting_objects2D[j].position[1])) / ((self.orbiting_objects2D[i].euclideanDistance(self.orbiting_objects2D[j].position))**3)

            self.orbiting_objects2D[i].updateAcceleration(acceleration_x, acceleration_y)

    def euler_cromer_method(self, stepsize, num_iterations):
        self.initialize_objects()
        self.euler_trajectory = [[[self.orbiting_objects2D[i].position_initial[0]], [self.orbiting_objects2D[i].position_initial[1]]] for i in range(len(self.orbiting_objects2D))]
        
        orbiting_x_prev = [self.orbiting_objects2D[i].position_initial[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_y_prev = [self.orbiting_objects2D[i].position_initial[1] for i in range(len(self.orbiting_objects2D))]
        orbiting_vx_prev = [self.orbiting_objects2D[i].velocity_initial[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_vy_prev = [self.orbiting_objects2D[i].velocity_initial[1] for i in range(len(self.orbiting_objects2D))]

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1

            for i in range(len(self.orbiting_objects2D)):
                self.orbiting_objects2D[i].updateVelocity(orbiting_vx_prev[i] + (self.orbiting_objects2D[i].acceleration[0] * stepsize), orbiting_vy_prev[i] + (self.orbiting_objects2D[i].acceleration[1] * stepsize))

            for i in range(len(self.orbiting_objects2D)):
                self.orbiting_objects2D[i].updatePosition(orbiting_x_prev[i] + (self.orbiting_objects2D[i].velocity[0] * stepsize), orbiting_y_prev[i] + (self.orbiting_objects2D[i].velocity[1] * stepsize))

            # It is working here as non-energy preserving euler method
            # for i in range(len(self.orbiting_objects2D)):
            #     self.orbiting_objects2D[i].updatePosition(orbiting_x_prev[i] + (orbiting_vx_prev[i] * stepsize), orbiting_y_prev[i] + (orbiting_vy_prev[i] * stepsize))
            
            self.updateAcceleration()

            for i in range(len(self.orbiting_objects2D)):
                self.euler_trajectory[i][0].append(self.orbiting_objects2D[i].position[0])
                self.euler_trajectory[i][1].append(self.orbiting_objects2D[i].position[1])

            orbiting_x_prev = [self.orbiting_objects2D[i].position[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_y_prev = [self.orbiting_objects2D[i].position[1] for i in range(len(self.orbiting_objects2D))]
            orbiting_vx_prev = [self.orbiting_objects2D[i].velocity[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_vy_prev = [self.orbiting_objects2D[i].velocity[1] for i in range(len(self.orbiting_objects2D))]

    def verlet_method (self, stepsize, num_iterations):
        self.initialize_objects()
        self.verlet_trajectory = [[[self.orbiting_objects2D[i].position_initial[0]], [self.orbiting_objects2D[i].position_initial[1]]] for i in range(len(self.orbiting_objects2D))]
        
        orbiting_x_prev = [self.orbiting_objects2D[i].position_initial[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_y_prev = [self.orbiting_objects2D[i].position_initial[1] for i in range(len(self.orbiting_objects2D))]
        orbiting_vx_prev = [self.orbiting_objects2D[i].velocity_initial[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_vy_prev = [self.orbiting_objects2D[i].velocity_initial[1] for i in range(len(self.orbiting_objects2D))]
        orbiting_ax_prev = [self.orbiting_objects2D[i].acceleration[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_ay_prev = [self.orbiting_objects2D[i].acceleration[1] for i in range(len(self.orbiting_objects2D))]
        
        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            
            for i in range(len(self.orbiting_objects2D)):
                self.orbiting_objects2D[i].updatePosition(orbiting_x_prev[i] + (orbiting_vx_prev[i] * stepsize) + ((orbiting_ax_prev[i] * (stepsize**2))/2), orbiting_y_prev[i] + (orbiting_vy_prev[i] * stepsize) + ((orbiting_ay_prev[i] * (stepsize**2))/2))    

            self.updateAcceleration()

            for i in range(len(self.orbiting_objects2D)):
                self.orbiting_objects2D[i].updateVelocity(orbiting_vx_prev[i] + (((orbiting_ax_prev[i] + self.orbiting_objects2D[i].acceleration[0]) * stepsize)/2), orbiting_vy_prev[i] + (((orbiting_ay_prev[i] + self.orbiting_objects2D[i].acceleration[1]) * stepsize)/2))    
            
            for i in range(len(self.orbiting_objects2D)):
                self.verlet_trajectory[i][0].append(self.orbiting_objects2D[i].position[0])
                self.verlet_trajectory[i][1].append(self.orbiting_objects2D[i].position[1])

            orbiting_x_prev = [self.orbiting_objects2D[i].position[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_y_prev = [self.orbiting_objects2D[i].position[1] for i in range(len(self.orbiting_objects2D))]
            orbiting_vx_prev = [self.orbiting_objects2D[i].velocity[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_vy_prev = [self.orbiting_objects2D[i].velocity[1] for i in range(len(self.orbiting_objects2D))]
            orbiting_ax_prev = [self.orbiting_objects2D[i].acceleration[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_ay_prev = [self.orbiting_objects2D[i].acceleration[1] for i in range(len(self.orbiting_objects2D))]

    def ronald_ruth_3rdorder_method (self, stepsize, num_iterations):
        self.initialize_objects()
        self.ronald_ruth_3rdorder_trajectory = [[[self.orbiting_objects2D[i].position_initial[0]], [self.orbiting_objects2D[i].position_initial[1]]] for i in range(len(self.orbiting_objects2D))]

        coefficients = [[1, -2/3, 2/3], [-1/24, 3/4, 7/24]]

        orbiting_x_prev = [self.orbiting_objects2D[i].position_initial[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_y_prev = [self.orbiting_objects2D[i].position_initial[1] for i in range(len(self.orbiting_objects2D))]
        orbiting_vx_prev = [self.orbiting_objects2D[i].velocity_initial[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_vy_prev = [self.orbiting_objects2D[i].velocity_initial[1] for i in range(len(self.orbiting_objects2D))]
        orbiting_ax_prev = [self.orbiting_objects2D[i].acceleration[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_ay_prev = [self.orbiting_objects2D[i].acceleration[1] for i in range(len(self.orbiting_objects2D))]

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            
            orbiting_x_array = [[orbiting_x_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]
            orbiting_y_array = [[orbiting_y_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]
            orbiting_vx_array = [[orbiting_vx_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]
            orbiting_vy_array = [[orbiting_vy_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]
            orbiting_ax_array = [[orbiting_ax_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]
            orbiting_ay_array = [[orbiting_ay_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]

            for i in range(len(coefficients[0])):

                for j in range(len(self.orbiting_objects2D)):
                    orbiting_vx_array[j][i+1] = orbiting_vx_array[j][i] + (coefficients[1][i] * orbiting_ax_array[j][i] * stepsize)
                    orbiting_vy_array[j][i+1] = orbiting_vy_array[j][i] + (coefficients[1][i] * orbiting_ay_array[j][i] * stepsize)    
                
                for j in range(len(self.orbiting_objects2D)):
                    orbiting_x_array[j][i+1] = orbiting_x_array[j][i] + (coefficients[0][i] * orbiting_vx_array[j][i+1] * stepsize)
                    orbiting_y_array[j][i+1] = orbiting_y_array[j][i] + (coefficients[0][i] * orbiting_vy_array[j][i+1] * stepsize) 

                for j in range(len(self.orbiting_objects2D)):
                    acceleration_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting_x_array[j][i+1] - self.reference_object2D.position[0])) / ((((orbiting_x_array[j][i+1] - self.reference_object2D.position[0])**2) + ((orbiting_y_array[j][i+1] - self.reference_object2D.position[1])**2))**(3/2))
                    acceleration_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting_y_array[j][i+1] - self.reference_object2D.position[1])) / ((((orbiting_x_array[j][i+1] - self.reference_object2D.position[0])**2) + ((orbiting_y_array[j][i+1] - self.reference_object2D.position[1])**2))**(3/2))

                    for k in range(len(self.orbiting_objects2D)):
                        if (k != j):
                            acceleration_x += -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_objects2D[j].name] * (orbiting_x_array[j][i+1] - orbiting_x_array[k][i+1])) / ((((orbiting_x_array[j][i+1] - orbiting_x_array[k][i+1])**2) + ((orbiting_y_array[j][i+1] - orbiting_y_array[k][i+1])**2))**(3/2))
                            acceleration_y += -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_objects2D[j].name] * (orbiting_y_array[j][i+1] - orbiting_y_array[k][i+1])) / ((((orbiting_x_array[j][i+1] - orbiting_x_array[k][i+1])**2) + ((orbiting_y_array[j][i+1] - orbiting_y_array[k][i+1])**2))**(3/2))

                    orbiting_ax_array[j][i+1] = acceleration_x
                    orbiting_ay_array[j][i+1] = acceleration_y

            for i in range(len(self.orbiting_objects2D)):
                self.orbiting_objects2D[i].updatePosition(orbiting_x_array[i][len(coefficients[0])], orbiting_y_array[i][len(coefficients[0])])
                self.orbiting_objects2D[i].updateVelocity(orbiting_vx_array[i][len(coefficients[0])], orbiting_vy_array[i][len(coefficients[0])])
                self.orbiting_objects2D[i].updateAcceleration(orbiting_ax_array[i][len(coefficients[0])], orbiting_ay_array[i][len(coefficients[0])])    

            for i in range(len(self.orbiting_objects2D)):
                self.ronald_ruth_3rdorder_trajectory[i][0].append(self.orbiting_objects2D[i].position[0])
                self.ronald_ruth_3rdorder_trajectory[i][1].append(self.orbiting_objects2D[i].position[1])

            orbiting_x_prev = [self.orbiting_objects2D[i].position[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_y_prev = [self.orbiting_objects2D[i].position[1] for i in range(len(self.orbiting_objects2D))]
            orbiting_vx_prev = [self.orbiting_objects2D[i].velocity[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_vy_prev = [self.orbiting_objects2D[i].velocity[1] for i in range(len(self.orbiting_objects2D))]
            orbiting_ax_prev = [self.orbiting_objects2D[i].acceleration[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_ay_prev = [self.orbiting_objects2D[i].acceleration[1] for i in range(len(self.orbiting_objects2D))]

    def ronald_ruth_4thorder_method (self, stepsize, num_iterations):
        self.initialize_objects()
        self.ronald_ruth_4thorder_trajectory = [[[self.orbiting_objects2D[i].position_initial[0]], [self.orbiting_objects2D[i].position_initial[1]]] for i in range(len(self.orbiting_objects2D))]

        coefficients = [[1/(2*(2-(2**(1/3)))), (1-(2**(1/3)))/(2*(2-(2**(1/3)))), (1-(2**(1/3)))/(2*(2-(2**(1/3)))), 1/(2*(2-(2**(1/3))))], [1/(2*(2-(2**(1/3)))), -(2**(1/3))/(2*(2-(2**(1/3)))), 1/(2*(2-(2**(1/3)))), 0]]

        orbiting_x_prev = [self.orbiting_objects2D[i].position_initial[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_y_prev = [self.orbiting_objects2D[i].position_initial[1] for i in range(len(self.orbiting_objects2D))]
        orbiting_vx_prev = [self.orbiting_objects2D[i].velocity_initial[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_vy_prev = [self.orbiting_objects2D[i].velocity_initial[1] for i in range(len(self.orbiting_objects2D))]
        orbiting_ax_prev = [self.orbiting_objects2D[i].acceleration[0] for i in range(len(self.orbiting_objects2D))]
        orbiting_ay_prev = [self.orbiting_objects2D[i].acceleration[1] for i in range(len(self.orbiting_objects2D))]

        iterations = 0
        while (iterations < num_iterations):
            iterations = iterations + 1
            
            orbiting_x_array = [[orbiting_x_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]
            orbiting_y_array = [[orbiting_y_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]
            orbiting_vx_array = [[orbiting_vx_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]
            orbiting_vy_array = [[orbiting_vy_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]
            orbiting_ax_array = [[orbiting_ax_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]
            orbiting_ay_array = [[orbiting_ay_prev[i] for j in range(len(coefficients[0]) + 1)] for i in range(len(self.orbiting_objects2D))]

            for i in range(len(coefficients[0])):

                for j in range(len(self.orbiting_objects2D)):
                    orbiting_vx_array[j][i+1] = orbiting_vx_array[j][i] + (coefficients[1][i] * orbiting_ax_array[j][i] * stepsize)
                    orbiting_vy_array[j][i+1] = orbiting_vy_array[j][i] + (coefficients[1][i] * orbiting_ay_array[j][i] * stepsize)    
                
                for j in range(len(self.orbiting_objects2D)):
                    orbiting_x_array[j][i+1] = orbiting_x_array[j][i] + (coefficients[0][i] * orbiting_vx_array[j][i+1] * stepsize)
                    orbiting_y_array[j][i+1] = orbiting_y_array[j][i] + (coefficients[0][i] * orbiting_vy_array[j][i+1] * stepsize) 

                for j in range(len(self.orbiting_objects2D)):
                    acceleration_x = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting_x_array[j][i+1] - self.reference_object2D.position[0])) / ((((orbiting_x_array[j][i+1] - self.reference_object2D.position[0])**2) + ((orbiting_y_array[j][i+1] - self.reference_object2D.position[1])**2))**(3/2))
                    acceleration_y = -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.reference_object2D.name] * (orbiting_y_array[j][i+1] - self.reference_object2D.position[1])) / ((((orbiting_x_array[j][i+1] - self.reference_object2D.position[0])**2) + ((orbiting_y_array[j][i+1] - self.reference_object2D.position[1])**2))**(3/2))

                    for k in range(len(self.orbiting_objects2D)):
                        if (k != j):
                            acceleration_x += -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_objects2D[j].name] * (orbiting_x_array[j][i+1] - orbiting_x_array[k][i+1])) / ((((orbiting_x_array[j][i+1] - orbiting_x_array[k][i+1])**2) + ((orbiting_y_array[j][i+1] - orbiting_y_array[k][i+1])**2))**(3/2))
                            acceleration_y += -1 * (CONVERSION_CONSTANT * GRAVITATIONAL_CONSTANT * MASS[self.orbiting_objects2D[j].name] * (orbiting_y_array[j][i+1] - orbiting_y_array[k][i+1])) / ((((orbiting_x_array[j][i+1] - orbiting_x_array[k][i+1])**2) + ((orbiting_y_array[j][i+1] - orbiting_y_array[k][i+1])**2))**(3/2))

                    orbiting_ax_array[j][i+1] = acceleration_x
                    orbiting_ay_array[j][i+1] = acceleration_y

            for i in range(len(self.orbiting_objects2D)):
                self.orbiting_objects2D[i].updatePosition(orbiting_x_array[i][len(coefficients[0])], orbiting_y_array[i][len(coefficients[0])])
                self.orbiting_objects2D[i].updateVelocity(orbiting_vx_array[i][len(coefficients[0])], orbiting_vy_array[i][len(coefficients[0])])
                self.orbiting_objects2D[i].updateAcceleration(orbiting_ax_array[i][len(coefficients[0])], orbiting_ay_array[i][len(coefficients[0])])    

            for i in range(len(self.orbiting_objects2D)):
                self.ronald_ruth_4thorder_trajectory[i][0].append(self.orbiting_objects2D[i].position[0])
                self.ronald_ruth_4thorder_trajectory[i][1].append(self.orbiting_objects2D[i].position[1])

            orbiting_x_prev = [self.orbiting_objects2D[i].position[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_y_prev = [self.orbiting_objects2D[i].position[1] for i in range(len(self.orbiting_objects2D))]
            orbiting_vx_prev = [self.orbiting_objects2D[i].velocity[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_vy_prev = [self.orbiting_objects2D[i].velocity[1] for i in range(len(self.orbiting_objects2D))]
            orbiting_ax_prev = [self.orbiting_objects2D[i].acceleration[0] for i in range(len(self.orbiting_objects2D))]
            orbiting_ay_prev = [self.orbiting_objects2D[i].acceleration[1] for i in range(len(self.orbiting_objects2D))]
            
    def plot_euler_trajectory(self):
        for i in range(len(self.orbiting_objects2D)):
            plt.plot(self.euler_trajectory[i][0], self.euler_trajectory[i][1], color = COLORS[i % len(COLORS)], label = self.orbiting_objects2D[i].name)    
        plt.legend(loc = 'best')
        plt.show()

    def plot_verlet_trajectory(self):
        for i in range(len(self.orbiting_objects2D)):
            plt.plot(self.verlet_trajectory[i][0], self.verlet_trajectory[i][1], color = COLORS[i % len(COLORS)], label = self.orbiting_objects2D[i].name)    
        plt.legend(loc = 'best')
        plt.show()

    def plot_ronald_ruth_3rdorder_trajectory(self):
        for i in range(len(self.orbiting_objects2D)):
            plt.plot(self.ronald_ruth_3rdorder_trajectory[i][0], self.ronald_ruth_3rdorder_trajectory[i][1], color = COLORS[i % len(COLORS)], label = self.orbiting_objects2D[i].name)    
        plt.legend(loc = 'best')
        plt.show()

    def plot_ronald_ruth_4thorder_trajectory(self):
        for i in range(len(self.orbiting_objects2D)):
            plt.plot(self.ronald_ruth_4thorder_trajectory[i][0], self.ronald_ruth_4thorder_trajectory[i][1], color = COLORS[i % len(COLORS)], label = self.orbiting_objects2D[i].name)    
        plt.legend(loc = 'best')
        plt.show()

    def visualize_trajectory(self):
        pygame.init()
        win = pygame.display.set_mode((720,720))
        radius=7.5
        time_stamp=0
        length=len(self.euler_trajectory[0][0])
        run=True
        pygame.draw.circle(win,(255,255,0),(360,360),radius)
        scaling=1

        for i in range(len(self.orbiting_objects2D)):
            scaling=max(scaling,DISTANCE["SUN_"+self.orbiting_objects2D[i].name])

        while run:
            pygame.time.delay(10)

            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    run = False
            
            for i in range(len(self.orbiting_objects2D)):
                pygame.draw.circle(win, COLORS_PYGAME[i%len(COLORS_PYGAME)], (self.euler_trajectory[i][0][time_stamp]*330/scaling + 360 , self.euler_trajectory[i][1][time_stamp]*330/scaling + 360), radius/3)
                
            pygame.display.update()
            if time_stamp<length-2:
                time_stamp+=1

        pygame.quit()
# Some Predefined Objects

OBJECTS = {
    "SUN": Object2D('SUN', MASS["SUN"],[0,0],[0,0]),
    "MERCURY": Object2D('MERCURY', MASS["MERCURY"],[DISTANCE["SUN_MERCURY"],0],[0, INITIAL_VELOCITY["MERCURY"]]),
    "VENUS": Object2D('VENUS', MASS["VENUS"],[DISTANCE["SUN_VENUS"],0],[0, INITIAL_VELOCITY["VENUS"]]),
    "EARTH": Object2D('EARTH', MASS["EARTH"],[DISTANCE["SUN_EARTH"],0],[0, INITIAL_VELOCITY["EARTH"]]),
    "MARS": Object2D('MARS', MASS["MARS"],[DISTANCE["SUN_MARS"],0],[0, INITIAL_VELOCITY["MARS"]]),
    "JUPITER": Object2D('JUPITER', MASS["JUPITER"],[DISTANCE["SUN_JUPITER"],0],[0, INITIAL_VELOCITY["JUPITER"]]),
    "SATURN": Object2D('SATURN',MASS["SATURN"],[DISTANCE["SUN_SATURN"], 0],[0, INITIAL_VELOCITY["SATURN"]]),
    "URANUS": Object2D('URANUS', MASS["URANUS"],[DISTANCE["SUN_URANUS"],0],[0, INITIAL_VELOCITY["URANUS"]]),
    "NEPTUNE": Object2D('NEPTUNE', MASS["NEPTUNE"],[DISTANCE["SUN_NEPTUNE"],0],[0, INITIAL_VELOCITY["NEPTUNE"]]),
    "PLUTO": Object2D('PLUTO', MASS["PLUTO"],[DISTANCE["SUN_PLUTO"],0],[0, INITIAL_VELOCITY["PLUTO"]]),
    "STAR1": Object2D('EARTH', MASS["EARTH"],[DISTANCE["SUN_EARTH"],0],[0, INITIAL_VELOCITY["EARTH"]]),
    "STAR2": Object2D('EARTH', MASS["EARTH"],[1 + DISTANCE["SUN_EARTH"],0],[0, INITIAL_VELOCITY["EARTH"]]),
}