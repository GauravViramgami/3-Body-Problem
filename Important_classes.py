from Global_constants import *

class Object2D:
    def __init__ (self, name, mass, x_initial, y_initial, vx_initial, vy_initial):
        self.name = name
        self.mass = mass
        self.x_initial = x_initial
        self.y_initial = y_initial
        self.vx_initial = vx_initial
        self.vy_initial = vy_initial
        self.x = x_initial
        self.y = y_initial
        self.vx = vx_initial
        self.vy = vy_initial
        
    def euclideanDistance (self, x, y):
        r = (((self.x - x)**2) + ((self.y - y)**2))**(1/2)
        return r
    
    def updatePosition (self, x_new, y_new):
        self.x = x_new
        self.y = y_new
        
    def updateVelocity (self, vx_new, vy_new):
        self.vx = vx_new
        self.vy = vy_new

class TwoBodySystem:
    def __init__ (self, reference_object2D, orbiting_object2D):
        self.reference_object2D = reference_object2D
        self.orbiting_object2D = orbiting_object2D
        self.euler_trajectory = [[], []]
        self.verlet_trajectory = [[], []]

    def euler_cromer_method(self, stepsize, num_iterations):
      self.euler_trajectory = [[self.orbiting_object2D.x_initial], [self.orbiting_object2D.y_initial]]
      GM = 4*(PI**2)
    
      x_prev, y_prev = self.orbiting_object2D.x_initial, self.orbiting_object2D.y_initial
      vx_prev, vy_prev = self.orbiting_object2D.vx_initial, self.orbiting_object2D.vy_initial

      iterations = 0
      while (iterations < num_iterations):
          iterations = iterations + 1
          r = self.orbiting_object2D.euclideanDistance(0, 0)

          self.orbiting_object2D.updateVelocity(vx_prev - ((GM * x_prev)/(r**3)) * stepsize, vy_prev - ((GM * y_prev)/(r**3)) * stepsize) 

          self.orbiting_object2D.updatePosition(x_prev + (self.orbiting_object2D.vx * stepsize), y_prev + (self.orbiting_object2D.vy * stepsize)) 

          self.euler_trajectory[0].append(self.orbiting_object2D.x)
          self.euler_trajectory[1].append(self.orbiting_object2D.y)

          x_prev, y_prev = (self.orbiting_object2D.x, self.orbiting_object2D.y)
          vx_prev, vy_prev = (self.orbiting_object2D.vx, self.orbiting_object2D.vy)

    