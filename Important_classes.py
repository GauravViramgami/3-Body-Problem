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
