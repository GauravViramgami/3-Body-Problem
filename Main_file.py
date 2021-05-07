from Important_classes import *
from Global_constants import *
import matplotlib.pyplot as plt


# Sun_earth=TwoBodySystem(OBJECTS["SUN"],OBJECTS["EARTH"])
# Sun_earth.ronald_ruth_4thorder_method(0.002,100000)
# Sun_earth.plot_ronald_ruth_4thorder_trajectory()

# sun_earth_jupiter_saturn = NBodySystem(OBJECTS["SUN"], [OBJECTS["MERCURY"], OBJECTS["VENUS"], OBJECTS["EARTH"], OBJECTS["MARS"], OBJECTS["JUPITER"], OBJECTS["SATURN"], OBJECTS["URANUS"], OBJECTS["NEPTUNE"], OBJECTS["PLUTO"]])
# sun_earth_jupiter_saturn.verlet_method(0.02, 10000)
# sun_earth_jupiter_saturn.plot_verlet_trajectory()

# new = NBodySystem(Object2D("EARTH", MASS["EARTH"], [0,0], [0,0]), [Object2D("SUN", MASS["SUN"], [3/2, 0], [0, INITIAL_VELOCITY["EARTH"]]), Object2D("SUN", MASS["SUN"], [-3/2, 0], [0, -INITIAL_VELOCITY["EARTH"]])])
# new.verlet_method(0.02, 100)
# new.visualize_trajectory()
# stepsize = 0.02
# num_iterations = 10000
# case3 = NBodySystem(OBJECTS["SUN"], [OBJECTS["MERCURY"], OBJECTS["VENUS"], OBJECTS["EARTH"], OBJECTS["MARS"], OBJECTS["JUPITER"], OBJECTS["SATURN"], OBJECTS["URANUS"], OBJECTS["NEPTUNE"], OBJECTS["PLUTO"]])
# case3.euler_cromer_method(stepsize, num_iterations)
# case3.plot_euler_cromer_trajectory()

# case3.verlet_method(stepsize, num_iterations)
# case3.plot_verlet_trajectory()
# stepsize = 0.002
# num_iterations = 1000
# case4 = ThreeBodySystem(OBJECTS["SUN"], OBJECTS["EARTH"], OBJECTS["JUPITER_1000"])
# case4.euler_cromer_method(stepsize, num_iterations)
# case4.plot_euler_cromer_trajectory()

# case4.verlet_method(stepsize, num_iterations)
# case4.plot_verlet_trajectory()

stepsize = 0.002
num_iterations = 100000
case4 = NBodySystem(OBJECTS["SUN"], [OBJECTS["MARS"], OBJECTS["ASTEROID1"], OBJECTS["ASTEROID2"], OBJECTS["JUPITER"]])
case4.verlet_method(stepsize, num_iterations)
case4.plot_verlet_trajectory()