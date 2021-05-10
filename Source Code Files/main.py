from Important_classes import *
from Global_constants import *
import matplotlib.pyplot as plt

## Case:1
# stepsize = 0.002
# num_iterations = 100000
# case1 = TwoBodySystem(OBJECTS["SUN"], OBJECTS["EARTH"])
# case1.euler_method(stepsize, num_iterations)
# case1.plot_euler_trajectory()

# case1.euler_cromer_method(stepsize, num_iterations)
# case1.plot_euler_cromer_trajectory()

# case1.verlet_method(stepsize, num_iterations)
# case1.plot_verlet_trajectory()

## Case:2
# stepsize = 0.002
# num_iterations = 100000
# case2 = ThreeBodySystem(OBJECTS["SUN"], OBJECTS["MARS"], OBJECTS["JUPITER"])

# case2.euler_method(stepsize, num_iterations)
# case2.plot_euler_trajectory()

# case2.euler_cromer_method(stepsize, num_iterations)
# case2.plot_euler_cromer_trajectory()

# case2.verlet_method(stepsize, num_iterations)
# case2.plot_verlet_trajectory()

## Case:3
# stepsize = 0.02
# num_iterations = 10000
# case3 = NBodySystem(OBJECTS["SUN"], [OBJECTS["MERCURY"], OBJECTS["VENUS"], OBJECTS["EARTH"], OBJECTS["MARS"], OBJECTS["JUPITER"], OBJECTS["SATURN"], OBJECTS["URANUS"], OBJECTS["NEPTUNE"], OBJECTS["PLUTO"]])

# case3.euler_method(stepsize, num_iterations)
# case3.plot_euler_trajectory()

# case3.euler_cromer_method(stepsize, num_iterations)
# case3.plot_euler_cromer_trajectory()

# case3.verlet_method(stepsize, num_iterations)
# case3.plot_verlet_trajectory()

## Case:4
# stepsize = 0.002
# num_iterations = 1000
# case4 = ThreeBodySystem(OBJECTS["SUN"], OBJECTS["MARS"], OBJECTS["JUPITER_1000"])

# case4.euler_method(stepsize, num_iterations)
# case4.plot_euler_trajectory()

# case4.euler_cromer_method(stepsize, num_iterations)
# case4.plot_euler_cromer_trajectory()

# case4.verlet_method(stepsize, num_iterations)
# case4.plot_verlet_trajectory()

## Case:5
stepsize = 0.002
num_iterations = 100000
case5 = NBodySystem(OBJECTS["SUN"], [OBJECTS["MARS"], OBJECTS["ASTEROID1"], OBJECTS["ASTEROID2"], OBJECTS["JUPITER"]])

# case5.euler_method(stepsize, num_iterations)
# case5.plot_euler_trajectory()

# case5.euler_cromer_method(stepsize, num_iterations)
# case5.plot_euler_cromer_trajectory()

case5.verlet_method(stepsize, num_iterations)
# case5.plot_verlet_trajectory()
case5.visualize_trajectory()

## Case:6
# stepsize = 0.002
# num_iterations = 100000
# case6 = ThreeBodySystem(OBJECTS["SUN"], OBJECTS["EARTH"], OBJECTS["MOON"])

# case6.euler_method(stepsize, num_iterations)
# case6.plot_euler_trajectory()

# case6.euler_cromer_method(stepsize, num_iterations)
# case6.plot_euler_cromer_trajectory()

# case6.verlet_method(stepsize, num_iterations)
# case6.plot_verlet_trajectory()
