from Important_classes import *
from Global_constants import *
import matplotlib.pyplot as plt


# Sun_earth=TwoBodySystem(OBJECTS["SUN"],OBJECTS["EARTH"])
# Sun_earth.ronald_ruth_4thorder_method(0.002,100000)
# Sun_earth.plot_ronald_ruth_4thorder_trajectory()

sun_earth_jupiter_saturn = NBodySystem(OBJECTS["SUN"], [OBJECTS["MERCURY"], OBJECTS["VENUS"], OBJECTS["EARTH"], OBJECTS["MARS"], OBJECTS["JUPITER"], OBJECTS["SATURN"], OBJECTS["URANUS"], OBJECTS["NEPTUNE"], OBJECTS["PLUTO"]])
sun_earth_jupiter_saturn.ronald_ruth_4thorder_method(0.02, 10000)
sun_earth_jupiter_saturn.visualize_trajectory()

# new = NBodySystem(OBJECTS["SUN"], [OBJECTS["STAR1"], OBJECTS["STAR2"]])
# new.verlet_method(0.02, 100)
# new.plot_verlet_trajectory()