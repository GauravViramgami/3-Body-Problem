from Important_classes import *
from Global_constants import *
import matplotlib.pyplot as plt


Sun_earth=ThreeBodySystem(OBJECTS["SUN"],OBJECTS["EARTH"],OBJECTS["JUPITER"])
Sun_earth.euler_cromer_method(0.002,100000)
Sun_earth.visualize_trajectory()

# sun_earth_jupiter_saturn = NBodySystem(OBJECTS["SUN"], [OBJECTS["MERCURY"], OBJECTS["VENUS"], OBJECTS["EARTH"], OBJECTS["MARS"], OBJECTS["JUPITER"], OBJECTS["SATURN"], OBJECTS["URANUS"], OBJECTS["NEPTUNE"], OBJECTS["PLUTO"]])
# sun_earth_jupiter_saturn.ronald_ruth_3rdorder_method(0.02, 10000)
# sun_earth_jupiter_saturn.plot_ronald_ruth_3rdorder_trajectory()

# new = NBodySystem(OBJECTS["SUN"], [OBJECTS["STAR1"], OBJECTS["STAR2"]])
# new.verlet_method(0.02, 100)
# new.plot_verlet_trajectory()