from Important_classes import *
from Global_constants import *
import matplotlib.pyplot as plt

sun = Object2D('SUN', MASS["SUN"],[0,0],[0,0])
earth = Object2D('EARTH', MASS["EARTH"],[DISTANCE["SUN_EARTH"],0],[0, INITIAL_VELOCITY["EARTH"]])
jupiter = Object2D('JUPITER', MASS["JUPITER"],[DISTANCE["SUN_JUPITER"],0],[0, INITIAL_VELOCITY["JUPITER"]])
saturn = Object2D('SATURN',MASS["SATURN"],[DISTANCE["SUN_SATURN"], 0],[0, INITIAL_VELOCITY["SATURN"]])
# sun_saturn_system = TwoBodySystem(sun,saturn)
# sun_saturn_system.euler_cromer_method(0.02,100000)
# sun_saturn_system.plot_euler_trajectory()
# sun_earth_jupiter = ThreeBodySystem(sun, earth, jupiter)
# sun_earth_jupiter.euler_cromer_method(0.002, 1000)
# sun_earth_jupiter.plot_euler_trajectory()
sun_saturn = NBodySystem(sun, [earth, jupiter, saturn])
sun_saturn.euler_cromer_method(0.02, 100000)
sun_saturn.plot_euler_trajectory()