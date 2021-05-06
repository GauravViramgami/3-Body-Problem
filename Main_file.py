from Important_classes import *
from Global_constants import *
import matplotlib.pyplot as plt

sun= Object2D('SUN', MASS["SUN"],[0,0],[0,0])
saturn= Object2D('SATURN',MASS["SATURN"],[DISTANCE["SUN_SATURN"],0],[0,2*PI*9.537/29.46])
sun_saturn_system = TwoBodySystem(sun,saturn)
sun_saturn_system.euler_cromer_method(0.02,100000)
sun_saturn_system.plot_euler_trajectory()