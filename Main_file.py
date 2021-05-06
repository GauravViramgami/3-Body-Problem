from Important_classes import *
from Global_constants import *
import matplotlib.pyplot as plt
sun= Object2D('Sun', MASS_SUN,0,0,0,0)
saturn= Object2D('Saturn',MASS_SATURN,DISTANCE_SUN_SATURN,0,0,2*PI*9.537/29.46)
sun_saturn_system = TwoBodySystem(sun,saturn)
sun_saturn_system.euler_cromer_method(0.002,100000)
plt.plot(sun_saturn_system.euler_trajectory[0],sun_saturn_system.euler_trajectory[1])
plt.show()