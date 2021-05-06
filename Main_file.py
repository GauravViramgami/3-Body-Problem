from Important_classes import *
from Global_constants import *
import matplotlib.pyplot as plt
sun= Object2D('SUN', MASS["SUN"],0,0,0,0)
saturn= Object2D('SATURN',MASS["SATURN"],DISTANCE["SUN_SATURN"],0,0,2*PI*9.537/29.46)
sun_saturn_system = TwoBodySystem(sun,saturn)
sun_saturn_system.euler_cromer_method(0.002,100000)
# print(sun_saturn_system.verlet_trajectory[0][:1000])
plt.plot(sun_saturn_system.euler_trajectory[0],sun_saturn_system.euler_trajectory[1])
plt.show()
# for i in range(100):
#     print(sun_saturn_system.verlet_trajectory[0][i], sun_saturn_system.verlet_trajectory[1][i])
# print(sun_saturn_system.orbiting_object2D.ax)
# print(sun_saturn_system.force_on_orbiting)
# print((GRAVITATIONAL_CONSTANT * MASS[sun_saturn_system.reference_object2D.name] * MASS[sun_saturn_system.orbiting_object2D.name] * sun_saturn_system.orbiting_object2D.x) / ((sun_saturn_system.orbiting_object2D.euclideanDistance(0, 0))**3))
# print(sun_saturn_system.orbiting_object2D.x)