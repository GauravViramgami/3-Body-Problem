# Global Constants

# Equation Constants (in SI)
PI = 3.141592
GRAVITATIONAL_CONSTANT = 6.674e-11

# Conversion of ((m^3)/(s^2)) to ((AU^3)/(yr^2))
CONVERSION_CONSTANT = ((6.68455e-12)**3)/((3.17098e-8)**2) 

# Mass (in kg)
MASS = {
    "SUN": 1.989e30, 
    "MERCURY": 3.285e23,
    "VENUS": 4.9e24,
    "EARTH": 6.0e24,
    "MARS" : 6.42e23,
    "JUPITER": 1.9e27,
    "SATURN": 5.7e26,
    "URANUS": 8.8e25,
    "NEPTUNE": 1.03e26,
    "PLUTO": 1.462e22,
    "JUPITER_1000": 1.9e30,
    "ASTEROID1": 2.8e21,
    "ASTEROID2": 3.2e21,
    "MOON": 0.07346e24,
}

# Distance (in AU)
DISTANCE = { 
    "SUN_SUN": 1.00,
    "SUN_MERCURY": 0.39,
    "SUN_VENUS": 0.72,
    "SUN_EARTH": 1.00,
    "SUN_MARS": 1.52,
    "SUN_JUPITER": 5.20,
    "SUN_SATURN": 9.54,
    "SUN_URANUS": 19.8,
    "SUN_NEPTUNE": 30.06,
    "SUN_PLUTO": 39.53,
    "SUN_JUPITER_1000": 5.20,
    "SUN_ASTEROID1": 2.7,
    "SUN_ASTEROID2": 3.1,
    "EARTH_MOON": 0.00256955529,
    "SUN_MOON": 1.00257,
}

# Time Period (in yr)
TIME_PERIOD = {
    "MERCURY": 0.2410,
    "VENUS": 0.6164,
    "EARTH": 1.0000,
    "MARS" : 1.8808,
    "JUPITER": 11.8600,
    "SATURN": 29.4571,
    "URANUS": 84.0205,
    "NEPTUNE": 164.8000,
    "PLUTO": 247.9400,
    "JUPITER_1000": 11.8600,
    "ASTEROID1": 4,
    "ASTEROID2": 5,
    "MOON": 27.322/365,
}

# Initial Velocity (in AU/yr)
INITIAL_VELOCITY = {
    "MERCURY": 2 * PI * DISTANCE["SUN_MERCURY"] / TIME_PERIOD["MERCURY"],
    "VENUS": 2 * PI * DISTANCE["SUN_VENUS"] / TIME_PERIOD["VENUS"],
    "EARTH": 2 * PI * DISTANCE["SUN_EARTH"] / TIME_PERIOD["EARTH"],
    "MARS" : 2 * PI * DISTANCE["SUN_MARS"] / TIME_PERIOD["MARS"],
    "JUPITER": 2 * PI * DISTANCE["SUN_JUPITER"] / TIME_PERIOD["JUPITER"],
    "SATURN": 2 * PI * DISTANCE["SUN_SATURN"] / TIME_PERIOD["SATURN"],
    "URANUS": 2 * PI * DISTANCE["SUN_URANUS"] / TIME_PERIOD["URANUS"],
    "NEPTUNE": 2 * PI * DISTANCE["SUN_NEPTUNE"] / TIME_PERIOD["NEPTUNE"],
    "PLUTO": 2 * PI * DISTANCE["SUN_PLUTO"] / TIME_PERIOD["PLUTO"],
    "ASTEROID1": 2 * PI * DISTANCE["SUN_ASTEROID1"] / TIME_PERIOD["ASTEROID1"],
    "ASTEROID2": 2 * PI * DISTANCE["SUN_ASTEROID2"] / TIME_PERIOD["ASTEROID2"],
    "JUPITER_1000": 2 * PI * DISTANCE["SUN_JUPITER_1000"] / TIME_PERIOD["JUPITER_1000"],
    "MOON": (2 * PI * DISTANCE["SUN_EARTH"] / TIME_PERIOD["EARTH"]) + (2 * PI * DISTANCE["EARTH_MOON"] / TIME_PERIOD["MOON"]),
}

# Colors for MatplotLib
COLORS = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

#Color for pygame
COLORS_PYGAME = [(219,206,202),	(248,226,176),(140,177,222),(240,231,231),(227,110,75),(197,171,110),(79,208,231),(75,112,221),(255, 241, 213)]
