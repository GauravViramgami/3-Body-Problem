import matplotlib.pyplot as plt

pi = 3.14
INF = 1e400

def eucliedienDistance (point1, point2):
    r = ((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2)**(1/2)
    return r

def sunSaturn (position, velocity, h, iterations):
    coordinates = [[position[0]], [position[1]]]
    GM = 4*(pi**2)
    
    x_prev, y_prev = position
    vx_prev, vy_prev = velocity
    x_new, y_new = position
    vx_new, vy_new = velocity

    numIterations = 0
    while (numIterations < iterations):
        numIterations = numIterations + 1
        r = eucliedienDistance((0, 0), (x_prev, y_prev))

        vx_new = vx_prev - ((GM * x_prev)/(r**3)) * h
        vy_new = vy_prev - ((GM * y_prev)/(r**3)) * h

        x_new = x_prev + (vx_new * h)
        y_new = y_prev + (vy_new * h)

        coordinates[0].append(x_new)
        coordinates[1].append(y_new)

        x_prev, y_prev = (x_new, y_new)
        vx_prev, vy_prev = (vx_new, vy_new)

    return coordinates

coordinates = sunSaturn((9.537, 0), (0, 2*pi*9.537/29.46), 0.002, 100000)
plt.plot(coordinates[0], coordinates[1])
plt.show()