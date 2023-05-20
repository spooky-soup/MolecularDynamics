import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Number of particles.
n = 1000
# timestamp
FPS = 1
dt = 1 / FPS
# Particle masses, scaled by some factor we're not using yet.
m = 1

# в данном случае рассмотрим молекулы гелия - He
sigma = 2.63
eps = 6.03
sigma_six = sigma ** 6


def distance(at1: 'Atom', at2: 'Atom') -> float:
    return abs((((at1.pos_x - at2.pos_x) ** 2) - (at1.pos_y - at2.pos_y) ** 2) ** 0.5)


def update_u(j, k) -> tuple[float, float]:
    rel_x = X[j] - X[k]
    rel_y = Y[j] - Y[k]
    r = (rel_x ** 2 + rel_y ** 2) ** 0.5
    r_six = r ** 6
    # Парный потенциал Леннард-Джонса
    pair_u = 4e3 * eps * (((sigma_six / r_six) ** 2) - (sigma_six / r_six))
    # Проекции потенциальной силы на оси x и y
    pair_u_x = float((rel_y * pair_u) / r)
    pair_u_y = float((rel_x * pair_u) / r)
    #аккумулируем силы для каждого атома
    uX[j] += pair_u_x
    uY[j] += pair_u_y
    uX[k] += -pair_u_x
    uX[j] += -pair_u_y
    return pair_u_x, pair_u_y


def calculate_acc(j) -> None:
    accX[j] = -(1 / m) * uX[j]
    accY[j] = -(1 / m) * uY[j]


def update_pos(j) -> None:
    new_pos_x = 2 * X[j] - X_prev[j] + accX[j] * dt * dt
    new_pos_y = 2 * Y[j] - Y_prev[j] + accY[j] * dt * dt
    if new_pos_x < 0 or new_pos_x > 1000:
        accX[j] = -accX[j]
        new_pos_x = 2 * X[j] - X_prev[j] + accX[j] * dt * dt
    if new_pos_y < 0 or new_pos_y > 1000:
        accY[j] = -accY[j]
        new_pos_y = 2 * Y[j] - Y_prev[j] + accY[j] * dt * dt
    X_prev[j] = X[j]
    Y_prev[j] = Y[j]
    X[j] = new_pos_x
    Y[j] = new_pos_y


uX = [0.0]*n
uY = [0.0]*n
accX = [0.0]*n
accY = [0.0]*n

X_prev = [random.random()*1000 for _ in range(n)]
Y_prev = [random.random()*1000 for _ in range(n)]

X = [X_prev[i] + 0.5 * accX[i] * dt * dt for i in range(n)]
Y = [Y_prev[i] + 0.5 * accY[i] * dt * dt for i in range(n)]


def advance():
    global uX, uY
    """Advance the simulation by dt seconds."""
    print("Advance")
    # Update the particles' potentials and acc according to their velocities.
    uX = [0.0]*n
    uY = [0.0]*n
    for j in range(n):
        for k in range(j + 1, n):
            if abs((((X[j] - X[k]) ** 2) - (Y[j] - Y[k]) ** 2) ** 0.5) < 2.5 * sigma:
                continue
            update_u(j, k)

        calculate_acc(j)
        update_pos(j)
    print("x = \f, y = \f", X[0], Y[0])


DPI = 100
width, height = 1000, 500

plt.close('all')

# animation
fig = plt.figure(2)
ax = plt.axes(xlim=(0, 1000), ylim=(0, 1000))
scat = ax.scatter(X, Y, s=10, c=np.arange(n))


def init_anim():
    """Initialize the animation"""
    scat.set_offsets(np.array([X, Y]))
    ax.set_title('Time = ' + '0.0 nsec')
    return scat


def animate(i):
    """Advance the animation by one step and update the frame."""
    global X, Y
    advance()
    x_np = np.array(X)
    y_np = np.array(Y)
    data = np.hstack((x_np[:, np.newaxis], y_np[:, np.newaxis]))
    ax.set_title('Time = ' + str(np.round(dt*i, decimals=2)) + ' nsec')
    scat.set_offsets(data)
    return scat


# Number of frames; set to None to run until explicitly quit.
frames = None
anim = FuncAnimation(fig, animate, frames=frames, interval=10, blit=False,
                     init_func=init_anim)

plt.show()
