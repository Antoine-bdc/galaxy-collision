# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import io
import pickle


suffix = "col"
data_folder = "../data/data" + "_" + suffix
save_film = 0
grid_multiplier = 2

with io.open(data_folder + "/parameters.txt", mode="r", encoding="utf-8") as f:
    for line in f:
        table = line.split()
        if table[0] != "#":
            exec(str(table[0]) + " = " + str(table[1]))
            print(table[0], table[1])

n_x *= grid_multiplier
n_y *= grid_multiplier
n_z *= grid_multiplier

dx = L_x / n_x
dy = L_y / n_y
dz = L_z / n_z

total_frames = 0


def get_density(r_stars, m_stars):
    density = np.zeros((n_x, n_y, n_z))
    for l in range(n_stars):
        i, j, k = get_coordinates(r_stars[l])
        density[i][j][k] += m_stars[l] / (dx * dy * dz)
    return density


def get_coordinates(r_star):
    i = int(r_star[0] / dx)
    j = int(r_star[1] / dy)
    k = int(r_star[2] / dz)
    return i % n_x, j % n_y, k % n_z


def generate_data():
    global total_frames
    data = []
    for i in range(0, 1):
        for j in range(0, 10):
            for k in range(0, 10):
                for l in range(0, 1):
                    DATA = np.loadtxt(data_folder + "/time_00" + str(i) + str(j) + str(k) + str(l) + ".txt").T
                    POS = DATA[[0, 1, 2]].T
                    M_stars = DATA[9]
                    density_2D = np.zeros((n_x, n_y))
                    density_3D = get_density(POS, M_stars)
                    for x in range(n_x):
                        for y in range(n_y):
                            for z in range(n_z):
                                density_2D[x, y] += density_3D[x, y, z]
                    for x in range(n_x):
                        for y in range(n_y):
                            if density_2D[x, y] > 1e-20:
                                density_2D[x, y] = np.log(density_2D[x, y])
                            else:
                                density_2D[x, y] = 1
                    min = np.min(density_2D)
                    for x in range(n_x):
                        for y in range(n_y):
                            if density_2D[x, y] == 1:
                                density_2D[x, y] = min - 1
                    data.append(density_2D)
                    total_frames += 1
                    print("Loading data :", total_frames)
    print(total_frames)
    return data


DATA = generate_data()
print("Data generated")

#pickle.dump(DATA, open( "DATA3000.p", "wb" ) )
#DATA = pickle.load( open( "DATA3000.p", "rb" ) )
#total_frames = 300
vmin, vmax = np.min(DATA[0]), np.max(DATA[0])
fig = plt.figure(figsize=(13, 13))


def f(data):
    global total_frames, iteration
    iteration = (iteration + 1) % total_frames
    print(iteration)
    return data[iteration]


iteration = 0

im = plt.imshow(f(DATA), animated=True)


def updatefig(*args):
    global iteration
    im = plt.imshow(f(DATA), animated=True, vmin=vmin, vmax=vmax, label=str(iteration))
    return im,


plt.legend()
plt.colorbar()
ani = animation.FuncAnimation(fig, updatefig, interval=(100), frames=total_frames, blit=True)
if save_film:
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=1800, extra_args=['-vcodec', 'libx264'])
    ani.save('../collision_animation.mp4', writer=writer)



plt.show()
