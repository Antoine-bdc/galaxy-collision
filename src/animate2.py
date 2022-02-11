# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import io
from time import sleep

suffix = ""
data_folder = "../data/kinematics"


with io.open(data_folder + "/parameters.txt", mode="r", encoding="utf-8") as f:
    for line in f:
        table = line.split()
        if table[0] != "#":
            exec(str(table[0]) + " = " + str(table[1]))
            print(table[0], table[1])


data_names = []
for i in range(0, 1):
    for j in range(0, 10):
        for k in range(0, 10):
            for l in range(0, 10):
                data_names.append(data_folder + "/kinematicstime_00" + str(i) + str(j) + str(k) + str(l) + ".txt")

fig = plt.figure()
fig.subplots_adjust(left=0, bottom=0, right=1, top=1)
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor((0.95, 0.95, 0.95))


colors = []
for i in range(5000):
    if i < 5000 * 2 / 3:
        colors.append('blue')
    if i >= 5000 * 2 / 3:
        colors.append('red')

# for i, name in enumerate(data_names):
#
#     axes.append(ax)

plt.axis('on')
ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.set_zlim(0, 2)
fig.set_size_inches(8, 8)


def update(i, fig, ax):
    DATA = np.loadtxt(data_names[i]).T
    X, Y, Z = DATA[[0, 1, 2]] / 1e17
    print(f"{round(i / len(data_names) * 100, 2)} % - ", X[0], Y[0], Z[0])

    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor((0.95, 0.95, 0.95))
    ax.scatter(X, Y, Z, alpha=0.7, s=100 / 2**(np.log2(len(X))), marker=',', linewidth=1, c=colors)

    ax.view_init(elev=110, azim=0)

    return fig, ax


anim = FuncAnimation(fig, update, frames=np.arange(0, len(data_names)), repeat=False, fargs=(fig, ax))

Writer = animation.writers['ffmpeg']
writer = Writer(fps=24, bitrate=2600, extra_args=['-vcodec', 'libx264'])
anim.save('../light_version/animations/scatter_plot.mp4', writer=writer)
