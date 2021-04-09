import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as c
import io
import time
from mpl_toolkits.mplot3d import Axes3D
import sys

suffix = "col"

data_folder = "../data/data" + "_" + suffix
potential_folder = "../data/potential"  + "_" + suffix
density_folder = "../data/density"  + "_" + suffix

# Importing simulation parameters
with io.open(data_folder + "/parameters.txt", mode="r", encoding="utf-8") as f:
    for line in f:
        table = line.split()
        if table[0] == "n_t":
            n_t = table[1]
        if table[0] != "#":
            exec(str(table[0]) + " = " + str(table[1]))
            #print(table[0], table[1])

# Additionnal variables
center = L_x / 2
dx = L_x / n_x
dy = L_y / n_y
dz = L_z / n_z
M = n_stars * m
G = c.G.value
PI = np.pi


R_max = (((n_t * dt)**2 * G * M / (4 * PI**2))**(2. / 3) - a1**2)**(1. / 2)
print("R_max =", R_max)
print("R_max/L_x =", R_max / L_x)
print("R_max/a =", R_max / a1)


numbers = []
for i in range(-100, 100):
    numbers.append(str(i))

plot_list = []
for i in range(1, len(sys.argv)):
    if sys.argv[i] in numbers:
        plot_list.append(int(sys.argv[i]))
    else:
        plot_list.append(sys.argv[i])


def norm(A):
    return (A[0]**2 + A[1]**2 + A[2]**2)**(0.5)


if len(plot_list) == 0:
    char_list = input("Which data would you like to plot ?   (Format [1,4,5] for data corresponding to 1, 4 and 5)\n\nList of data :\n7 - Evolution of time over all time steps \n6 - Measured density vs theoretical\n5 - Measured potential vs theoretical\n4 - Star distribution in 2D projection at time 0\n3 - Energy evolution over time\n2 - Measured acceleration, speed and density vs theoretical\n1 - Trajectories of all stars in one plot\n")
    print("here", type(char_list))
    exec("plot_list = " + char_list)
    print("here", plot_list, type(plot_list))
will_plot = []
if (7 in plot_list or 'time' in plot_list):
    will_plot.append("- Evolution of time over all time steps (only for dynamical dt)")
if (6 in plot_list or 'density' in plot_list):
    will_plot.append("- Measured density vs theoretical")
if (5 in plot_list or 'potential' in plot_list):
    will_plot.append("- Measured potential vs theoretical")
if (3 in plot_list or 'energy' in plot_list):
    will_plot.append("- Energy evolution over time")
if (2 in plot_list or 'acceleration' in plot_list):
    will_plot.append("- Measured acceleration, speed and density vs theoretical")
if (1 in plot_list or 'trajectory' in plot_list):
    will_plot.append("- Trajectories of all stars in one plot (small number of particles only)")
if (4 in plot_list or 'scatter' in plot_list):
    will_plot.append("- Star distribution in 2D projection at time 0")

print("The routine will plot :")
for i in range(len(will_plot)):
    print(will_plot[i])
print("")


def draw_plots(plot_list, n_t):
    # Evolution of dt
    if (7 in plot_list or 'time' in plot_list):
        dt_data = np.loadtxt(data_folder + "/dt.txt")
        print(np.shape(dt_data))
        plt.plot(range(n_t - 1), dt_data[1:, 1])
        plt.yscale("log")
        plt.show()
        plt.plot(dt_data[1:, 0], dt_data[1:, 1])
        plt.yscale("log")
        plt.show()
        plt.plot(range(n_t - 1), dt_data[1:, 0])
        #plt.yscale("log")
        plt.show()

    # Test du calcul de la densité
    if (6 in plot_list or 'density' in plot_list):
        for i in range(0, 31, 3):
            date = ""
            for j in range(6, 0, -1):
                if j == int(np.log10(n_t)):
                    date += str(i)

                else:
                    date += "0"
            if len(date) == 7:
                date2 = ""
                print("too long")
                for k in range(1, 7):
                    date2 += str(date[k])
                date = date2
            n_date = int(date)
            if n_date == 0:
                date = "000001"
                n_date = 1
            density = np.loadtxt(density_folder + "/density_" + date + ".txt")
            th = np.copy(density)
            for i in range(n_x):
                x = dx * (i - 1. / 2)
                for j in range(n_y):
                    y = dx * (j - 1. / 2)
                    k = n_x // 2
                    z = dx * (k - 1. / 2)
                    th[i][j] = 3 * M * a1**2 / (4 * np.pi * ((x - center)**2 + (y - center)**2 + (z - center)**2 + a1**2)**(5. / 2))


            _min, _max = np.amin(th), np.amax(th)

            fig, (ax1, ax2, ax3) = plt.subplots(figsize=(13, 3), ncols=3)
            theorique = ax1.imshow(th, vmin=_min, vmax=_max)
            fig.colorbar(theorique, ax=ax1)
            ax1.set_title("Densité théorique")

            data = ax2.imshow(density, vmin=_min, vmax=_max)
            fig.colorbar(data, ax=ax2)
            ax2.set_title("Densité obtenue")

            difference = ax3.imshow(th - density, vmin=_min, vmax=_max)
            fig.colorbar(difference, ax=ax3)
            ax3.set_title("Différence (théorique - calculée)")
            plt.suptitle("Time iteration = " + str(n_date))
            plt.show()

            fig, (ax1) = plt.subplots(figsize=(9, 7), ncols=1)
            ax1.plot(th[n_x // 2, :], label="Densité théorique")
            ax1.plot(density[n_x // 2, :], label="Densité calculée")
            ax1.legend()
            ax1.set_title("Coupe de la densité à y = z = 1/2")
            for i in range(n_x):
                plt.axvline(i, linewidth=0.1)
            plt.suptitle("Time iteration = " + str(n_date))
            plt.show()

            plt.imshow(density)
            plt.tight_layout()
            plt.title("Densité du profil de Plummer")
            plt.colorbar()
            plt.show()

            # Shell density
            fig, ax = plt.subplots()
            DATA = np.loadtxt(data_folder + "/time_"  + date + ".txt").T
            POS = DATA[[0, 1, 2]]
            VEL = DATA[[3, 4, 5]]
            N = n_stars

            center_x = np.mean(DATA[:int(fraction_stars * n_stars)][0])
            center_y = np.mean(DATA[:int(fraction_stars * n_stars)][1])
            center_z = np.mean(DATA[:int(fraction_stars * n_stars)][2])
            #center_x, center_y, center_z = L_x/2, L_y/2, L_z/2,
            print("center at :", center_x / L_x, center_y / L_y, center_z / L_z)

            for i in range(N):
                POS[0][i] -= center_x
                POS[1][i] -= center_y
                POS[2][i] -= center_z

            ACC = DATA[[6, 7, 8]]
            scalar_acc = np.zeros(N)
            scalar_vel = np.zeros(N)
            th_acc = np.zeros(N)
            radial_pos = np.zeros(N)

            for i in range(N):
                radial_pos[i] = norm(POS.T[i])

            # Plot radial density distribution
            n, bins, paches = plt.hist(radial_pos, bins=100, label="Radial distribution")
            r_dist = np.zeros(len(bins) - 1)
            X = range(len(n))
            rho_th = np.zeros(len(n))
            for i in range(len(n)):
                r_dist[i] = n[i] * m / (bins[i + 1]**3 * 4 * np.pi / 3 - bins[i]**3 * 4 * np.pi / 3)
                r = np.mean([bins[i + 1], bins[i]])
                rho_th[i] = (3 * M / (4 * np.pi)) * (a1**2 / (r**2 + a1**2)**(5. / 2))
            plt.clf()
            plt.plot(X, rho_th, label="Density théorique")
            plt.plot(X, r_dist, label="Densité obtenue sous forme de bins de distance radiale")
            #plt.xscale("log")
            plt.title("Densité en fonction du rayon r")
            plt.yscale("log")
            plt.xscale("log")
            plt.legend()
            plt.suptitle("Iteration : " + str(n_date))
            plt.show()



    # Testing calculation of the potential
    if (5 in plot_list or 'potential' in plot_list):
        for i in range(0, 31, 6):
            date = ""
            for j in range(6, 0, -1):
                if j == int(np.log10(n_t)):
                    date += str(i)
                else:
                    date += "0"
            if len(date) == 7:
                date2 = ""
                print("too long")
                for k in range(1, 7):
                    print(k, len(date))
                    date2 += str(date[k])
                date = date2
            n_date = int(date)
            if n_date == 0:
                date = "000001"
                n_date = 1
            potential = np.loadtxt(potential_folder + "/potential_" + date + ".txt")
            th = np.copy(potential)
            for i in range(n_x):
                x = dx * (i - 1. / 2)
                for j in range(n_y):
                    y = dx * (j - 1. / 2)
                    k = n_x // 2
                    z = dx * (k - 1. / 2)
                    th[i][j] = - M * G / ((x - center)**2 + (y - center)**2 + (z - center)**2 + a1**2)**(1. / 2)  # Th pour Plummer
                    #th[i][j] =  -G * M / ((x-center)**2 + (y-center)**2 + (z-center)**2)**(1./2)    # Th pour 1 corps
            diff = th - potential

            fig, (ax1, ax2, ax3) = plt.subplots(figsize=(13, 4), ncols=3)
            theorique = ax1.imshow(th, vmax=0)
            fig.colorbar(theorique, ax=ax1)
            ax1.set_title("Potentiel theorique")

            data = ax2.imshow(potential, vmax=0)
            fig.colorbar(data, ax=ax2)
            ax2.set_title("Potentiel calculé")

            difference = ax3.imshow(diff)
            fig.colorbar(difference, ax=ax3)
            ax3.set_title("Différence théorique - calculé")
            plt.suptitle("Time iteration = " + str(n_date))
            plt.show()


            plt.imshow(potential)
            plt.tight_layout()
            plt.title("Potentiel gravitationnel du profil de Plummer")
            plt.colorbar()
            plt.show()


            fig, (ax1) = plt.subplots(figsize=(9, 7), ncols=1)
            ax1.plot(-th[n_x // 2, :], label="potentiel théorique")
            ax1.plot(-potential[n_x // 2, :], label="potentiel calculé")
            ax1.legend()
            plt.yscale("log")
            ax1.set_title("Coupe du potentiel phi à z = y = 1/2")
            for i in range(n_x):
                plt.axvline(i, linewidth=0.1)
            plt.suptitle("Time iteration = " + str(n_date))
            plt.show()


    if (4 in plot_list or 'scatter' in plot_list):
        for i in range(0, 1):

            fig, ax = plt.subplots()
            DATA = np.loadtxt(data_folder + "/time_00000" + str(i) + ".txt").T
            plt.xlim(0, L_x)
            plt.ylim(0, L_x)
            plt.scatter(DATA[0, :], DATA[1, :], s=2, marker="o", color="blue")
            OLD_DATA = np.array(DATA)
            plt.show()

    # Energies
    if (3 in plot_list or 'energy' in plot_list):
        DATA = np.loadtxt(data_folder + "/energy.txt").T

        R = range(len(DATA[0]))
        DIV = np.copy(DATA[1, :])
        for i in range(len(DIV)):
            DIV[i] *= -2 / DATA[2, i]
        plt.plot(DATA[0, :], DATA[1, :], label="kinetic_energy")
        plt.plot(DATA[0, :], DATA[2, :], label="potential_energy")
        plt.plot(DATA[0, :], DATA[3, :], label="total_energy")
        plt.legend()
        #plt.yscale("log")
        plt.show()

        plt.plot(DATA[0, :], DIV, label="-2 E_k / E_p")
        plt.ylim(0, 2)
        plt.legend()
        plt.show()




        # Acceleration tests
    if (2 in plot_list or 'acceleration' in plot_list):
        for i in range(0, 10, 1):
            date = ""
            for j in range(6, 0, -1):
                if j == int(np.log10(n_t)):
                    date += str(i)

                else:
                    date += "0"
            if len(date) == 7:
                date2 = ""
                print("too long")
                for k in range(1, 7):
                    date2 += str(date[k])
                date = date2
            n_date = int(date)
            if n_date == 0:
                date = "000001"
                n_date = 1
            fig, ax = plt.subplots()
            DATA = np.loadtxt(data_folder + "/time_" + date + ".txt").T
            POS = DATA[[0, 1, 2]]
            VEL = DATA[[3, 4, 5]]
            N = n_stars

            center_x = np.mean(DATA[:int(fraction_stars * n_stars)][0])
            center_y = np.mean(DATA[:int(fraction_stars * n_stars)][1])
            center_z = np.mean(DATA[:int(fraction_stars * n_stars)][2])
            #center_x, center_y, center_z = L_x/2, L_y/2, L_z/2,
            print("center at :", center_x / L_x, center_y / L_y, center_z / L_z)

            for i in range(N):
                POS[0][i] -= center_x
                POS[1][i] -= center_y
                POS[2][i] -= center_z

            ACC = DATA[[6, 7, 8]]
            scalar_acc = np.zeros(N)
            scalar_vel = np.zeros(N)
            th_acc = np.zeros(N)
            radial_pos = np.zeros(N)

            for i in range(N):
                radial_pos[i] = norm(POS.T[i])
                ve = (2 * G * M)**(1. / 2) * (radial_pos[i]**2 + a1**2)**(-1. / 4)
                scalar_vel[i] = norm(VEL.T[i]) / ve
                scalar_acc[i] = norm(ACC.T[i])
                th_acc[i] = -(G * M / (radial_pos[i]**2)) / (1 + a1**2 / radial_pos[i]**2)**(3. / 2)
            tableau = np.zeros(100)
            tableau_th = np.zeros(100)
            for j in range(100):
                r_repart = 1 + j * L_x / 200
                for i in range(N):
                    if radial_pos[i] < r_repart:
                        tableau[j] += m
                tableau_th[j] = M / (1 + (a1 / r_repart)**2)**(3. / 2)
            plt.plot(np.linspace(0, L_x, 100), tableau, label="Distribution obtenue M(r)")
            plt.plot(np.linspace(0, L_x, 100), tableau_th, label="Distribution de Plummer M(r)")
            plt.legend()
            plt.title("Distrubution de M(r) : la masse contenue dans la sphere de rayon r")
            plt.suptitle("Iteration : " + str(n_date))
            plt.show()
            plt.xscale("log")
            plt.yscale("log")

            plt.plot(radial_pos, scalar_acc, ls='', marker='o', ms=1, label="|a|(r) obtenue via la méthode PM")
            plt.plot(radial_pos, -th_acc, ls='', marker='o', ms=1, label="|a|(r) de Plummer")

            for i in range(n_x):
                plt.axvline(i * dx, linewidth=0.1)
            plt.legend()
            plt.title("Norme de l'accélération en fonction de la distance au centre r")
            plt.suptitle("Iteration : " + str(n_date))
            plt.show()

            # Plot velocity distribution
            N = 10000
            Xp = np.linspace(0, 1, N)
            Z = np.zeros(N)

            for i in range(N):
                x = Xp[i] - 1 / (2 * N)
                Z[i] = (1 - x**2)**(7. / 2) * x**2 / 0.042951
            plt.plot(Xp[1:], Z[1:], label="Theoretical distribution")
            plt.hist(scalar_vel, bins=100, label="Velocity distribution", density=True)

            plt.title("Distribution de la norme des vitesses")
            plt.legend()
            plt.suptitle("Iteration" + str(n_date))
            plt.show()

        #Trajectories
    if (1 in plot_list or 'trajectory' in plot_list):
        size = n_stars / 50000
        print(size)
        plt.xlim(0, L_x)
        plt.ylim(0, L_x)
        n_t -= 1

        i1_max = int(n_t / 1000)
        if (i1_max != 0):
            i2_max, i3_max, i4_max = 9, 9, 9
        else:
            i2_max = int((n_t - i1_max * 1000) / 100)
            if (i2_max != 0):
                i3_max, i4_max = 9, 9
            else:
                i3_max = int((n_t - i1_max * 1000  - i2_max * 100) / 10)
                if (i3_max != 0):
                    i4_max = 9
                else:
                    i4_max = int((n_t - i1_max * 1000  - i2_max * 100 - i3_max * 10))
        print(i1_max, i2_max, i3_max, i4_max)


        #for i1 in range(0, 1 + i1_max):
        #    for i2 in range(0, 1 + i2_max):
        #        for i3 in range(0, 1 + i3_max):
        #            for i4 in range(0, 1 + i4_max):
        dim3 = 0
        fig = plt.figure(figsize=(13, 13))
        col_tab = np.zeros((n_stars, 3))
        for i in range(n_stars):
            col_tab[i][0] = np.random.random()
            col_tab[i][1] = np.random.random()
            col_tab[i][2] = np.random.random()
        for i1 in range(0, 1):
            for i2 in range(0, 1):
                for i3 in range(0, 1):
                    for i4 in range(0, 1):
                        t = 1000 * i1 + 100 * i2 + 10 * i3 + i4
                        print(str(i1) + str(i2) + str(i3) + str(i4))
                        PM_DATA = np.loadtxt(data_folder + "/time_00" + str(i1) + str(i2) + str(i3) + str(i4) + ".txt").T
                        #PP_DATA = np.loadtxt(data_folder + "/time_00"+str(i1)+str(i2)+str(i3)+str(i4)+".txt").T
                        if not dim3:
                            for k in range(n_stars):
                                if (str(i1) + str(i2) + str(i3) + str(i4) == "0001"):
                                    plt.scatter(PM_DATA[0][k], PM_DATA[1][k], s=size, marker=".", color=col_tab[k], label="PM")
                                #    plt.scatter(PP_DATA[0,:], PP_DATA[1,:], s=size, marker="o", color = "blue", label="PP") #color = [(t / n_t, 1 - t/n_t, 1)]
                                else:
                                    plt.scatter(PM_DATA[0][k], PM_DATA[1][k], s=size, marker=".", color=col_tab[k])
                                #    plt.scatter(PP_DATA[0,:], PP_DATA[1,:], s=size, marker="o", color = "blue") #color = [(t / n_t, 1 - t/n_t, 1)]
                            for i in range(n_x):
                                plt.axvline(i * dx, linewidth=0.1)
                                plt.axhline(i * dx, linewidth=0.1)
                        else:
                            ax = fig.add_subplot(111, projection='3d')
                            for k in range(n_stars):
                                ax.scatter(PM_DATA[0][k], PM_DATA[1][k], PM_DATA[2][k], marker='.', s=1, color=col_tab[k])
                            ax.set_xlabel('X Label')
                            ax.set_ylabel('Y Label')
                            ax.set_zlabel('Z Label')
        plt.show()
        n_t += 1


draw_plots(plot_list, n_t)
