## Introduction :

This code was written to execute celestial mechanics simulations with a high number of objects using particle mesh computational methods.
In particular it has been designed to simulate the collision of a galaxy with an impactor which can either be a point mass or another galaxy.
A few templates for constants.c with initial parameters running a typical galaxy collision simulation can be found in the /collisions folder.
The general structure of the code is the following :
- constants.c contains the global simulations parameters
- main.c initiates and executes simulations using parameters found in constants.c and functions found in physics.c. Additionally it writes simulation data in files
- physics.c contains all functions performing physical computations
- utils.c contains utility functions performing basic tasks


## Requirements :

GCC compiler to run the main program
Python libraries for plotting :
- numpy
- matplotlib
- astropy.constants
- mpl_toolkits
Doxygen to build the documentation


## Initial conditions files :

3 constants.c files can be found in the /collisions folder.
To use any of them one can just copy all their content and paste to the file named "constants.c" or just rename the file constants.c. Either way, only the file constants.c will be loaded by the simulation.
standard_simulation : runs a not too heavy standard simulation
static_plummer_galaxy : runs the (rather heavy) made to test all methods of the PM method
plummer_blackhole_collision : runs a simulation that will perform a collision between a Plummer galaxy and a SMBH.
plummer_plummer_high_b : runs a collision between two plummer galaxies with a high impact parameter
long_plummer_plummer : runs a long simulation of the collision between two plummer galaxies
## Run simulations :


- Set initial conditions and simulation parameters in the constants.c file. See documentation or in file comments for the description of each variable.
- compile using 'make' command
- execute using './main_exe' command


## Plot data :

To plot data execute the python file visualisation.py.
The program will then prompt a question for which plots you want to get, answer with the format [X1, X2, X3] where Xi is the data you want to plot.
It is also possible to type the command "python3 visualisation.py 2 3 5 6" to plot the data 2 3 5 6
By default the python program will plot the data from the folders ..\data\, ..\potential\  and ..\density\
If one wants to plot files from different data they can open the "visualisation.py" file and change the suffix. For instance suffix="plummer" will plot the data from the file ..\data_plummer\, ..\potential_plummer\  and ..\density_plummer\

Similarly, to plot an animation one can use the command "python3 animate_density.py". You can set the "save_film" variable to True in the file to save the animation as ../collision_animation.mp4


## Other :

-To clean all data created by the program run command 'make clean_data'
-The documentation can be found in HTML format in the folder ../html and can be accessed by opening index.html in any web browser
