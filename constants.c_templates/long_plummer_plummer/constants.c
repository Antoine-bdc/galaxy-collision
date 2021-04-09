/**
* @file constants.c
* @brief File containing all initial conditions and simulation parameters
*/
#include <stdlib.h>

// Numerical constants:
double L_x = 2.e17;                      /// X-dimension physical length of the grid
double L_y = 2.e17;                      /// Y-dimension physical length of the grid
double L_z = 2.e17;                      /// Z-dimension physical length of the grid
double n_x = 64;                         /// Grid number of X-cells
double n_y = 64;                         /// Grid number of Y-cells
double n_z = 64;                         /// Grid number of Z-cells
double C_vel = 1.e-2;                    /// Defines precision required for dynamical dt in terms of particles velocity
double C_dyn = 1.e-2;                    /// Defines precision required for dynamical dt in terms of grid density
double t = 0.;                           /// Initial time

// Physical constants:
double G = 6.67e-11;                     /// Universal gravity constant
double M_tot = 0;                        /// Total mass in case of central point mass acceleration only

// Methods :
int acceleration_method = 2;             /// Acceleration method (0 = Central), (1 = PP), (2 = PM),
int differential_method = 2;             /// Differential method (0 = EULER), (1 = LEAPFROG), (2 = KICK_DRIFT_KICK)
int potential_method = 0;                /// Potential calculation method (0 = Jacobi), (1 = Gauss_Seidel), (2 = SOR)
int dynamical_dt = 0;                    /// Dynamical evolution of time step (0 = False), (1 = True)
int second_object = 1;                   /// Presence of additionnal object (0 = None), (1 = galaxy 2), (2 = SMBH)

// Simulation parameters:
double dt = 1e10;                        /// Time step delta t
int n_t = 10001  ;                        /// Total number of time steps
double jacobi_precision = 1.e-4;         /// Potential precision for jacobi method
double gauss_seidel_precision = 1.e-4;   /// Potential precision for Gauss Seidel and SOR method

// initial conditions
double m = 1.5e31;                       /// Mass for each star
double m_SMBH = 5.e34;                  /// Mass of the supermassive black hole if applicable
int n_stars = 5000;                         /// Total number of stars
double fraction_stars = 0.66;              /// Fraction of stars in the galaxy 1
double a1 = 0.05;                         /// Galaxy 1 typical size (in grid length fraction)
double a2 = 0.025;                        /// Galaxy 2 typical size (in grid length fraction) if applicable

double x_1 = 1./2;                       /// initial position x of galaxy 1 (in grid length fraction)
double y_1 = 1./2;                       /// initial position y of galaxy 1 (in grid length fraction)
double z_1 = 1./2;                       /// initial position z of galaxy 1 (in grid length fraction)

double x_2 = 1.25/2;                      /// initial position x of galaxy 2 / SMBH (in grid length fraction) if applicable
double y_2 = 0.1/2;                      /// initial position y of galaxy 2 / SMBH (in grid length fraction) if applicable
double z_2 = 1./2;                       /// initial position z of galaxy 2 / SMBH (in grid length fraction) if applicable

double vx_1 = 0 ;                        /// initial global velocity x of galaxy 1 (in grid length fraction)
double vy_1 = 0.;                        /// initial global velocity y of galaxy 1 (in grid length fraction)
double vz_1 = 0.;                        /// initial global velocity z of galaxy 1 (in grid length fraction)

double vx_2 = 0.;                 /// initial global velocity x of galaxy 2 / SMBH(in grid length fraction) if applicable
double vy_2 = 5;                        /// initial global velocity y of galaxy 2 / SMBH(in grid length fraction) if applicable
double vz_2 = 0.;                        /// initial global velocity z of galaxy 2 / SMBH(in grid length fraction) if applicable
