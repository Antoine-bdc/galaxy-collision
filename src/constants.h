// Numerical constants:
double L_x;
double L_y;
double L_z;
double n_x;
double n_y;
double n_z;
double C_vel;
double C_dyn;
double t;

// Physical constants:
double G;
double M_tot;

// Methods :
int acceleration_method;         // (0 = Central), (1 = PP), (2 = PM),
int differential_method;         // (0 = EULER), (1 = LEAPFROG), (2 = KICK_DRIFT_KICK)
int potential_method;            // (0 = Jacobi), (1 = Gauss_Seidel), (2 = SOR)
int dynamical_dt;                // (0 = False), (1 = True)

// Simulation parameters:
double dt;
int n_t;
double jacobi_precision;
double gauss_seidel_precision;
int second_object;               // (0 = None), (1 = galaxy 2), (2 = SMBH)

// initial conditions
double m;
double m_SMBH;
int n_stars;
double fraction_stars;               // fraction of stars in the galaxy 1
double a1;                           // galaxy 1 typical size (in grid length fraction)
double a2;                           // galaxy 1 typical size (in grid length fraction)
// initial position of galaxy 1 (in grid length fraction)
double x_1;
double y_1;
double z_1;
// initial position of galaxy 2 (in grid length fraction)
double x_2;
double y_2;
double z_2;
// initial velocity of galaxy 1 (in grid length fraction)
double vx_1;
double vy_1;
double vz_1;
// initial velocity of galaxy 2 (in grid length fraction)
double vx_2;
double vy_2;
double vz_2;
