/**
* @file physics.c
* @author Antoine
* @brief This file contains all the functions related to the physical computations of the system.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "constants.h"


void reposition(double** r_stars){
  /**
* @brief moves stars that escaped from the grid to a location inside the grid using periodical conditions
* @param [in] r_stars Position of all stars
*/

  int l;
  for (l=0; l<n_stars; l++){
    r_stars[l][0] = modulo_double(r_stars[l][0], L_x);
    r_stars[l][1] = modulo_double(r_stars[l][1], L_y);
    r_stars[l][2] = modulo_double(r_stars[l][2], L_z);
  }
}

double* central_acceleration(int i_star, double** r_stars, double* m_stars){
  /**
* @brief Computes the acceleration on star of index i_star from a central potential of mass M_star (Global parameter)
* @param [in] i_star Star index
* @param [in] r_stars double** array containing all star positions
* @param [in] m_stars double* array containing all star masses
* @return Returns the acceleration for the star of index i_star as a double* array
*/

  int i;
  double d;
  double* acc = (double*)malloc(3*sizeof(double));
  double* center = (double*)malloc(3*sizeof(double));
  center[0] = L_x/2;
  center[1] = L_y/2;
  center[2] = L_z/2;
  acc[0] = 0.;
  acc[1] = 0.;
  acc[2] = 0.;

  d = dist(r_stars[i_star], center);
  // Central acceleration from Newtonian physics
  acc[0] += (r_stars[i_star][0] - center[0]) * (-G * m_stars[0]) / (d*d*d);
  acc[1] += (r_stars[i_star][1] - center[1]) * (-G * m_stars[0]) / (d*d*d);
  acc[2] += (r_stars[i_star][2] - center[2]) * (-G * m_stars[0]) / (d*d*d);
  free(center);
  return acc;
}


// Computes acceleration on star i in the field of all other masses
double* PP_acceleration(int i_star, double** r_stars, double* m_stars){
  /**
* @brief Computes the acceleration on star of index i_star from the interaction from all other stars using PM method
* @param [in] i_star Star index
* @param [in] r_stars double** array containing all star positions
* @param [in] m_stars double* array containing all star masses
* @return Returns the acceleration for the star of index i_star as a double* array
*/
  int i;
  double* acc = (double*)malloc(3*sizeof(double));
  double d;

  acc[0] = 0.;
  acc[1] = 0.;
  acc[2] = 0.;

  // PP potential calculation
  for (i=0; i<n_stars; i++){
      if (i!=i_star){
      	d = dist(r_stars[i_star], r_stars[i]);
      	acc[0] += (r_stars[i_star][0] - r_stars[i][0]) * (-G * m_stars[i]) / (d*d*d);
        acc[1] += (r_stars[i_star][1] - r_stars[i][1]) * (-G * m_stars[i]) / (d*d*d);
        acc[2] += (r_stars[i_star][2] - r_stars[i][2]) * (-G * m_stars[i]) / (d*d*d);
      }
    }
  return acc;
}


double* PM_acceleration(double* r_star, double**** potential_gradient){
  /**
* @brief Computes the acceleration on star located at position r_star using PM méthod
* @param [in] r_star Star position
* @param [in] potential_gradient Gradient of the gravitationnal potential
* @return Returns the acceleration for the star in position r_star as a double* array
*/

  // Star iterator
  int l, i, j, k;
  int i_star, j_star, k_star;
  // Grid iterator
  double x, y, z;
  double intersect;
  double d_x = L_x/n_x, d_y = L_y/n_y, d_z = L_z/n_z;
  double* cell_r_min = (double*)malloc(3*sizeof(double));
  double* cell_r_max = (double*)malloc(3*sizeof(double));
  double* star_r_min = (double*)malloc(3*sizeof(double));
  double* star_r_max = (double*)malloc(3*sizeof(double));
  double* acc = (double*)malloc(3*sizeof(double));
  acc[0] = 0.;
  acc[1] = 0.;
  acc[2] = 0.;

  int* coord = get_coordinates(r_star);
  i_star = coord[0];
  j_star = coord[1];
  k_star = coord[2];
  free(coord);
  // Set cubic star sides coordinates
  star_r_min[0] = modulo_double(r_star[0] - d_x/2., L_x);
  star_r_min[1] = modulo_double(r_star[1] - d_y/2., L_y);
  star_r_min[2] = modulo_double(r_star[2] - d_z/2., L_z);
  star_r_max[0] = modulo_double(r_star[0] + d_x/2., L_x);
  star_r_max[1] = modulo_double(r_star[1] + d_y/2., L_y);
  star_r_max[2] = modulo_double(r_star[2] + d_z/2., L_z);

  // Looking at all 27 adjacent cells
  for (i=i_star; i<=i_star+2; i++){
    for (j=j_star; j<=j_star+2; j++){
      for (k=k_star; k<=k_star+2; k++){
        // Set cell sides coordinates
        cell_r_min[0] = modulo_double(((double)(i) - 1) * d_x, L_x);
        cell_r_min[1] = modulo_double(((double)(j) - 1) * d_y, L_y);
        cell_r_min[2] = modulo_double(((double)(k) - 1) * d_z, L_z);
        cell_r_max[0] = modulo_double(((double)(i)) * d_x, L_x);
        cell_r_max[1] = modulo_double(((double)(j)) * d_y, L_y);
        cell_r_max[2] = modulo_double(((double)(k)) * d_z, L_z);
        intersect = cube_intersect(cell_r_min, cell_r_max, star_r_min, star_r_max);
        acc[0] += - 1 * intersect * potential_gradient[modulo(i, n_x)][modulo(j, n_y)][modulo(k, n_z)][0];
        acc[1] += - 1 * intersect * potential_gradient[modulo(i, n_x)][modulo(j, n_y)][modulo(k, n_z)][1];
        acc[2] += - 1 * intersect * potential_gradient[modulo(i, n_x)][modulo(j, n_y)][modulo(k, n_z)][2];
      }
    }
  }

  free(cell_r_max);
  free(cell_r_min);
  free(star_r_min);
  free(star_r_max);
  return acc;
}


void LEAP_FROG(double local_dt, double** r_stars, double** v_stars, double** a_stars, double* m_stars, double**** potential_gradient){
  /**
  * @brief Updates positions and velocities of all stars after the time step dt using the LEAPFROG differential method
  * @param [in] local_dt time step for the LEAPFROG algorithm
  * @param [in] r_stars double** array containing all star positions
  * @param [in] v_stars double** array containing all star velocities
  * @param [in] a_stars double* array containing all star accelerations
  * @param [in] m_stars double* array containing all star masses
  * @param [in] potential_gradient Gradient of the gravitationnal potential
*/

  int i;
  for (i=0; i<n_stars; i++){
    r_stars[i] = sum_array(r_stars[i], sum_array(mult_array(v_stars[i], local_dt), mult_array(a_stars[i], local_dt*local_dt/2)));
    v_stars[i] = sum_array(v_stars[i], mult_array(a_stars[i], local_dt/2));
  }
  for (i=0; i<n_stars; i++){
    if (acceleration_method == 0){
      a_stars[i] = central_acceleration(i, r_stars, m_stars);
    }
    if (acceleration_method == 1){
      a_stars[i] = PP_acceleration(i, r_stars, m_stars);
    }
    if (acceleration_method == 2){
      a_stars[i] = PM_acceleration(r_stars[i], potential_gradient);
    }
  }
  for (i=0; i<n_stars; i++){
    v_stars[i] = sum_array(v_stars[i], mult_array(a_stars[i], local_dt/2));
  }
}


void KICK_DRIFT_KICK(FILE* dt_file, double** r_stars, double** v_stars, double** a_stars, double* m_stars, double**** potential_gradient){
  /**
  * @brief Updates positions and velocities of all stars after the time step dt using the KICK DRIFT KICK differential method
  * @param [in] dt_file File variable to write the evolution of dt over time
  * @param [in] r_stars double** array containing all star positions
  * @param [in] v_stars double** array containing all star velocities
  * @param [in] a_stars double* array containing all star accelerations
  * @param [in] m_stars double* array containing all star masses
  * @param [in] potential_gradient Gradient of the gravitationnal potential
  */
  int i;
  double v, v_max = 0.;
  // Kick
  for (i=0; i<n_stars; i++){
    v_stars[i] = sum_array(v_stars[i], mult_array(a_stars[i], dt/2));
  }
  // Drift
  for (i=0; i<n_stars; i++){
    r_stars[i] = sum_array(r_stars[i], mult_array(v_stars[i], dt));
  }

  for (i=0; i<n_stars; i++){
    if (acceleration_method == 0){
      a_stars[i] = central_acceleration(i, r_stars, m_stars);
    }
    if (acceleration_method == 1){
      a_stars[i] = PP_acceleration(i, r_stars, m_stars);
    }
    if (acceleration_method == 2){
      a_stars[i] = PM_acceleration(r_stars[i], potential_gradient);
    }
  }
  // Kick
  for (i=0; i<n_stars; i++){
    v_stars[i] = sum_array(v_stars[i], mult_array(a_stars[i], dt/2));
  }
  // Updating dt
  for (i=1; i<n_stars; i++){
    v = v_stars[i][0]*v_stars[i][0] + v_stars[i][1]*v_stars[i][1] + v_stars[i][2]*v_stars[i][2];
    //printf("HERE :   %g %g %g %g \n",dt, v_stars[i][0], v, v_max);
    if (v_max < v){
      v_max = v;
    }
  }
  printf("vmax %g => dt = %g\n", v_max, C_vel * L_x / (n_x * v_max));
  if (dynamical_dt==1){
    dt = C_vel * L_x / (n_x * v_max);
    printf("dt = %g\n", dt);
  }
  fprintf(dt_file, "%g %g \n",t , dt);

}

void EULER(double local_dt, double** r_stars, double** v_stars, double** a_stars, double* m_stars, double**** potential_gradient){
    /**
    * @brief Updates positions and velocities of all stars after the time step dt using the LEAPFROG differential method
    * @param [in] local_dt time step for the LEAPFROG algorithm
    * @param [in] r_stars double** array containing all star positions
    * @param [in] v_stars double** array containing all star velocities
    * @param [in] a_stars double* array containing all star accelerations
    * @param [in] m_stars double* array containing all star masses
    * @param [in] potential_gradient Gradient of the gravitationnal potential
  */

  int i;
  for (i=0; i<n_stars; i++){
    v_stars[i] = sum_array(v_stars[i], mult_array(a_stars[i], local_dt));
    r_stars[i] = sum_array(r_stars[i], mult_array(v_stars[i], local_dt));
  }
  for (i=0; i<n_stars; i++){
    if (acceleration_method == 0){
      a_stars[i] = central_acceleration(i, r_stars, m_stars);
    }
    if (acceleration_method == 1){
      a_stars[i] = PP_acceleration(i, r_stars, m_stars);
    }
    if (acceleration_method == 2){
      a_stars[i] = PM_acceleration(r_stars[i], potential_gradient);
    }
  }
}

double g_dist(double q){
  /**
* @brief Returns the distribution function of velocity
* @param [in] q value of q
* @return g(q)
*/

  return q*q*pow(1 - q*q, 7./2.);
}

void initiate_SMBH(double* center, double* v_SMBH, double** r_stars, double** v_stars, double** a_stars, double* m_stars){
  /**
  * @brief Sets the 1st star to be a supermassive black hole with position center and velocity v_SMBH
  * @param [in] center position of the SMBH
  * @param [in] v_SMBH velocity vector of SMBH
  * @param [in] r_stars double** array containing all star positions
  * @param [in] v_stars double** array containing all star velocities
  * @param [in] a_stars double* array containing all star accelerations
  * @param [in] m_stars double* array containing all star masses
*/

  r_stars[0][0] = center[0];
  r_stars[0][1] = center[1];
  r_stars[0][2] = center[2];
  v_stars[0][0] = v_SMBH[0];
  v_stars[0][1] = v_SMBH[1];
  v_stars[0][2] = v_SMBH[2];
  m_stars[0] = m_SMBH;
  int j;
  for (j=0; j<n_stars; j++){
    if (acceleration_method == 0){
      a_stars[j] = central_acceleration(j, r_stars, m_stars);
    }
    if (acceleration_method == 1){
      a_stars[j] = PP_acceleration(j, r_stars, m_stars);
    }
    if (acceleration_method == 2){
      a_stars[j] = central_acceleration(j, r_stars, m_stars);
    }
  }
}

void initiate_plummer(double** r_stars, double** v_stars, double** a_stars, double* m_stars){
  /**
  * @brief Initiates the Plummer model galaxy distribition
  * @param [in] r_stars double** array containing all star positions
  * @param [in] v_stars double** array containing all star velocities
  * @param [in] a_stars double* array containing all star accelerations
  * @param [in] m_stars double* array containing all star masses
*/
  int i, j, reject;
  double a = L_x/5;
  double m, d, phi, v, rdn, M, ve, theta;
  M = n_stars * m;
  int n_reject = 0;

  // Initiate positions
  for (i=0; i<n_stars; i++){
    m_stars[i] = m;
    rdn = random_double(0., 1.);
    d = a / sqrt(pow(rdn, -2./3.) - 1);
    phi = 2 * M_PI * random_double(0., 1.);
    theta = acos(1 - 2 * random_double(0., 1.));
    r_stars[i][0] = L_x/2 + d * sin(theta) * cos(phi);
    r_stars[i][1] = L_y/2 + d * sin(theta) * sin(phi);
    r_stars[i][2] = L_z/2 + d * cos(theta);
  }

  reposition(r_stars);

  // Initiate speed :
  for (i=0; i<n_stars; i++){
    ve = sqrt(2*G*M) * pow(d*d + a*a, -1./4.);
    reject = 1;
    while(reject==1){
      rdn = random_double(0, 1);
      v = random_double(0, 1);
      if (g_dist(rdn)>v){
        v = rdn * ve; // * ve
        reject = 0;
      }
      else{
        n_reject += 1;
      }
    }
    phi = 2 * M_PI * random_double(0., 1.);
    theta = acos(1 - 2 * random_double(0., 1.));

    v_stars[i][0] = v * sin(theta) * cos(phi);
    v_stars[i][1] = v * sin(theta) * sin(phi);
    v_stars[i][2] = v * cos(theta);

  }

  printf("total rejects : %d out of %d\n", n_reject, n_stars);

  for (j=0; j<n_stars; j++){
    if (acceleration_method == 0){
      a_stars[j] = central_acceleration(j, r_stars, m_stars);
    }
    if (acceleration_method == 1){
      a_stars[j] = PP_acceleration(j, r_stars, m_stars);
    }
    if (acceleration_method == 2){
      a_stars[j] = central_acceleration(j, r_stars, m_stars);
    }
  }
}

void initiate_a_plummer(double* center, int starting_i, int galaxy_stars, double* v_global, double a, double** r_stars, double** v_stars, double** a_stars, double* m_stars){
  /**
  * @brief Initiates the Plummer model galaxy distribition with center "center" and global velocity "v_global" for a total of "galaxy_stars" stars
  * @param [in] center position of the Plummer distribution
  * @param [in] starting_i index of the first star of this galaxy
  * @param [in] galaxy_stars Total number of stars of the Plummer distribution
  * @param [in] v_global global velocity for all stars in the galaxy
  * @param [in] a parameter a of Plummer distribution ~ typical size of the galaxy
  * @param [in] r_stars double** array containing all star positions
  * @param [in] v_stars double** array containing all star velocities
  * @param [in] a_stars double* array containing all star accelerations
  * @param [in] m_stars double* array containing all star masses
*/

  // random distance and phase for each star
  int i, j, reject;
  double d, phi, v, rdn, M, ve, theta;
  M = galaxy_stars * m;
  int n_reject = 0;
  // Initiate positions
  for (i=starting_i; i<starting_i + galaxy_stars; i++){
    rdn = random_double(0., 1.);
    d = a / sqrt(pow(rdn, -2./3.) - 1);
    phi = 2 * M_PI * random_double(0., 1.);
    theta = acos(1 - 2 * random_double(0., 1.));
    m_stars[i] = m;
    r_stars[i][0] = center[0] + d * sin(theta) * cos(phi);
    r_stars[i][1] = center[1] + d * sin(theta) * sin(phi);
    r_stars[i][2] = center[2] + d * cos(theta);

    ve = sqrt(2*G*M) * pow(d*d + a*a, -1./4.);
    reject = 1;
    while(reject==1){
      rdn = random_double(0., 1.);
      v = random_double(0., 1.);
      if (v < g_dist(rdn)){
        //printf("%15.5g %15.5g %15.5g\n", v, rdn, rdn * ve);
        v = rdn * ve;
        reject = 0;
      }
      else{
        n_reject += 1;
      }
    }
    phi = 2 * M_PI * random_double(0., 1.);
    theta = acos(1 - 2 * random_double(0., 1.));
    v_stars[i][0] = v_global[0] + v * sin(theta) * cos(phi);
    v_stars[i][1] = v_global[1] + v * sin(theta) * sin(phi);
    v_stars[i][2] = v_global[2] + v * cos(theta);
  }

  reposition(r_stars);

  printf("total rejects : %d out of %d\n", n_reject, n_stars);

  for (j=starting_i; j<starting_i + galaxy_stars; j++){
    if (acceleration_method == 0){
      a_stars[j] = central_acceleration(j, r_stars, m_stars);
    }
    if (acceleration_method == 1){
      a_stars[j] = PP_acceleration(j, r_stars, m_stars);
    }
    if (acceleration_method == 2){
      a_stars[j] = central_acceleration(j, r_stars, m_stars);
    }
  }
}

double potential_energy(double** r_stars, double** v_stars, double* m_stars){
  /**
* @brief Rerturns the potential energy of the system using newtonian PP formula
* @param [in] r_stars double** array containing all star positions
* @param [in] v_stars double** array containing all star velocities
* @param [in] m_stars double* array containing all star masses
* @return double with the total potential energy
*/

  double E_p = 0.;
  int i, j;
  for (i=1; i<n_stars-1; i++){
    for(j=i+1; j<n_stars; j++){
      E_p += -G * m_stars[i] * m_stars[j] / dist(r_stars[i], r_stars[j]);
    }
  }
  return E_p;
}

double kinetic_energy(double** r_stars, double** v_stars, double* m_stars){
  /**
* @brief Returns the kinetic energy of the system
* @param [in] r_stars double** array containing all star positions
* @param [in] v_stars double** array containing all star velocities
* @param [in] m_stars double* array containing all star masses
* @return double with the total potential energy
*/

  double E_k = 0.;
  double* zero = (double*)malloc(3*sizeof(double));
  zero[0] = 0.;
  zero[1] = 0.;
  zero[2] = 0.;
  int i;
  for (i=1; i<n_stars; i++){
    E_k += 0.5 * m_stars[i] * (dist(v_stars[i], zero) * dist(v_stars[i], zero));
  }
  free(zero);
  return E_k;
}


void update_density_0(double*** density, double** r_stars, double* m_stars){
  /**
* @brief Updates the density grid using star positions and masses.  Order 0 of density computation (Next Grid Point)
* @param [in] density double** array containing all star positions
* @param [in] v_stars double** array containing all star velocities
* @param [in] m_stars double* array containing all star masses
*/
  // star iterator
  int l, i, j, k;
  // grid iterator
  int* coord;
  double x, y, z;

  for (i=0;i<n_x;i++){
      for (j=0;j<n_y;j++){
          for (k=0;k<n_z;k++){
            density[i][j][k] = 0;
          }
      }
  }
  for (l=0; l<n_stars; l++){
    coord = get_coordinates(r_stars[l]);
    i = coord[0];
    j = coord[1];
    k = coord[2];
    density[i][j][k] += m_stars[l] / (L_x/n_x)*(L_y/n_y)*(L_z/n_z);
    free(coord);
  }

  // printf("Total out : %d out of %d coordinates", nb_out, 3*n_stars);
}

void update_density_1(double*** density, double** r_stars, double* m_stars){
  /**
* @brief Updates the density grid using star positions and masses.  Order 1 of density computation (Cloud grid cell)
* @param [in] density double** array containing all star positions
* @param [in] v_stars double** array containing all star velocities
* @param [in] m_stars double* array containing all star masses
*/
  // Star iterator
  int l, i, j, k;
  // Grid iterator
  int i_star, j_star, k_star;
  int* coord;
  double x, y, z;
  double intersect;
  double d_x = L_x/n_x, d_y = L_y/n_y, d_z = L_z/n_z;
  double* cell_r_min = (double*)malloc(3*sizeof(double));
  double* cell_r_max = (double*)malloc(3*sizeof(double));
  double* star_r_min = (double*)malloc(3*sizeof(double));
  double* star_r_max = (double*)malloc(3*sizeof(double));

  for (i=0;i<n_x;i++){
      for (j=0;j<n_y;j++){
          for (k=0;k<n_z;k++){
            density[i][j][k] = 0;
          }
      }
  }
  for (l=0; l<n_stars; l++){
    coord = get_coordinates(r_stars[l]);
    i_star = coord[0];
    j_star = coord[1];
    k_star = coord[2];
    free(coord);
    // Set cubic star sides coordinates
    star_r_min[0] = modulo_double(r_stars[l][0] - d_x/2., L_x);
    star_r_min[1] = modulo_double(r_stars[l][1] - d_y/2., L_y);
    star_r_min[2] = modulo_double(r_stars[l][2] - d_z/2., L_z);
    star_r_max[0] = modulo_double(r_stars[l][0] + d_x/2., L_x);
    star_r_max[1] = modulo_double(r_stars[l][1] + d_y/2., L_y);
    star_r_max[2] = modulo_double(r_stars[l][2] + d_z/2., L_z);

    // Looking at all 27 adjacent cells
    for (i=i_star; i<=i_star+2; i++){
      for (j=j_star; j<=j_star+2; j++){
        for (k=k_star; k<=k_star+2; k++){
          // Set cell sides coordinates
          cell_r_min[0] = modulo_double(((double)(i) - 1) * d_x, L_x);
          cell_r_min[1] = modulo_double(((double)(j) - 1) * d_y, L_y);
          cell_r_min[2] = modulo_double(((double)(k) - 1) * d_z, L_z);
          cell_r_max[0] = modulo_double(((double)(i)) * d_x, L_x);
          cell_r_max[1] = modulo_double(((double)(j)) * d_y, L_y);
          cell_r_max[2] = modulo_double(((double)(k)) * d_z, L_z);
          intersect = cube_intersect(star_r_min, star_r_max, cell_r_min, cell_r_max) * m_stars[l]  / (d_x*d_y*d_z);
          //printf("intersect : %g \nm_star : %g \n(d_x*d_y*d_z) : %g \nproportion : %g \n\n", intersect, m_stars[l], (d_x*d_y*d_z), cube_intersect(star_r_min, star_r_max, cell_r_min, cell_r_max));
          density[modulo(i, n_x)][modulo(j, n_y)][modulo(k, n_z)] += intersect;
        }
      }
    }
  }

  free(cell_r_max);
  free(cell_r_min);
  free(star_r_min);
  free(star_r_max);
}

double*** initiate_potential(double*** density){
  /**
* @brief First iteration of the gravitationnal potential array. Initial condition set by the Poisson equation with Laplacian(phi) = phi/dx²
* @param [in] density Density 3D array
* @return potential Gravitationnal potential 3D  array
*/
  int i, j, k;
  double*** potential = initiate_space();
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      for (k=0; k<n_z; k++){
        potential[i][j][k] = 4 * M_PI * G * density[i][j][k] * (L_x/n_x)*(L_x/n_x);
      }
    }
  }
  return potential;
}

void jacobi_potential(double*** potential, double*** density){
  /**
  * @brief Updates the potential until precision "jacobi_precision" is met for the set density using jacobi relaxation method
  * @param [in] potential Gravitationnal potentail 3D array
  * @param [in] density Density 3D array
  */
  double*** previous_potential = initiate_space();
  double eps0 = 1., eps1 = 100.;
  int n=0, i, j, k;
  while (fabs(eps1/eps0 - 1.) > jacobi_precision){
    copy_space(potential, previous_potential);
    eps0 = eps1;
    eps1 = 0.;
    for (i=0; i<n_x; i++){
      for (j=0; j<n_y; j++){
        for (k=0; k<n_z; k++){
          potential[i][j][k] = (1./6) * (previous_potential[modulo(i+1, n_x)][j][k] + previous_potential[i][modulo(j+1, n_y)][k] + previous_potential[i][j][modulo(k+1, n_z)] + previous_potential[modulo(i-1, n_x)][j][k] + previous_potential[i][modulo(j-1, n_y)][k] + previous_potential[i][j][modulo(k-1, n_z)] - 4*M_PI*G * density[i][j][k] * (L_x/n_x)*(L_x/n_x));
          eps1 += potential[i][j][k]*potential[i][j][k];
        }
      }
    }
    n+=1;
    eps1 = sqrt(eps1);
    // printf("%d \n", n);
  }
  free_space(previous_potential);
}

void gauss_seidel_potential(double*** potential, double*** density){
  /**
  * @brief Updates the potential until precision "gauss_seidel_precision" is met for the set density using Gauss Seidel relaxation method
  * @param [in] potential Gravitationnal potentail 3D array
  * @param [in] density Density 3D array
  */
  double eps0 = 1., eps1 = 100.;
  int n=0, i, j, k;
  while (fabs(eps1/eps0 - 1.) > gauss_seidel_precision){
    eps0 = eps1;
    eps1 = 0.;
    // Black cells
    for (i=0; i<n_x; i++){
      for (j=0; j<n_y; j++){
        for (k=0; k<n_z; k++){
            potential[i][j][k] = (1./6) * (potential[modulo(i+1, n_x)][j][k] + potential[i][modulo(j+1, n_y)][k] + potential[i][j][modulo(k+1, n_z)] + potential[modulo(i-1, n_x)][j][k] + potential[i][modulo(j-1, n_y)][k] + potential[i][j][modulo(k-1, n_z)] - 4*M_PI*G * density[i][j][k] * (L_x/n_x)*(L_x/n_x));
            if (modulo(i+j+k, 2) == 0){
            eps1 += potential[i][j][k] * potential[i][j][k];
          }
        }
      }
    }
    // Red cells
    for (i=0; i<n_x; i++){
      for (j=0; j<n_y; j++){
        for (k=0; k<n_z; k++){
          if (modulo(i+j+k, 2) == 1){
            potential[i][j][k] = (1./6) * (potential[modulo(i+1, n_x)][j][k] + potential[i][modulo(j+1, n_y)][k] + potential[i][j][modulo(k+1, n_z)] + potential[modulo(i-1, n_x)][j][k] + potential[i][modulo(j-1, n_y)][k] + potential[i][j][modulo(k-1, n_z)] - 4*M_PI*G * density[i][j][k] * (L_x/n_x)*(L_x/n_x));
            eps1 += potential[i][j][k] * potential[i][j][k];
          }
        }
      }
    }
    eps1 = sqrt(eps1);

    if (n>0){
      // printf("%d (%g)\n", n, eps1/eps0 - 1);
    }
    n+=1;
  }
}

void SOR_potential(double*** potential, double*** density){
  /**
  * @brief Updates the potential until precision "gauss_seidel_precision" is met for the set density using Successive Over Relaxation method
  * @param [in] potential Gravitationnal potentail 3D array
  * @param [in] density Density 3D array
  */
  int i, j, k;
  double omega = 2 /(1 + 2 * M_PI * pow(n_x * n_y * n_z, 1./3));
  double*** previous_potential = initiate_space();
  copy_space(potential, previous_potential);
  gauss_seidel_potential(potential, density);
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      for (k=0; k<n_z; k++){
        potential[i][j][k] = previous_potential[i][j][k] + omega * (potential[i][j][k] - previous_potential[i][j][k]);
      }
    }
  }
  free_space(previous_potential);
}
