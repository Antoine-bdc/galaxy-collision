/**
* @file main.c
* @author Cereale
* @date 5 Feb 2020
* @brief File main.c of the galaxy collision project. Contains initiation of most of the useful variable and integrates the system over n_t time steps
* Initial conditions and methods used are defined by the values contained in the constant.c file
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "physics.h"
#include "utils.h"
#include <time.h>

int main(){

  /**
* @brief Main function executing all the routines necessary for PM integration
* @details The main function will perform the following tasks :
* - initiation of all variables including : files, arrays containing system quantities,
* - initiation of the initial conditions : Plummer galaxy and potentially an impactor (SMBH, another Plummer)
* - initiation of the grid potential, potential laplacian and gradient
* - updating of grid values using PM method and integrating the system using the corresponding differential method
* - saves dynamical data and grid data
*/

  // Saving seed for reproducibility
  srand(42);

  // Numerical variables
  int i;

  // Initiating writing variables
  FILE* kinematics_file;
  FILE* density_file;
  FILE* potential_file;
  FILE* laplacian_file;
  FILE* energy_file;
  FILE* dt_file;
  char* directory_kinematics = "../data/kinematics/";
  char* file_name = "time_";
  char* directory_potential = "../data/potential/";
  char* directory_density = "../data/density/";
  char* file_name_potential = "potential_";
  char* file_name_density = "density_";
  char* file_type = ".txt";
  char full_name[100] = { };
  char full_name_potential[100] = { };
  char full_name_density[100] = { };


  // Create positions velocities accelerations and mass of stars
  double*** density = initiate_space();
  double*** potential = initiate_potential(density);
  double*** laplacian_field = initiate_space();
  double**** potential_gradient = initiate_vector_field();
  double** r_stars = (double**)malloc(n_stars*sizeof(double*));
  double** v_stars = (double**)malloc(n_stars*sizeof(double*));
  double** a_stars = (double**)malloc(n_stars*sizeof(double*));
  double* m_stars = (double*)malloc(n_stars*sizeof(double));
  double E_k, E_p;

  for (i=0; i<n_stars; i++){
    r_stars[i] = (double*)malloc(3*sizeof(double));
    v_stars[i] = (double*)malloc(3*sizeof(double));
    a_stars[i] = (double*)malloc(3*sizeof(double));
  }

  // Galaxies initiation
  double* center = (double*)malloc(3*sizeof(double));
  double* v_global = (double*)malloc(3*sizeof(double));
  int n1 = fraction_stars*n_stars;
  int n2 = n_stars - n1;

  center[0] = x_1 * L_x;
  center[1] = y_1 * L_y;
  center[2] = z_1 * L_z;
  v_global[0] = vx_1;
  v_global[1] = vy_1;
  v_global[2] = vz_1;

  initiate_a_plummer(center, 0, n1, v_global, a1 * L_x, r_stars, v_stars, a_stars, m_stars);
  if (second_object==1){
    center[0] = x_2 * L_x;
    center[1] = y_2 * L_y;
    center[2] = z_2 * L_z;
    v_global[0] = vx_2 * L_x/(n_t*dt);
    v_global[1] = vy_2 * L_y/(n_t*dt);
    v_global[2] = vz_2 * L_z/(n_t*dt);
    initiate_a_plummer(center, n1, n2, v_global, a2 * L_x, r_stars, v_stars, a_stars, m_stars);
  }
  if (second_object==2){
    center[0] = x_2 * L_x;
    center[1] = y_2 * L_y;
    center[2] = z_2 * L_z;
    v_global[0] = vx_2 * L_x/(n_t*dt);
    v_global[1] = vy_2 * L_y/(n_t*dt);
    v_global[2] = vz_2 * L_z/(n_t*dt);
    initiate_SMBH(center, v_global, r_stars, v_stars, a_stars, m_stars);
  }

  free(center);
  free(v_global);

  laplacian_file = fopen("../data/laplacian.txt", "w");
  energy_file = fopen("../data/energy.txt", "w");
  dt_file = fopen("../data/dt.txt", "w");

  // Initiation of the potential array with corresponding stars
  if (acceleration_method == 2){
    update_density_1(density, r_stars, m_stars);
    if (potential_method==0){
      //jacobi_potential(potential, density);
    }
    if (potential_method==1){
      gauss_seidel_potential(potential, density);
    }
    if (potential_method == 2){
      SOR_potential(potential, density);
    }
    laplacian_field = laplacian(potential);
    potential_gradient = gradient(potential);
  }


  if (dynamical_dt==1){
    double v_max = sqrt(2 * G * m / (a1));
    dt = ((L_x/n_x) / v_max)*C_vel;
    printf("dt initial :%g\n", dt);
  }
  // Save simulation parameters
  write_parameters();

  // Saving data at time 0
  sprintf(full_name, "%s%s%06d%s", directory_kinematics, file_name, 0, file_type);
  sprintf(full_name_potential, "%s%s%06d%s", directory_potential, file_name_potential, 0, file_type);
  sprintf(full_name_density, "%s%s%06d%s", directory_density, file_name_density, 0, file_type);

  kinematics_file = fopen(full_name, "w");
  density_file = fopen(full_name_density, "w");
  potential_file = fopen(full_name_potential, "w");

  write_data(kinematics_file, r_stars, v_stars, a_stars, m_stars);
  print_field(density, density_file);
  print_field(potential, potential_file);

  fclose(kinematics_file);
  fclose(density_file);
  fclose(potential_file);

  // Integrating using differential_method
  for (i=0; i<n_t; i++){
    t += dt;
    printf("Time : %g  (%d / %d)      ---      ", t, i, n_t);

    // Integrating dynamical variables
    if (differential_method==1){
      LEAP_FROG(dt, r_stars, v_stars, a_stars, m_stars, potential_gradient);
    }
    if (differential_method==0){
      EULER(dt, r_stars, v_stars, a_stars, m_stars, potential_gradient);
    }
    if (differential_method==2){
      KICK_DRIFT_KICK(dt_file, r_stars, v_stars, a_stars, m_stars, potential_gradient);
    }
    reposition(r_stars);

    // Computing physical quantities
    if (acceleration_method == 2){

      update_density_1(density, r_stars, m_stars);
      if (potential_method==0){
        jacobi_potential(potential, density);
      }
      if (potential_method==1){
        gauss_seidel_potential(potential, density);
      }
      if (potential_method == 2){
        SOR_potential(potential, density);
      }
      free_vector_field(potential_gradient);
      potential_gradient = gradient(potential);
    }

    if (n_stars<100001){ // If too many stars then PP calculation time of potential energy explodes
      E_k = kinetic_energy(r_stars, v_stars, m_stars);
      E_p = potential_energy(r_stars, v_stars, m_stars);
      fprintf(energy_file, "%15.5e %15.5e %15.5e %15.5e\n", t, E_k, E_p, E_p + E_k);
    }

    // Saving data (phase space, potential, density)
    sprintf(full_name, "%s%s%06d%s", directory_kinematics, file_name, i, file_type);
    sprintf(full_name_potential, "%s%s%06d%s", directory_potential, file_name_potential, i, file_type);
    sprintf(full_name_density, "%s%s%06d%s", directory_density, file_name_density, i, file_type);

    kinematics_file = fopen(full_name, "w");
    density_file = fopen(full_name_density, "w");
    potential_file = fopen(full_name_potential, "w");

    print_field(density, density_file);
    print_field(potential, potential_file);
    write_data(kinematics_file, r_stars, v_stars, a_stars, m_stars); // rx, ry, rz, vx, vy, vz, ax, ay, az, m

    fclose(kinematics_file);
    fclose(density_file);
    fclose(potential_file);
  }

  print_field(laplacian_field, laplacian_file);

  // Free memory :
  free_vector_field(potential_gradient);
  free_space(density);
  free_space(potential);
  free_space(laplacian_field);
  free_vector(r_stars);
  free_vector(v_stars);
  free_vector(a_stars);

  free(m_stars);
  free(v_global);
  free(center);
  //fclose(energy_file);
  fclose(density_file);
  fclose(potential_file);
  fclose(laplacian_file);
  fclose(dt_file);

  return 0;
}
