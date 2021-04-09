/**
* @file utils.c
* @author Cereale
* @brief File containing utility functions
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"

double dist(double* r_1, double* r_2){
  /**
* @brief Cartesian distance between two vectors r_1 and r_2
* @param [in] r_1 3-components vector
* @param [in] r_2 3-components vector
* @return d Distance between r_1 and r_2
*/

  double d = sqrt((r_1[0]-r_2[0])*(r_1[0]-r_2[0]) + (r_1[1]-r_2[1])*(r_1[1]-r_2[1]) + (r_1[2]-r_2[2])*(r_1[2]-r_2[2]));
  return d;
}

double random_double(double min, double max){
/**
* @brief returns an uniform random value in interval [min, max]
* @param [in] min Lower boundary
* @param [in] max Higher boundary
* @return random uniform value
*/
  double range = max-min;
  double div = RAND_MAX/range;
  return (min + (rand()/div));
}

double norm(double* vector){
  /**
* @brief Computes the norm of the 3-components vector "vector"
* @param [in] vector double* array
* @return norm of vector
*/

  double v0 = vector[0], v1 = vector[1], v2 = vector[2];
  return sqrt(v0*v0 + v1*v1 + v2*v2);
}

double* sum_array(double* r_1, double* r_2){
/**
* @brief returns the sum of two 3 components vectors r_1 and r_2
* @param [in] r_1 3-components vector
* @param [in] r_2 3-components vector
* @return out double* array (r_1 + r_2)
*/
  double* out = (double*)malloc(3*sizeof(double));
  out[0] = r_1[0] + r_2[0];
  out[1] = r_1[1] + r_2[1];
  out[2] = r_1[2] + r_2[2];
  free(r_1);
  free(r_2);
  return(out);
}


double* sub_array(double* r_1, double* r_2){
/**
* @brief returns the sum difference two 3 components vectors r_1 and r_2
* @param [in] r_1 3-components vector
* @param [in] r_2 3-components vector
* @return out double* array (r_1 - r_2)
*/
  double* out = (double*)malloc(3*sizeof(double));
  out[0] = r_1[0] - r_2[0];
  out[1] = r_1[1] - r_2[1];
  out[2] = r_1[2] - r_2[2];
  free(r_1);
  free(r_2);
  return(out);
}


double* mult_array(double* r_1, double lambda){
  /**
  * @brief returns the product of a 3 components vectors r_1 with the scalar lambda
  * @param [in] r_1 3-components vector
  * @param [in] lambda double scalar
  * @return out double* array lambda * r_1
  */
  double* out = (double*)malloc(3*sizeof(double));
  out[0] = r_1[0] * lambda;
  out[1] = r_1[1] * lambda;
  out[2] = r_1[2] * lambda;
  return(out);
}

int modulo(int i, int n){
  /**
* @brief returns the value of integer i modulo n
* @param [in] i integer value
* @param [in] n integer value
*/

  int j = i;
  while(j<0){
    j += n;
  }
  while(j>=n){
    j -= n;
  }
  return j;
}

double modulo_double(double i, double n){
  /**
* @brief returns the value of real i modulo n
* @param [in] i double value
* @param [in] n double value
*/
  double j = i;
  while(j<0){
    j += n;
  }
  while(j>=n){
    j -= n;
  }
  return j;
}


void write_parameters(){
  /**
  * @brief Writes all simulation parameters in the file ../data/parameters.txt
*/

  FILE* fichier = fopen("../data/parameters.txt", "w");

  fprintf(fichier, "# Grid :\n");
  fprintf(fichier, "L_x %g\nL_y %g\nL_z %g\nn_x %g\nn_y %g\nn_z %g\ndt %g\nn_t %d\n", L_x, L_y, L_z, n_x, n_y, n_z, dt, n_t);

  fprintf(fichier, "# Methods :\n");
  fprintf(fichier, "acceleration_method %d\ndifferential_method %d\n", acceleration_method, differential_method);
  fprintf(fichier, "second_object %d\n", second_object);
  fprintf(fichier, "m_SMBH %d\n", m_SMBH);

  fprintf(fichier, "# Physics :\n");
  fprintf(fichier, "n_stars %d\n", n_stars);
  fprintf(fichier, "m %g\n", m);
  fprintf(fichier, "fraction_stars %g\n", fraction_stars);
  fprintf(fichier, "a1 %g\n", a1 * L_x);
  fprintf(fichier, "a2 %g\n", a2 * L_x);
  fprintf(fichier, "x_1 %g\ny_1 %g\nz_1 %g\n", x_1 * L_x, y_1 * L_y, z_1 * L_z);
  fprintf(fichier, "x_2 %g\ny_2 %g\nz_2 %g\n", x_2 * L_x, y_2 * L_y, z_2 * L_z);
  fprintf(fichier, "vx_1 %g\nvy_1 %g\nvz_1 %g\n", vx_1, vy_1, vz_1);
  fprintf(fichier, "vx_2 %g\nvy_2 %g\nvz_2 %g\n", vx_2, vy_2, vz_2);



  fclose(fichier);
}


void write_data(FILE* fichier, double** r_stars, double** v_stars, double** a_stars, double* m_stars){
  /**
  * @brief Writes all dynamical data to the file "fichier" to folder ../data/
  * @param [in] fichier FILE* variable through which the file
  * @param [in] r_stars double** array containing all star positions
  * @param [in] v_stars double** array containing all star velocities
  * @param [in] a_stars double* array containing all star accelerations
  * @param [in] m_stars double* array containing all star masses
  */
  int i;
  for (i=0; i<n_stars; i++){
    fprintf(fichier, "%15.5e %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e\n", r_stars[i][0], r_stars[i][1], r_stars[i][2], v_stars[i][0], v_stars[i][1], v_stars[i][2], a_stars[i][0], a_stars[i][1], a_stars[i][2], m_stars[i]);
  }
}

int* get_coordinates(double* r_star){
  /**
* @brief Returns the i, j, k coordinates on the grid of the star at position "r_star"
* @param [in] r_star cartesian coordinates coordinates of a given star
* @return int* array containing the i, j, k coordinates of the star
*/

  int i, j, k;
  int* coordinates = (int*)malloc(3*sizeof(int));
  i = (int)(r_star[0]*n_x/L_x);
  j = (int)(r_star[1]*n_y/L_y);
  k = (int)(r_star[2]*n_z/L_z);
  coordinates[0] = modulo(i, n_x);
  coordinates[1] = modulo(j, n_y);
  coordinates[2] = modulo(k, n_z);
  return coordinates;
}

void copy_space(double*** old_space, double*** new_space){
    /**
  * @brief Assigns to the 3D array new_space values from old_space
  * @param [in] old_space double*** array of size n_x*n_y*n_z
  * @param [in] new_space double*** array of size n_x*n_y*n_z
*/

  int i, j, k;
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      for (k=0; k<n_z; k++){
        new_space[i][j][k] = old_space[i][j][k];
      }
    }
  }
}

double*** initiate_space(){
  /**
* @brief Creates and returns a double*** array of size n_x*n_y*n_z
*/
  int i, j, k;
  double*** space = (double***)malloc(n_x*sizeof(double**));
  for (i=0; i<n_x; i++){
    space[i] = (double**)malloc(n_y*sizeof(double*));
  }
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      space[i][j] = (double*)malloc(n_z*sizeof(double));
    }
  }
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      for (k=0; k<n_z; k++){
        space[i][j][k] = 0.;
      }
    }
  }
  return space;
}

double**** initiate_vector_field(){
  /**
* @brief Creates and returns a double**** array of size n_x*n_y*n_z * 3
*/
  int i, j, k;
  double**** space = (double****)malloc(n_x*sizeof(double***));
  for (i=0; i<n_x; i++){
    space[i] = (double***)malloc(n_y*sizeof(double**));
  }
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      space[i][j] = (double**)malloc(n_z*sizeof(double*));
    }
  }
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      for (k=0; k<n_z; k++){
        space[i][j][k] = (double*)malloc(3*sizeof(double));
      }
    }
  }

  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      for (k=0; k<n_z; k++){
        space[i][j][k][0] = 0.;
        space[i][j][k][1] = 0.;
        space[i][j][k][2] = 0.;
      }
    }
  }
  return space;
}

void free_vector_field(double**** vector_field){
  /**
* @brief Frees all dynamically alocated arrays contained in vector_field and then vector_field itself
* @param [in] vector_field double**** array of size n_x*n_y*n_z * 3
*/
  int i, j, k;
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      for (k=0; k<n_z; k++){
        free(vector_field[i][j][k]);
      }
    }
  }
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      free(vector_field[i][j]);
    }
  }
  for (i=0; i<n_x; i++){
    free(vector_field[i]);
  }
  free(vector_field);
}

void free_space(double*** space){
  /**
  * @brief Frees all dynamically alocated arrays contained in space and then space itself
  * @param [in] space double*** array of size n_x*n_y*n_z
  */
  int i, j;
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      free(space[i][j]);
    }
  }
  for (i=0; i<n_x; i++){
    free(space[i]);
  }
  free(space);
}

void free_vector(double** space){
  /**
* @brief Frees dynamically alocated arrays contained in space and then space itself
* @param [in] space double** array of size n_stars * 3
*/
  int i, j;
  for (i=0; i<n_stars; i++){
    free(space[i]);
  }
  free(space);
}


double max(double a, double b){
  /**
* @brief Returns the maximum between a and b
* @param [in] a real number
* @param [in] b real number
* @return max(a,b)
*/
  if (a>b){
    return a;
  }
  else{
    return b;
  }
}

double min(double a, double b){
  /**
  * @brief Returns the minimum between a and b
  * @param [in] a real number
  * @param [in] b real number
  * @return min(a,b)
  */
  if (a<b){
    return a;
  }
  else{
    return b;
  }
}

double cube_intersect(double* rmin_1, double* rmax_1,double* rmin_2, double* rmax_2){
  /**
* @brief Returns the proportion of cube volume between two cubes of side d_x = L_x/n_x with caracterised each with two opposite cube summits
* @param [in] rmin_1 vector characterising the summit closest to coordinate (0, 0, 0) of cube 1
* @param [in] rmax_1 vector characterising the summit closest to coordinate (L_x, L_y, L_z) of cube 1
* @param [in] rmin_2 vector characterising the summit closest to coordinate (0, 0, 0) of cube 2
* @param [in] rmax_2 vector characterising the summit closest to coordinate (L_x, L_y, L_z) of cube 2
* @return Real value of the area intersecting both cubes
*/

  double O_x, O_y, O_z;
  double d_x = L_x/n_x, d_y = L_y/n_y, d_z = L_z/n_z;
  O_x = max(min(rmax_1[0], rmax_2[0]) - max(rmin_1[0], rmin_2[0]), 0);
  O_y = max(min(rmax_1[1], rmax_2[1]) - max(rmin_1[1], rmin_2[1]), 0);
  O_z = max(min(rmax_1[2], rmax_2[2]) - max(rmin_1[2], rmin_2[2]), 0);
  //printf("--> %15.5g %15.5g %15.5g %15.5g\n", O_x, O_y, O_z, O_x * O_y * O_z / (d_x * d_y * d_z));
  return (O_x*O_y*O_z) / (d_x * d_y * d_z);
}

double**** gradient(double*** field){
/**
* @brief Computes and returns the gradient of the scalar field "field"
* @param [in] field double*** array of size n_x*n_y*n_z
* @return returns double**** array of size n_x*n_y*n_z * 3
*/

  int i, j, k;
  double**** gradient = initiate_vector_field();
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      for (k=0; k<n_z; k++){
        gradient[i][j][k][0] = ((field[modulo(i+1, n_x)][j][k] - field[modulo(i-1, n_x)][j][k]) / (2 * L_x / n_x));
        gradient[i][j][k][1] = (field[i][modulo(j+1, n_y)][k] - field[i][modulo(j-1, n_y)][k]) / (2 * L_y / n_y);
        gradient[i][j][k][2] = (field[i][j][modulo(k+1, n_z)] - field[i][j][modulo(k-1, n_z)]) / (2 * L_z / n_z);
      }
    }
  }
  return gradient;
}

double*** laplacian(double*** field){
  /**
  * @brief Computes and returns the laplacian of the scalar field "field"
  * @param [in] field double*** array of size n_x*n_y*n_z
  * @return returns double*** array of size n_x*n_y*n_z
  */

  int i, j, k;
  double*** laplacian = initiate_space();

  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){
      for (k=0; k<n_z; k++){
        laplacian[i][j][k] += (field[modulo(i+1, n_x)][j][k] - 2 * field[modulo(i, n_x)][j][k] + field[modulo(i-1, n_x)][j][k]) / ((L_x*L_x) / (n_x*n_x));
        laplacian[i][j][k] += (field[i][modulo(j+1, n_y)][k] - 2 * field[i][modulo(j, n_y)][k] + field[i][modulo(j-1, n_y)][k]) / ((L_y*L_y) / (n_y*n_y));
        laplacian[i][j][k] += (field[i][j][modulo(k+1, n_z)] - 2 * field[i][j][modulo(k, n_z)] + field[i][j][modulo(k-1, n_z)]) / ((L_z*L_z) / (n_z*n_z));
      }
    }
  }
  return laplacian;
}

void print_field(double*** field, FILE* field_file){
  /**
* @brief Saves to the file "field_file" all values of the field "field" at z=n_z/2
* @param [in] field double*** array of size n_x*n_y*n_z
*/

  int i, j, k;
  double tot = 0, slice = 0;
  double column;
  int half = 1 + n_x/2;
  for (i=0; i<n_x; i++){
    for (j=0; j<n_y; j++){

      column = 0.;
      for (k=0; k<n_z; k++){
        column += field[i][j][k];
      }
      fprintf(field_file, "%15.5e", field[i][j][half]);
      tot += column;
      slice += field[i][j][half];
    }
    fprintf(field_file, "\n");
  }
  printf("total field at time %g is %g. On slice %g\n", t, tot, slice);
}
