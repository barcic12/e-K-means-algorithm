#ifndef FINALPROJECT_SPKMEANS_H
#define FINALPROJECT_SPKMEANS_H
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

struct p_matrix{
    int i;
    int j;
    double c;
    double s;
};

/* functions declaration: */
static int apply_requested_algo(int goal, double **dots,double **matrix, double* array, int n, int d, int k, int print);
static int spk_goal(double **dots, double **matrix, double *dots_location, int n, int d, int k);
static void wam_goal(double **dots, double **matrix, int n, int d);
static void ddg_goal (double **dots, double **matrix, double* array, int n, int d);
static void lnorm_goal (double **dots, double **matrix, double *array, int n, int d);
static void convert_diagonal_to_matrix(double **matrix, double *array, int n);
static void jacobi_goal (double **dots, double **matrix, double* array, int n);
static void calc_matrix_multi_p(double **matrix, struct p_matrix p, int n);
static double jacobi_step(int cnt, struct p_matrix p_matrices[], double** matrix, int n);
static void calc_p_transpose_multi_matrix_multi_p(double **matrix,struct p_matrix p, int n);
static void calc_eigen_vector_matrix(int cnt, struct p_matrix p_matrices[], double** matrix, int n);
static int calc_k(double *eigen_values, int n);
static void kmeans_calc(double **dots,double **centroids, double *dots_location, int n, int d, int k);
static void sort_eigenvalues_and_eigenvectors(double *array, double **matrix, int n);
static double calc_distance(double *dot, double *centroid, int d);
static int find_nearest_centroid(double *dot, double **centroids, int k, int d);
static void update_centroids(double **dots,double *dots_location, double **new_centroids, int n, int d, int k);
static int check_equals_2d_list(double **list1, double **list2, int row_num, int col_num);
static void print_2d_array(double **array, int row_num, int col_num);
static double** get_new_matrix(int n, int d);

#endif


