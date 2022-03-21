
#include "spkmeans.h"

#define MAX_ITER 300
/* check e */
#define  MAX_JACOBI_ITER 100
#define MAX_BUFFER 1000
#define INPUT_ERROR                                     \
        {                                               \
        printf("Invalid Input!");                        \
        assert(0);                                      \
        }
#define ERROR                                            \
        {                                                \
        printf("An Error Has Occured");                  \
        assert(0);                                       \
        }

/* ----------- main ----------------
 * implement the requested algorithm while running this c file directly */

int main(int argc, char *argv[]){
    /* variables declaration: */
    int k, n=0, d=0, i, j, goal=0;
    double **dots, **matrix;
    double *array;
    char *line, *tokens, *file_name;
    const char* conversion[] = {"spk","wam","ddg","lnorm","jacobi"};
    FILE *input_file, *input_file2;

    if(argc!=4) INPUT_ERROR
    line = malloc(MAX_BUFFER*sizeof(char));
    if(!line) ERROR

    /* getting the arguments: */
    k = strtol(argv[1], NULL, 10);
    for(i=0; i<5; i++)
        if(strcmp((const char*)argv[2],(const char*)(conversion[i]))==0) {
            goal = i + 1;
        }

    if (goal==0) INPUT_ERROR

    file_name = argv[3];

    input_file = fopen(file_name,"r");

    /*scanning the file to get n and d:*/
    while (fscanf(input_file, "%s", line)>0 && *line != EOF){
        if (n==0) {
            i = 0;
            d=1;
            while (*(line + i) != '\0'){
                if (*(line + i) == ',')
                    d += 1;
                i++;
            }
        }
        n+=1;
    }
    fclose(input_file);

    /*Allocate memory  :*/
    dots = (double**)get_new_matrix(n,d);
    matrix =  (double**)get_new_matrix(n,n);
    array = (double*)malloc(n*sizeof(double));
    if(!array) ERROR;

    /* extract the data from the file to the dots array: */
    input_file2 = fopen(file_name,"r");
    for(i = 0; i < n && *line != EOF; i++) {
        fscanf(input_file2, "%s", line);
        tokens = strtok(line,",");
        j=0;
        while(tokens!= NULL){
            dots[i][j] = atof(tokens);
            tokens = strtok(NULL,",");
            j++;
        }
    }
    fclose(input_file2);

    /* applying the algorithm: */
    apply_requested_algo(goal, dots ,matrix, array, n,d,k,1);

    /* frees */
    for(i=0;i<n;i++) {
        free(matrix[i]);
        free(dots[i]);
    }
    free(matrix);
    free(dots);
    free(array);

    return 0;
}

/*------------------------- applying the requested algorithm according to the the goal: ------------------------------*/

static int apply_requested_algo(int goal, double **dots,double **matrix, double* array, int n, int d, int k, int print) {
    int i, j;
    switch (goal) {
        case 1:
            k = spk_goal(dots, matrix, array, n, d, k);
            if (print)
                print_2d_array(matrix, k, k);
            break;

        case 2:
            wam_goal(dots, matrix, n, d);
            if (print)
                print_2d_array(matrix, n, n);
            break;

        case 3:
            ddg_goal(dots, matrix, array, n, d);
            convert_diagonal_to_matrix(matrix, array, n);
            print_2d_array(matrix, n, n);
            break;

        case 4:
            lnorm_goal(dots, matrix, array, n, d);
            if (print)
                print_2d_array(matrix, n, n);
            break;

        case 5:
            jacobi_goal(dots, matrix, array, n);
            if (print) {
                for (i = 0; i < n; i++) {
                    if (array[i] <= 0 && array[i] > -0.00005)
                        printf("%.4f", 0.0000);
                    else
                        printf("%.4f", array[i]);
                    if (i != n - 1)
                        printf(",");
                    else
                        printf("\n");
                }
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++) {
                        if (matrix[j][i] <= 0 && matrix[j][i] > -0.00005)
                            printf("%.4f", 0.0000);
                        else
                            printf("%.4f", matrix[j][i]);
                        if (j != n - 1)
                            printf(",");
                        else
                            printf("\n");
                    }
                break;
            }
    }
    return k;
}


/*----------- wam - Weighted Adjacency Matrix calc function: -----------
 * assign the wam matrix of dots in matrix */

static void wam_goal(double **dots, double **matrix, int n, int d){
    int i, j;
    for(i=0; i< n; i++) {
        for (j = 0; j < i+1; j++)
            if (i == j)
                matrix[i][j] = 0;
            else {
                matrix[i][j] = exp(-((double) sqrt(calc_distance(dots[i], dots[j], d))) / 2);
                matrix[j][i]=matrix[i][j];
            }
    }
}

/* ----------- ddg -  Diagonal Degree Matrix: -----------*/
/* Assign the W matrix into matrix
 * and calc an array with the diagonal of the D matrix in array */

static void ddg_goal (double **dots, double **matrix, double* array, int n, int d){
    /* variables declaration*/
    int i, j;

    /*calc wam*/
    wam_goal(dots, matrix, n, d);

    /*calc D in array*/
    for(i=0; i<n; i++) {
        array[i] = 0;
        for (j = 0; j < n; j++)
            array[i] += matrix[i][j];
    }
}

/*----------- lnorm -  Normalized Graph Laplacian calc function: -----------
 * assign the I-D^(-1/2)WD^(-1/2) matrix in the matrix pointer */

static void lnorm_goal (double **dots, double **matrix, double *array, int n, int d){
    int i, j;
    ddg_goal(dots, matrix,array, n, d);
    /* assign I-D^(-1/2)WD^(-1/2) in matrix*/
    for(i=0; i<n;i++)
        for(j=0; j<i+1; j++) {
            /* calc D^(-1/2)WD^(-1/2)*/
            matrix[i][j] *= 1/sqrt(array[i]);
            matrix[i][j] *= 1/sqrt(array[j]);
            /*calc I-D^(-1/2)WD^(-1/2)*/
            matrix[i][j]*=-1;
            if(i==j)
                matrix[i][j]+=1;
            matrix[j][i] = matrix[i][j];
        }
}


/*----------- convert_diagonal_to_matrix: -----------*/
/* assign the diagonal in array in matrix as a diagonal matrix*/

static void convert_diagonal_to_matrix(double **matrix, double *array, int n){
    int i, j;
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            if(i==j)
                matrix[i][j]=array[i];
            else
                matrix[i][j]=0.0;
}

/*----------- jacobi_goal function: -----------*/
/* calc the eigenvectors (of L_norm of W of dots) in matrix
 * and a list of the eigenvalues in array */

static void jacobi_goal (double **dots, double **matrix, double* array, int n) {

    /*variables declaration*/
    int cnt=0, i, j;
    double off_new_matrix=0,off_matrix=0, diff=0.1;
    struct p_matrix p_matrices[MAX_JACOBI_ITER];

    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            matrix[i][j] = dots[i][j];

    /*calculate the off^2 of matrix */
    for(i=0; i<n; i++)
        for(j=i+1; j<n; j++)
            off_matrix+=matrix[i][j]*matrix[i][j];
    off_matrix*=2;

    /* evaluate the Jacobi method to find eigenvalues and the P's matrices*/
    while (cnt < MAX_JACOBI_ITER && diff > pow(10,-15)) {
        off_new_matrix =  jacobi_step(cnt, p_matrices, matrix ,n);
        diff = off_matrix - off_new_matrix;
        off_matrix = off_new_matrix;
        cnt+=1;
    }

    /*  assign the eigenvalues in a list*/
    for (i = 0; i < n; i++)
        array[i] = matrix[i][i];

    /*calc the matrices of the eigenvectors into matrix*/
    calc_eigen_vector_matrix(cnt, p_matrices, matrix,n);
}

/* ----------- jacobi_step: -----------
 * assign in the given matrix (size n),
 * a new matrix with values of P^t*matrix*P
 * and return the off(A)^2 of the new matrix*/

static double jacobi_step(int cnt, struct p_matrix p_matrices[], double** matrix, int n){
    /* variables declaration */
    int i, j, ind_i=0, ind_j=1, sign_theta;
    double theta, t, c, s, biggest_off_diagonal=matrix[0][1], sum_sqr_matrix=0.0;

    /* calculation necessary values */
    for(i=0; i < n; i++)
        for (j=0 ; j<n ; j++){
            if(i!=j && fabs(matrix[i][j])>biggest_off_diagonal){
                biggest_off_diagonal = fabs(matrix[i][j]);
                ind_i = i;
                ind_j =  j;
            }
        }

    theta = (matrix[ind_j][ind_j]-matrix[ind_i][ind_i])/(2*matrix[ind_i][ind_j]);
    sign_theta = (theta>=0)? 1: -1;
    t = (theta==0)? 1: sign_theta/(fabs(theta)+sqrt(theta*theta+1));
    c = 1/sqrt(t*t+1);
    s = t*c;

    /*assign values in the p_matrices*/
    p_matrices[cnt].i = ind_i;
    p_matrices[cnt].j = ind_j;
    p_matrices[cnt].s = s;
    p_matrices[cnt].c = c;

    /*assign values of P^t*Matrix*P in matrix */
    calc_p_transpose_multi_matrix_multi_p(matrix,p_matrices[cnt], n);

     /* calc sum of square of off diagonal of new matrix*/
    for(i=0; i<n; i++)
        for(j=i+1; j<n; j++)
                sum_sqr_matrix+=matrix[i][j]*matrix[i][j];

    return sum_sqr_matrix*2;
}

/* ----------- calc_p_transpose_multi_matrix_multi_p: -----------
 * assign in matrix a matrix of p^t*matrix*p */

static void calc_p_transpose_multi_matrix_multi_p(double **matrix,struct p_matrix p, int n){

    double *col_i, *col_j;
    double matrix_i_j = matrix[p.i][p.j], matrix_j_j=matrix[p.j][p.j], matrix_i_i=matrix[p.i][p.i];
    int row;
    col_i = (double*)malloc(n*sizeof(double));
    col_j = (double*)malloc(n*sizeof(double));

    for(row=0; row<n; row++){
        col_i[row] = matrix[row][p.i];
        col_j[row] = matrix[row][p.j];
    }

    for(row = 0; row<n; row++){
        matrix[row][p.i] = p.c * col_i[row] - p.s * col_j[row];
        matrix[p.i][row] = matrix[row][p.i];
        matrix[row][p.j] = p.c * col_j[row] + p.s * col_i[row];
        matrix[p.j][row] = matrix[row][p.j];
    }
    matrix[p.i][p.i] = p.c * p.c * matrix_i_i + p.s * p.s * matrix_j_j - 2 * p.c * p.s * matrix_i_j;
    matrix[p.j][p.j] = p.s * p.s * matrix_i_i + p.c * p.c * matrix_j_j + 2 * p.c * p.s * matrix_i_j;
    matrix[p.i][p.j] = 0;
    matrix[p.j][p.i] = 0;

    free(col_i);
    free(col_j);
}

/* ----------- calc_eigen_vector_matrix: -----------
 * assign in the given matrix (size n),
 * a matrix of eigenvectors in the colomns from array of p's matrices (size cnt) */

static void calc_eigen_vector_matrix(int cnt, struct p_matrix p_matrices[], double **matrix, int n){
    int i,j;

    /*assign p1 in matrix*/
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            if (i == j) {
                if (i == p_matrices[0].i || i == p_matrices[0].j)
                    matrix[i][j] = p_matrices[0].c;
                else
                    matrix[i][j] = 1;
            } else {
                if (i == p_matrices[0].i && j == p_matrices[0].j)
                    matrix[i][j] = p_matrices[0].s;
                else if (i == p_matrices[0].j && j == p_matrices[0].i)
                    matrix[i][j] = - p_matrices[0].s;
                else
                    matrix[i][j] = 0;
            }

    /* calc p_1*p_2*...*p_cnt*/
    for(i=1; i<cnt; i++)
        calc_matrix_multi_p(matrix, p_matrices[i], n);
}

/* ----------- calc_matrix_multi_p: -----------
 * assign in the given matrix (size n),
 * a matrix of the duplicity of matrix with the p given matrix */

static void calc_matrix_multi_p(double **matrix, struct p_matrix p, int n) {
    double temp_row_i, temp_row_j;
    int row;
    for(row=0;row<n;row++) {
        temp_row_i = matrix[row][p.i];
        temp_row_j  = matrix[row][p.j];
        matrix[row][p.i] = p.c * temp_row_i - p.s * temp_row_j;
        matrix[row][p.j] = p.c * temp_row_j + p.s * temp_row_i;
    }
}

/*----------- spk - Full spectral kmeans calc function: -----------
 * assigning in matrix the k centroids in dimension k
 * and in array the clasters indexes of the n points */

static int spk_goal(double **dots, double **matrix, double *array, int n, int d, int k){
    double** centroids;
    double sum_of_sqr, temp;
    int i, j;

    /*assign A=lnorm  in matrix*/
    lnorm_goal(dots, matrix, array, n, d);
    jacobi_goal (matrix, matrix, array, n);
    sort_eigenvalues_and_eigenvectors(array, matrix, n);

    if(k==0)
        k = calc_k(array, n);

    /*normalized the matrix U - k first eigenvectors in the columns of matrix */
    for(i=0; i<n;i++) {
        sum_of_sqr = 0;
        for (j = 0; j < k; j++)
            sum_of_sqr += matrix[i][j] * matrix[i][j];
        temp = sqrt(sum_of_sqr);
        sum_of_sqr = temp;
        for(j=0;j<k;j++)
            matrix[i][j]/=sum_of_sqr;
    }

    centroids =  (double**)get_new_matrix(k,k);

    /* assigning the k first vectors as initialized centroids*/
    for(i=0; i<k;i++)
        for(j=0;j<k;j++)
            centroids[i][j] = matrix[i][j];

    /*applying the kmeans algorithem on the n rows of matrix as a nXk matrix*/
    kmeans_calc(matrix,centroids, array, n,k,k);

    /*assigning final centroids in matrix*/
    for(i=0; i<k; i++){
        for(j=0;j<k; j++)
            matrix[i][j] = centroids[i][j];
    }

    /*frees*/
    for(i=0;i<k;i++)
        free(centroids[i]);
    free(centroids);

    return k;
}

/*-------- sort_eigenvalues_and_eigenvectors--------
 * sort the eigenvalues in array and the columns of matrix in respect of the eigenvalues*/

static void sort_eigenvalues_and_eigenvectors(double *array, double **matrix, int n){
    int i, j, m;
    double temp;
    for(i=0;i<n;i++)
        for(j=i+1; j < n ;j++)
            if(array[j] < array[i]) {
                temp = array[i];
                array[i] = array[j];
                array[j] = temp;
                for(m=0; m < n ; m++){
                    temp = matrix[m][j];
                    matrix[m][j] = matrix[m][i];
                    matrix[m][i] = temp;
                }
            }
}

/*----------- calc_k : -----------
 * calc c by finding the maximum delta between the eigen_values*/
static int calc_k(double *eigen_values, int n){
    int i, k = 0;
    double max_delta = 0.0;
    for(i=0; i<floor(n/2);i++)
        if(fabs(eigen_values[i]-eigen_values[i+1]) > max_delta) {
            max_delta = fabs(eigen_values[i]-eigen_values[i+1]);
            k=i;
        }
    return k+1;
}

/*----------- kmeans_calc function: -----------*/
/* update centroids with final k means centroids
 * and dots_location[i] with the index of the final clusters of the dot i*/

static void kmeans_calc(double **dots,double **centroids, double *dots_location, int n, int d, int k){

    /* variables declaration: */
    int i, j, iter_counter = 1; /*try for kmeans++ parpose*/
    double **old_centroids;

    /* applying the algorithm: */
    /* Initialize a centroids matrix, and a dots cluster location array :*/
    old_centroids =  (double**)get_new_matrix(k,d);

    /* the algorithm flow :*/
    while(iter_counter < MAX_ITER){
        for(i=0; i<n; i++) {
            dots_location[i] = find_nearest_centroid(dots[i], centroids, k, d);
        }
        for(i=0; i<k; i++)
            for(j=0;j<d; j++)
                old_centroids[i][j]=centroids[i][j];

        update_centroids(dots, dots_location, centroids, n, d, k);

        if(check_equals_2d_list(old_centroids, centroids, k, d))
            break;

        iter_counter+=1;
    } /* end of while */

    /* frees */
    for(i=0; i<k; i++)
        free(old_centroids[i]);
    free(old_centroids);
}

/*----------- calc_distance: -----------
 * calc distance between two d dimension vectors*/

static double calc_distance(double *dot, double *centroid, int d){
    double sum = 0;
    int i;
    for(i=0 ; i<d ; i++)
        sum += ((dot[i]-centroid[i])*(dot[i]-centroid[i]));
    return sum;
}

/*----------- find_nearest_centroid: -----------*/
/* returns the index of the nearest_centroid of the dot*/

static int find_nearest_centroid(double *dot, double **centroids, int k, int d){
    double min_distance = calc_distance(dot, centroids[0], d);
    double distance;
    int min_index = 0, i;
    for (i = 1 ; i < k ; ++i){
        distance = calc_distance(dot, *(centroids+i), d);
        if (distance < min_distance){
            min_distance = distance;
            min_index = i;
        }
    }
    return min_index;
}

/*----------- check_equals_2d_list: ----------- */
/*compare between two 2 dimensions lists of the same dimensions:*/

static int check_equals_2d_list(double **list1, double **list2, int row_num, int col_num){
    int i, j;
    for(i = 0; i < row_num; i++) {
        for (j = 0; j < col_num; j++) {
            if (list1[i][j] != list2[i][j])
                return 0;
        }
    }
    return 1;
}

/*----------- update_centroid: -----------*/

static void update_centroids(double **dots,double *dots_location, double **centroids, int n, int d, int k) {

    /* space allocation */
    int i, j;
    int *num_of_dots_per_cluster;
    double **sum_of_dots_per_cluster =  (double**)get_new_matrix(k,k);

    num_of_dots_per_cluster = (int*) malloc(k * sizeof(int));
    if(!num_of_dots_per_cluster)  ERROR


    /* initializing num_of_dots_per_cluster */
    for (i = 0; i < k; i++)
        num_of_dots_per_cluster[i] = 0;

    /* assignments in the sum and num arrays */
    for (i=0; i<n; i++){
        num_of_dots_per_cluster[(int)dots_location[i]] += 1;
        for(j=0; j<d; j++)
            sum_of_dots_per_cluster[(int)dots_location[i]][j] += dots[i][j];
    }

    /* assignments in the centroids array */
    for (i=0; i<k; i++)
        for(j=0; j<d; j++)
            centroids[i][j] = (sum_of_dots_per_cluster[i][j])/num_of_dots_per_cluster[i];

    /* free space */
    free(num_of_dots_per_cluster);
    for(i = 0; i<k; i++)
        free(sum_of_dots_per_cluster[i]);
    free(sum_of_dots_per_cluster);
}

/* ----------- print_2d_array: -----------*/

static void print_2d_array(double **array, int row_num, int col_num){
    int i, j;
    for(i=0; i < row_num; i++)
    {
        for(j = 0; j < col_num; j++) {
            if (array[i][j]<= 0 && array[i][j] > -0.00005)
                printf("%.4f",0.0000);
            else
                printf("%.4f", array[i][j]);
            if(j!= col_num-1)
                printf(",");
        }
        printf("\n");
    }
}

/*------get_new_matrix-----
 * return pointer to n*d zeros matrix */
static double** get_new_matrix(int n, int d){
    int i;
    double **matrix = (double**) malloc(n*sizeof(double*));
    if(!matrix) ERROR
    for(i=0; i<n; i++){
        matrix[i] = (double*)calloc(d, sizeof(double));
        if(!matrix[i]) ERROR
    }
    return matrix;
}

