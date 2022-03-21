import sys
import numpy as np
import spkmeansmodule
import math

np.set_printoptions(formatter={'float': lambda x: "{0:0.4f}".format(x)})
np.random.seed(0)


# kmeans_pp function

def kmeans_pp(points, k):
    np.random.seed(0)
    k_init_centroids = np.empty((k, points.shape[1]))
    k_init_centroids[0, :] = points[np.random.choice(points.shape[0]), :]  # choose random row
    d_values = np.empty(points.shape[0])  # temp array-like to contain d value for each point during the iteration
    for z in range(k - 1):  # loop z=1 until k
        d_temp = np.empty((z + 1, points.shape[1]))  # temp matrix-like to contain points[i]-centroid[i]
        for i in range(points.shape[0]):  # calculate d'i for each point
            d_temp = points[i, :] - k_init_centroids[0:z + 1, :]  # using broadcasting
            d_temp **= 2
            d_values[i] = np.min(np.sum(d_temp[:, 1:], axis=1))
        d_values /= np.sum(d_values, axis=0)
        k_init_centroids[z + 1, :] = points[np.random.choice(points.shape[0], p=d_values), :]
    return k_init_centroids


# print_matrix in the demanding format

def print_matrix(matrix):
    output_str = ""
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if (matrix[i][j]<=0) and (matrix[i][j]>-0.00005):
                output_str += "{:.4f}".format(0)
            else:
                output_str += "{:.4f}".format(matrix[i][j])
            if j != (matrix.shape[1] - 1):
                output_str += ','
        output_str+= '\n'
    print(output_str, end = "")


# main:

# getting the arguments:
assert (len(sys.argv) == 4), "Invalid Input!"

goal_dict = {"spk": 1, "wam" : 2 , "ddg" : 3, "lnorm" : 4, "jacobi" : 5}

# define k and goal
k = 0
goal = None
try:
    k = int(sys.argv[1])
    goal = goal_dict[sys.argv[2]]
except ValueError:
    print("Invalid Input!")

# get the files
points = np.genfromtxt(sys.argv[3], delimiter=',')
if(len(points.shape)==1):
    temp = np.zeros((points.shape[0],1))
    for i in range(points.shape[0]):
        temp[i][0] = points[i]
    points=temp
n = points.shape[0]

matrix = np.zeros((n,n))

## goal != spk
if goal in range(2,6):
    # sending to c to implement the requested algorithem corresponding to goal
    final_matrix = spkmeansmodule.fit(points.tolist(), matrix.tolist(), goal)
    py_final_matrix = np.array(final_matrix)
    if goal!=3 and goal!=5:
        print_matrix(py_final_matrix)
    if(goal == 5):
        for i in range (n):
            print("{:.4f}".format(py_final_matrix[0][i]), end = "")
            if(i!=n-1):
                print(",", end = "")
        print("")
        print_matrix(py_final_matrix[1:n+1,:].transpose())
## goal == spk
else:
    # sending to c to implement the lnorm on dots
    final_lnorm_matrix = spkmeansmodule.fit(points.tolist(), matrix.tolist(), 4)
    py_lnorm_matrix = np.array(final_lnorm_matrix)
    # sending to c to implement the jacobi on the result of lnorm
    final_jacobi_matrix = spkmeansmodule.fit(py_lnorm_matrix.tolist(), matrix.tolist(), 5)
    py_jacobi_matrix = np.array(final_jacobi_matrix)

    #sort the jacobi matrix
    for i in range(n):
        for j in range(i+1,n):
            if py_jacobi_matrix[0][j] < py_jacobi_matrix[0][i]:
                for m in range (n+1):
                    temp = py_jacobi_matrix[m][i]
                    py_jacobi_matrix[m][i] = py_jacobi_matrix[m][j]
                    py_jacobi_matrix[m][j] = temp

    #calc k if k==0
    if (k==0):
        max_diff = abs(py_jacobi_matrix[0][0]-py_jacobi_matrix[0][1])
        k=1
        for i in range(1,math.floor(n/2)):
            if abs(py_jacobi_matrix[0][i]-py_jacobi_matrix[0][i+1])>max_diff:
                max_diff = abs(py_jacobi_matrix[0][i]-py_jacobi_matrix[0][i+1])
                k = i+1;


    # leave the k first colomns and delete row of eigenvalues
    k_vector_matrix = py_jacobi_matrix[1:n+1,0:k]


    # normalization
    for i in range(n):
        sum = 0
        for j in range(k):
            sum = sum + k_vector_matrix[i][j]*k_vector_matrix[i][j]
        sum = math.sqrt(sum)
        for j in range(k):
            k_vector_matrix[i][j] = k_vector_matrix[i][j]/sum

    ## add index colomn
    vector_matrix_with_ind = np.zeros((n,k+1))
    for i in range (n):
        for j in range(k + 1):
            if j==0:
                vector_matrix_with_ind[i][0] = i
            else:
                vector_matrix_with_ind[i][j] = k_vector_matrix[i][j-1]

    ## find the kmeans++ k initial centroids
    matrix = kmeans_pp(vector_matrix_with_ind,k)

    # sending to c to implement the kmeans on the k first eigenvectors in matrix
    final_centroids = spkmeansmodule.fit(k_vector_matrix.tolist(), matrix[:,1:k+1].tolist(), 6)
    py_final_centroids = np.array(final_centroids)

    # output the data
    # first line
    for i in range (k):
        print(f"{int(matrix[i][0])}",end="")
        if i!=k-1:
            print(",",end="")
    print("\n", end="")

    # second line
    print_matrix(py_final_centroids)
