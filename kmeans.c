#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>


static void KMeansClustering(int, int, int, int, double *, double *);


//Calc Euc distance
static double EuclideanDistance(double *x1 , double *x2, int d){
    double res = 0;
    int i;
    for (i = 0; i < d; i++) {
        res = res + (x1[i]-x2[i])*(x1[i]-x2[i]); }
    return res;
}

//Make a new matrix (vectors with each entry being a vector)
static double** makeANewMatrix(int p, int n) {
    int i;
    double **newMat = (double**)calloc(p, sizeof(double*));
    assert(newMat != NULL);

    for(i = 0; i < p; ++i) {
        newMat[i] = (double*)calloc(n, sizeof(double));
        assert(newMat[i] != NULL);
    }
    return newMat;
}

//find min Euc distance
static int FindMinDistance(double *x, double **vecs, int k, int d){
    int i=0;
    int resIndex = 0;
    double resDistance = 0;
    double distance = 0;
    resDistance = EuclideanDistance(x, vecs[0], d);
    for (i=1; i<k; i++){
        distance = EuclideanDistance(x, vecs[i], d);
        if (distance < resDistance){
            resIndex = i;
            resDistance = distance;
        }
    }
    return resIndex;
}

//fill the matrix with proper values from python
static void fillMatrix(int N, int d, double **vecs, PyObject *vectors_Obj) {
    
    PyObject *thisVector;
    PyObject *thisVal;


    int i, j;

    for (i = 0; i < N; i++) {
        thisVector = PyList_GetItem(vectors_Obj, i);
        assert(PyList_Check(thisVector));
        
        for (j = 0; j < d; j++)
        {
            thisVal = PyList_GetItem(thisVector, j);
            vecs[i][j] = PyFloat_AsDouble(thisVal);
        }
    }
}

//free all memory from matrix
static void freeMatFromMem(double **mat, int length) {
    int j;
    for (j = 0; j < length; ++j)
    {
        free(mat[j]);
    }

    free(mat);
}



//the algorithm
static PyObject* fit(PyObject *self, PyObject *args){
    
    double **vectors, **centroids;
    double *vectorsB, *centroidsB;
    int k, N, d, MAX_ITER;
    

    PyObject *PyVecs;
    PyObject *PyStartCents;

    if (!PyArg_ParseTuple(args, "O!O!iiii", &PyList_Type, &PyVecs, &PyList_Type, &PyStartCents, &k, &d, &MAX_ITER,&N))
        return NULL;
    
    

    //initiating all the needed pointers and lists
    centroidsB=malloc(k * d * sizeof(double));
    vectorsB=malloc(N * d * sizeof(double));

    int i, j;

    vectors = makeANewMatrix(N, d);
    fillMatrix(N, d, vectors, PyVecs);

    centroids = makeANewMatrix(k, d);
    fillMatrix(k, d, centroids, PyStartCents);

    //making the matrix one big array in order to use our old kmeans
    for (i = 0; i < N; i++) {
        for (j = 0; j < d; j++)
        {
            vectorsB[i*d+j] = vectors[i][j];
            if(i<k){
                centroidsB[i*d+j]=centroids[i][j];
            }
        }
    }






    //HERE WE HAVE ALL IN C ONLY COMPUTE!
    KMeansClustering(k, N, d, MAX_ITER, vectorsB, centroidsB);

    PyObject *PyAllRows;
    PyObject *PyOneRow;
    PyObject *PyFloat;

    PyAllRows = PyList_New(N);
    for(i = 0; i < N; i++){
		PyOneRow = PyList_New(d);

		for(j = 0; j < d; j++){
			PyFloat = Py_BuildValue("d", centroidsB[i*d+j]);

			PyList_SetItem(PyOneRow, j, PyFloat);
		}
		PyList_SetItem(PyAllRows, i, PyOneRow);
	}

    //free memory
    freeMatFromMem(vectors,N);
    freeMatFromMem(centroids,k);
    return PyAllRows;

}


static PyMethodDef capiMethods[] = {
        {"fit",
                (PyCFunction) fit,
                     METH_VARARGS,
                PyDoc_STR("Sort N vectors to their k fitting clusters")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        capiMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void){
    return PyModule_Create(&moduledef);
}

//the old hw1 algorithm
static void KMeansClustering(int k, int N, int d, int MAX_ITER, double* vectors,double* centroids){
    /*stating variables and arrays*/
    double **centroids_p;
    double **vectors_p;
    double *new_cents;
    double **new_cents_p;
    int *sNums;
    int centChanged, clustNum, i, j;
    int currIter = 0;
    double epsilon = 0.000001;

    
    



    /*init arrays*/
    sNums = calloc(k, sizeof(double));
    new_cents = calloc(k * d, sizeof(double));
    new_cents_p = calloc(k, sizeof(double*));
    assert(new_cents != NULL);
    assert(new_cents_p != NULL);
    for (i = 0; i < k ; i++ ){
        new_cents_p[i] = new_cents + i * d;}

    //centroids = malloc(k * d * sizeof(double));
    centroids_p = calloc(k,sizeof(double*));
    for(i = 0; i < k; i++)
        centroids_p[i] = centroids + i * d;

    //centroids = malloc(k * d * sizeof(double));
    vectors_p = calloc(N,sizeof(double*));
    for(i = 0; i < N; i++)
        vectors_p[i] = vectors + i * d;


    //finished init of arrays

    


    

    //compute each observation's cluster and calculate new centroids*/
    centChanged = 1;
    while (centChanged == 1  && currIter < MAX_ITER){
        centChanged = 0;
        //for each observation, find closest cluster*/

        for (i = 0; i < N; i++) {
            clustNum = FindMinDistance(vectors_p[i], centroids_p, k, d);
            for (j = 0; j < d; j++) {
                new_cents_p[clustNum][j] = (new_cents_p[clustNum][j] * sNums[clustNum] + vectors_p[i][j])
                                           / (sNums[clustNum]+1);
            }
            sNums[clustNum] = sNums[clustNum] + 1;
        }

        //check if centroid group changed*/
        for (i = 0; i < k; ++i) {
            double dist = (EuclideanDistance(centroids_p[i], new_cents_p[i], d));
            if (dist > epsilon && centChanged == 0){
                centChanged = 1;
            }
        }
        currIter++;
        for (i=0; i < k; ++i){
            for (j=0; j<d; ++j){
                centroids_p[i][j] = new_cents_p[i][j];
                new_cents_p[i][j] = 0;
            }
        }


        for (i=0; i < k; i++){
            sNums[i] = 0;
        }
    }
    free(centroids_p);
    free(vectors);
    free(sNums);
    free(vectors_p);
    free(new_cents);
    free(new_cents_p);
}


