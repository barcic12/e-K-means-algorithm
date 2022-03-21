#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int KMeansClustering(int, int);
double EuclideanDistance(double [], double [], int);
int FindMinDistance(double [], double* [], int k, int d);


int KMeansClustering(int k, int MAX_ITER){
    double *centroids;
    double **centroids_p;
    double *observations;
    double **observations_p;
    double *s;
    double **S_p;
    int *sNums;
    int centChanged, clustNum, i, j;
    int currIter = 0;
    char line_cnt[1000];
    int N=0;
    int d=1;
    char ch;

   
    while (!feof(stdin)){
    if (fgets( line_cnt,1000,stdin) != NULL){
        N++;
    }
}
rewind(stdin);
    while((ch = fgetc(stdin))!=EOF){
        if(ch=='\n'||ch=='\0'){
            break;
        }
        else if(ch==','){ 
            d++;
        }
    }

rewind(stdin);
assert(N>k);

    sNums = calloc(k, sizeof(double));
    s = calloc(k*d, sizeof(double));
    S_p = calloc(k, sizeof(double*));
    assert(s != NULL);
    assert(S_p != NULL);
    for (i = 0; i < k ; i++ ){
        S_p[i] = s + i * d;}

    centroids = malloc(k * d * sizeof(double));
    centroids_p = malloc(k * sizeof(double*));
    for( i=0 ; i < k ; i++ )
        centroids_p[i] = centroids+ i * d;

    observations = malloc(N * d * sizeof(double));
    observations_p = malloc(N * sizeof(double*));
    for( i=0 ; i < N ; i++ ){
        observations_p[i] = observations+ i * d;}

    i=0;
    while((scanf("%lf,", &observations[i]))==1) {
        i++;


    }
   
    for (i = 0; i < k*d; ++i) {
        centroids[i] = observations[i];
    }

    centChanged = 1;
    while (centChanged == 1  && currIter < MAX_ITER){
        centChanged = 0;

        for (i = 0; i < N; ++i) {

            clustNum = FindMinDistance(observations_p[i], centroids_p, k, d);

            for (j = 0; j < d; ++j) {
                S_p[clustNum][j] = (S_p[clustNum][j] * sNums[clustNum] + observations_p[i][j])
                        / (sNums[clustNum]+1);
            }
            sNums[clustNum] = sNums[clustNum] + 1;

        }


      
        for (i = 0; i < k; ++i) {
            if (EuclideanDistance(centroids_p[i], S_p[i], d) != 0.0 && centChanged == 0){

                centChanged = 1;
            }
        }
        currIter++;

        for (i = 0; i < k*d; ++i) {
            centroids[i] = s[i];

            s[i] = 0;

        }

        for (i=0; i < k; i++){
            sNums[i] = 0;

        }

    }

    for (i = 0; i < k; ++i) {
        for (j = 0; j < d - 1; ++j) {
            printf("%.4f,", centroids_p[i][j]);
        }
        printf("%.4f\n", centroids_p[i][d-1]);
    }
    free(centroids);
    free(centroids_p);
    free(observations);
    free(observations_p);
    free(sNums);
    free(S_p);
    free(s);
    return 1;
}



double EuclideanDistance(double *x1 , double *x2, int d){
    double res = 0;
    int i;
    for (i = 0; i < d; i++) {
        res = res + (x1[i]-x2[i])*(x1[i]-x2[i]); }
    return res;
}

int FindMinDistance(double *x, double **vecs, int k, int d){
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


int main(int argc, char* argv[]) {
    int i;
    long val;
    char *next;
    assert((argc == 3) || (argc == 2));
    val =0;
    if(val!=0){
        return 0;
    }

    if (argc==3) {
        for (i = 1; i < 3; i++) {
            val = strtol (argv[i], &next, 10);
             if ((next == argv[i]) || (*next != '\0')) {
                 printf("the argumnts are not in the right format\n");
                 return 0;
        }    
    }
     KMeansClustering(atoi(argv[1]), atoi(argv[2]));   
    }
    else {
        for (i = 1; i < 2; i++) {


            val = strtol (argv[i], &next, 10);

            if ((next == argv[i]) || (*next != '\0')) {
                 printf("the argumnts are not in the right format\n");
                 return 0;
            } 
        
        }
        KMeansClustering(atoi(argv[1]),  200);
    }

    return 0;

}
