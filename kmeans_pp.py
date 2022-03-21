import mykmeanssp
import numpy as np
import pandas as pd
import sys

def minChance(prob, d):
    if(len(prob) == 0):
        return d
    for i in range(len(prob)):
        prob[i] = min(prob[i], d[i])
    return prob

def calculate_probabilty(prob, Sum):
    lst = []
    for i in range(len(prob)):
        lst.append(prob[i]/Sum)
    return lst

def toCommaSeparatedString(vec):
    return ",".join([repr(comp) for comp in vec])

def printCentroids(cents,k):
    rounded = np.round(cents, 4)
    i=0
    for vec in rounded:
        if i==k:
            break
        i+=1
        print(toCommaSeparatedString(vec))     



def merge(file1, file2):
    table_left = pd.read_csv(file1, header=None)
    table_right = pd.read_csv(file2, header=None)
    table = pd.merge(table_left, table_right, on=0)
    table = table.sort_values(by=[0], ascending=True)
    table = table.drop(columns=[0])

    N = table.shape[0] # N = number od data points given
    D = table.shape[1] #dimensions of the vectors

    vectors = table.to_numpy()
    return N,D,vectors

def computeInitialCentroids(K, N, vectors):
    np.random.seed(0)
    arr = np.arange(N) #uses only to get a random int - [0-N]
    selectedIndices = []
    initialCentroids = []

    #the algorithm
    rand = np.random.choice(arr)
    selectedIndices.append(rand)
    initialCentroids.append(vectors[rand])

    min_distances = []
    for i in range(K-1):
        d = []
        
        for vector in vectors:
            d.append(np.sum(np.power((vector - initialCentroids[i]), 2), axis=0)) 
        min_distances = minChance(min_distances, d) 
        probability_list = calculate_probabilty(min_distances, sum(min_distances)) 
        rand = np.random.choice(arr, p = probability_list)
        initialCentroids.append(vectors[rand])
        selectedIndices.append(rand)
    return selectedIndices, initialCentroids
    


def start(K,max_iter,file1,file2):
    N, D, vectors = merge(file1, file2)
    if(K >= N):
        raise("number of data points must be greater than K")

    selected_indexes, initialCentroids = computeInitialCentroids(K, N, vectors)

    initialCentroids = [cen.tolist() for cen in initialCentroids] #convert all centroids from np to list for API
    vectorsList = vectors.tolist() #convert vectors to list item to be consumed by C API
    kmeansResult = mykmeanssp.fit(vectorsList, initialCentroids, K, D, max_iter,N)

    print(toCommaSeparatedString(selected_indexes))
    printCentroids(kmeansResult,K)

if (__name__ == "__main__"):
    if (len(sys.argv) == 5):
        try:
            int(sys.argv[1])
            int(sys.argv[2])
            start(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4])
        except ValueError:
            print("the arguments are not in the right format")
    elif (len(sys.argv) == 4):
        try:
            int(sys.argv[1])
            start(int(sys.argv[1]),300, sys.argv[2], sys.argv[3])
        except ValueError:
            print("the arguments are not in the right format")
    else:
        print("error in the number of arguments")
