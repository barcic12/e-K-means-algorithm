import sys

def euclidean_distance(x1, x2):
    result = 0
    for i in range(len(x1)):
        result += (float(x1[i]) - float(x2[i])) ** 2
    return result


def min_distance(x, vec):
    result = 0
    res_distance = euclidean_distance(x, vec[0])
    for i in range(1, len(vec)):
        distance = euclidean_distance(x, vec[i])
        if distance < res_distance:
            result = i
            res_distance = distance
    return result


def k_means_clustering(k, max_iter):
    centroids = []
    observations = []
    line_list=[]
    i = 0

    while True:
        try:
            line = input()
            line_list = [float(x) for x in line.strip().split(",")]
            observations.append(line_list)
            if i < k:
                centroids.append(line_list)
            i += 1
        except EOFError:
            break

    change = True
    iter = 0
    N=i
    assert(N>k)
    d = len(line_list)
    while change and iter < max_iter:
        change = False
        new_cents = [[0 for y in range(d)] for x in range(k)]
        s_nums = [0 for i in range(k)]
        for i in range(N):
            cluster = min_distance(observations[i], centroids)
            for j in range(d):
                new_cents[cluster][j] = ((new_cents[cluster][j] * s_nums[cluster] +
                                                 observations[i][j])) / (s_nums[cluster] + 1)
            s_nums[cluster] = s_nums[cluster] + 1
        for i in range(k):
            if (euclidean_distance(centroids[i], new_cents[i]) != 0):
                change = True
                break
        centroids = new_cents[:]
        iter += 1

    for i in range(k):
        for j in range(d - 1):
            print(str("%0.4f" % centroids[i][j]) + ",", end="")
        print("%0.4f" % centroids[i][d - 1])


if (len(sys.argv) == 3):
    try:
        int(sys.argv[1])
        int(sys.argv[2])
        k_means_clustering(int(sys.argv[1]), int(sys.argv[2]))
    except ValueError:
        print("the argumnts are not in the right format")
elif (len(sys.argv) == 2):
    try:
        int(sys.argv[1])
        k_means_clustering(int(sys.argv[1]), 200)
    except ValueError:
        print("the argumnts are not in the right format")
else:
    print("error in the number of argument")











    

