import sys
import os
import pandas as pd
import numpy as np
import math
import mykmeanssp

def distance(vector1, vector2, d):
    sum = 0
    for i in range(d):
        sum+=(vector1[i]-vector2[i])**2
    return math.sqrt(sum)

def validate_params(N):
    try:
        K = int(sys.argv[1])
    except ValueError:
        print("Invalid number of clusters!")
        return
    
    iter = 300 
    if len(sys.argv) == 6:
        try:
            iter = int(sys.argv[2])
        except:
            print("Invalid maximum iteration!")
            return
    

    if len(sys.argv) == 6:
        eps = float(sys.argv[3])
    else:
        eps = float(sys.argv[2])
    
    if eps < 0:
      print("Invalid epsilon!")

    if (K <= 1 or K >=N):
        print("Invalid number of clusters!")
        return
    if (iter <= 1 or iter >= 1000):
        print("Invalid maximum iteration!")
        return
    return K, iter, eps
    
def validate_files():
    if len(sys.argv) == 6:
        file_name_1 = str(sys.argv[4])
    else:
        file_name_1 = str(sys.argv[3])

    if len(sys.argv) == 6:
        file_name_2 = str(sys.argv[5])
    else:
        file_name_2 = str(sys.argv[4])

    file1_path_Ext = os.path.splitext(file_name_1)
    file2_path_Ext = os.path.splitext(file_name_2)

    if (file1_path_Ext[1] != ".txt" and \
        file1_path_Ext[1] != ".csv") or \
        (file2_path_Ext[1] != ".txt" and \
        file2_path_Ext[1] != ".csv"):
        raise
    
    return file_name_1, file_name_2

def init_centroids(vectors, centroids_num):
  centroids = []
  centroids_indexes = []

  vectors_num = len(vectors)

  chosen = np.random.choice(np.arange(0, vectors_num))
  centroids.append(vectors[chosen])
  centroids_indexes.append(chosen)
  for i in range(1, centroids_num):
    probabilities = []

    for vector in vectors:
      min_distance = distance(vector, centroids[0], len(vector))
      for k in range(1, i):
        min_distance = min(min_distance, distance(vector, centroids[k], len(vector)))
      probabilities.append(min_distance)
    
    distances_sum = sum(probabilities)
    probabilities = [p / distances_sum for p in probabilities]

    chosen = np.random.choice(vectors_num, p=probabilities)

    centroids.append(vectors[chosen].copy())
    centroids_indexes.append(chosen)

  return centroids, centroids_indexes

def call_c_kmeans(K, N, d, iter, eps, vectors, centroids):
    centroids = mykmeanssp.fit(K, N, d, iter, eps, vectors, centroids)
    if centroids == None:
      raise
    return centroids

def get_vectors_from_files(file_name_1, file_name_2):
    data1 = pd.read_csv(file_name_1, header=None)
    data2 = pd.read_csv(file_name_2, header=None)

    merged_data = pd.merge(data1, data2, on=data1.columns[0], how='inner')

    merged_data_sorted = merged_data.sort_values(by=merged_data.columns[0])
    merged_data_sorted = merged_data_sorted.astype(float)

    merged_data_sorted_without_first_column = merged_data_sorted.iloc[:, 1:]
    vectors = merged_data_sorted_without_first_column.values.tolist()

    return vectors

def main(K, iter, eps, vectors): 
    dimension = len(vectors[0])
    np.random.seed(0)
    
    centroids, centroids_indexes = init_centroids(vectors, K)
    d = len(vectors[0])
    new_centroids = call_c_kmeans(K, len(vectors), d, iter, eps, vectors, centroids)
    print(','.join([str(c) for c in centroids_indexes]))

    for item in new_centroids:
        print(','.join(["%.4f" % num for num in item]))

if __name__ == "__main__":
    try:
        if len(sys.argv) != 5 and len(sys.argv) != 6:
            raise
        files_result = validate_files()
        if not files_result == None:
            file_name_1, file_name_2 = files_result
            vectors = get_vectors_from_files(file_name_1, file_name_2)
            params_result = validate_params(len(vectors))
            if not params_result == None:
                K, iter, eps = params_result
                main(K, iter, eps, vectors)
    except:
         print("An Error Has Occurred")