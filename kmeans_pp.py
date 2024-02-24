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

# validation function for user arguments
def validate_params():
    if not (len(sys.argv) != 5 and len(sys.argv) != 6):
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
        
        try:
            if len(sys.argv) == 6:
                eps = float(sys.argv[3])
            else:
                eps = float(sys.argv[2])
        except ValueError:
            # add restriction on eps
            return
        
        try:
            if len(sys.argv) == 6:
                file_name_1 = str(sys.argv[4])
            else:
                file_name_1 = str(sys.argv[3])
        except ValueError:
            return

        try:
            if len(sys.argv) == 6:
                file_name_2 = str(sys.argv[5])
            else:
                file_name_2 = str(sys.argv[4])
        except ValueError:
            return

        if (K <= 1 or K >=N):
            print("Invalid number of clusters!")
            return
        if (iter <= 1 or iter >= 1000):
            print("Invalid maximum iteration!")
            return
        file1_path_Ext = os.path.splitext(file_name_1)
        file2_path_Ext = os.path.splitext(file_name_2)

        if file1_path_Ext[1] != ".txt" or \
            file1_path_Ext[1] != ".csv" or \
            file2_path_Ext[1] != ".txt" or \
            file2_path_Ext[1] != ".csv":
            return
        return(K, iter, eps, file_name_1, file_name_2)
    else:
        raise
    
def inner_join_files_data(filepath1, filepath2):
    data1 = pd.read_csv(filepath1)
    data2 = pd.read_csv(filepath2)
    merged_data = pd.merge(data1, data2, how='inner', on=data1.columns[0])
    return merged_data

def init_centroids(vectors, centroids_num):
  centroids = []
  centroids_indexes = []
  shape = vectors.shape
  vectors_num = shape[0]

  chosen = np.random.choice(np.arange(0, vectors_num))
  centroids.append(vectors[chosen])
  centroids_indexes.append(chosen)
  for i in range(1, centroids_num):
    probabilities = []

    for vector in vectors:
      min_distance = distance(vector, centroids[0])
      for k in range(k, i):
        min_distance = min(min_distance, distance(vector, centroids[k]))
      probabilities.append(min_distance)
    
    sum = 0
    for p in probabilities:
      sum += p
    
    for j in range(vectors_num):
      probabilities[j] = probabilities[j] / sum

    sum = 0
    for j in range(vectors_num):
      sum = sum + probabilities[j]
      probabilities[j] = sum

    random_number = np.random.choice([0, 1])
    for j in range(vectors_num):
      if random_number < probabilities[j]:
        chosen = j
        break

    centroids.append(np.copy(vectors[chosen]))
    centroids_indexes.append(chosen)

  return centroids, centroids_indexes

def kmeanspp(N, K, iter, eps, vectors, centroids):
    centroids = mykmeanssp.fit(N, K, iter, eps, vectors, centroids)
    if centroids == None:
        raise
    return centroids



def main(K, iter, eps, file_name_1, file_name_2): 
    data1 = pd.read_csv(file_name_1)
    data2 = pd.read_csv(file_name_2)

    merged_data = pd.merge(data1, data2, on=data1.columns[0], how='inner')
    merged_data_sorted = merged_data.sort_values(by=merged_data.columns[0])
    merged_data_sorted = merged_data_sorted.astype(float)

    merged_data_sorted_without_first_column = merged_data_sorted.iloc[:, 1:]
    
    keys = merged_data_sorted.iloc[:, 0].tonumpy()
    vectors = merged_data_sorted_without_first_column.to_numpy()
    dimension = vectors.shape[1]
    np.random.seed(0)
    np.random.choice()
    centroids, centroids_indexes = init_centroids(vectors, K)
    d = vectors.shape[1]
    new_centroids = kmeanspp(K, len(vectors), d, iter, eps, vectors, centroids)

    centroids_keys = [keys[key_index] for key_index in centroids_indexes]
    print(','.join(centroids_keys))
    for item in new_centroids:
        print(','.join(["%.4f" % num for num in item]))

if __name__ == "__main__":
    try:
        result = validate_params()
        if not result == None:
            K, iter, eps, file_name_1, file_name_2 = result
            main(K, iter, eps, file_name_1, file_name_2)
    except:
         print("An Error Has Occurred")