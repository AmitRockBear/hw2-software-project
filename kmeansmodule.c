
#define _GNU_SOURCE
#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Python.h>

void free_array_of_pointers(double** arr, int length) {
  int i;
  for (i=0; i<length; i++) {
    free(arr[i]);
  }
  free(arr);
}

double calculate_distance(double* vec1, double* vec2, int size) {
  int i;
  double sum;
  
  sum = 0;
  for (i=0; i<size; i++) {
    sum += pow(vec1[i]-vec2[i], 2);
  }

  return sqrt(sum);
}

int find_closest_centroid_to_vector_index(double* vector, int vector_size, double** centroids, int centroids_num) {
  int j, min_j;
  double min_distance, distance_from_centroid;

  min_j = 0;
  min_distance = calculate_distance(vector, centroids[0], vector_size);

  for (j=0; j<centroids_num; j++) {
    distance_from_centroid = calculate_distance(vector, centroids[j], vector_size);
    if (min_distance > distance_from_centroid) {
      min_distance = distance_from_centroid;
      min_j = j;
    }
  }

  return min_j;
}

double* create_new_centroid(double* centroid_sum, int centroid_counter, int centroid_size) {
  int p;
  double *new_centroid;

  new_centroid = calloc(centroid_size, sizeof(double));
  if (new_centroid == NULL) {
      return NULL;
  }

  for (p=0; p<centroid_size; p++) {
    new_centroid[p] = centroid_sum[p] / centroid_counter;
  }

  return new_centroid;
}

int calculate_centroids_convergence(double** centroids, double** vectors, int centroids_num, int centroid_size, int vectors_num, int max_iterations, double eps) {
  int iter_couter, closest_centroid_index, i, j, p;
  double max_distance, centroids_distance, **centroids_sum, *counters, *new_centroid_j;
  max_distance = eps + 1;
  iter_couter = 0;

  while (max_distance >= eps && iter_couter < max_iterations) {
    max_distance = 0;
    centroids_sum = (double**)calloc(centroids_num, sizeof(double *));
    if (centroids_sum == NULL) {
      return 1;
    }
    counters = (double*)calloc(centroids_num, sizeof(double));
    if (counters == NULL) {
      free(centroids_sum);
      return 1;
    }

    for (j=0; j<centroids_num; j++) {
      centroids_sum[j] = calloc(centroid_size, sizeof(double));

      if (centroids_sum[j] == NULL) {
        free_array_of_pointers(centroids_sum, j);
        free(counters);
        return 1;
      }
    }

    for (i=0; i<vectors_num; i++) {
      closest_centroid_index = find_closest_centroid_to_vector_index(vectors[i], centroid_size, centroids, centroids_num);
      counters[closest_centroid_index]++;

      for (p=0; p<centroid_size; p++) {
        centroids_sum[closest_centroid_index][p]+=vectors[i][p];
      } 
    }

    for (j=0; j<centroids_num; j++) {
      if (counters[j] > 0) {
        new_centroid_j = create_new_centroid(centroids_sum[j], counters[j], centroid_size);

        if (new_centroid_j == NULL) {
          free_array_of_pointers(centroids_sum, centroids_num);
          free(counters);
          return 1;
        }
        
        centroids_distance = calculate_distance(centroids[j], new_centroid_j, centroid_size);
        if (centroids_distance > max_distance) {
          max_distance = centroids_distance;
        }
        free(centroids[j]);
        centroids[j]= new_centroid_j;
      }
    }

    free_array_of_pointers(centroids_sum, centroids_num);
    free(counters);

    iter_couter++;
  }

  return 0;
}

double** kmeans(int K, int N, int d, int iter, double eps, double** vectors, double** centroids) {
    int res;

    res = calculate_centroids_convergence(centroids, vectors, K, d, N, iter, eps);
    if (res == 1) {
      free_array_of_pointers(vectors, N);
      free_array_of_pointers(centroids, K);
      return NULL;
    }

    free_array_of_pointers(vectors, N);

    return centroids;
}

static PyObject* fit(PyObject *self, PyObject *args) {
    int K, N, d, iter;
    double eps, **vectors, **centroids;
    if(!PyArg_ParseTuple(args, "iiiiddoo", &K, &N, &d, &iter, &eps, &vectors, &centroids)) {
        return NULL;
    }

    return Py_BuildValue("o", kmeans(K, N, d, iter, eps, vectors, centroids));
}

static PyMethodDef mykmeansspMethods[] = {
    {
      "fit",
      (PyCFunction) fit,
      METH_VARARGS,
      PyDoc_STR(
        "Perform some operation using K-means algorithm. \
        Input: \
        K : int - The number of centroids. \
        N : int - The number of data points. \
        d : int - The dimensionality of the data points. \
        iter : int - The maximum number of iterations for the algorithm. \
        eps : float - The threshold for convergence. \
        vectors : double** - A 2D list representing the input data points. Each inner list represents a data point and should have 'd' elements. \
        centroids : double** - A 2D list representing the initial centroids for the clusters. It should have 'K' inner lists, each representing a centroid point and should have 'd' elements."
      )
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mykmeansspModule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    mykmeansspMethods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&mykmeansspModule);
    if (!m) {
        return NULL;
    }
    return m;
}

