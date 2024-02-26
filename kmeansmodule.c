
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
  printf("calculate_distance starting");
  
  sum = 0;
  for (i=0; i<size; i++) {
    sum += pow(vec1[i]-vec2[i], 2);
  }
  printf("calculate_distance ending");

  return sqrt(sum);
}

int find_closest_centroid_to_vector_index(double* vector, int vector_size, double** centroids, int centroids_num) {
  int j, min_j;
  double min_distance, distance_from_centroid;
  printf("find_closest_centroid_to_vector_index starting");

  min_j = 0;
  min_distance = calculate_distance(vector, centroids[0], vector_size);

  for (j=0; j<centroids_num; j++) {
    distance_from_centroid = calculate_distance(vector, centroids[j], vector_size);
    if (min_distance > distance_from_centroid) {
      min_distance = distance_from_centroid;
      min_j = j;
    }
  }
  printf("find_closest_centroid_to_vector_index ending");

  return min_j;
}

double* create_new_centroid(double* centroid_sum, int centroid_counter, int centroid_size) {
  int p;
  double *new_centroid;
  printf("create_new_centroid starting");

  new_centroid = calloc(centroid_size, sizeof(double));
  if (new_centroid == NULL) {
  printf("create_new_centroid new_centroid == NULL error");

      return NULL;
  }

  for (p=0; p<centroid_size; p++) {
    new_centroid[p] = centroid_sum[p] / centroid_counter;
  }
  printf("create_new_centroid ending");

  return new_centroid;
}

int calculate_centroids_convergence(double** centroids, double** vectors, int centroids_num, int centroid_size, int vectors_num, int max_iterations, double eps) {
  int iter_couter, closest_centroid_index, i, j, p;
  double max_distance, centroids_distance, **centroids_sum, *counters, *new_centroid_j;
  max_distance = eps + 1;
  iter_couter = 0;
  printf("calculate_centroids_convergence starting");

  while (max_distance >= eps && iter_couter < max_iterations) {
    max_distance = 0;
    centroids_sum = (double**)calloc(centroids_num, sizeof(double *));
    if (centroids_sum == NULL) {
      printf("calculate_centroids_convergence centroids_sum == NULL failed");
      return 1;
    }
    counters = (double*)calloc(centroids_num, sizeof(double));
    if (counters == NULL) {
      printf("calculate_centroids_convergence counters == NULL failed");
      free(centroids_sum);
      return 1;
    }

    for (j=0; j<centroids_num; j++) {
      centroids_sum[j] = calloc(centroid_size, sizeof(double));

      if (centroids_sum[j] == NULL) {
        printf("calculate_centroids_convergence centroids_sum[j] == NULL failed");
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
          printf("calculate_centroids_convergence new_centroid_j == NULL failed");
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
    printf("calculate_centroids_convergence freeing memory in end of iteration for loop");
    free_array_of_pointers(centroids_sum, centroids_num);
    free(counters);

    iter_couter++;
  }

  return 0;
}

double** kmeans(int K, int N, int d, int iter, double eps, double** vectors, double** centroids) {
    int res;
    printf("kmeans starting");

    res = calculate_centroids_convergence(centroids, vectors, K, d, N, iter, eps);
    if (res == 1) {
      printf("kmeans error");
      free_array_of_pointers(vectors, N);
      free_array_of_pointers(centroids, K);
      return NULL;
    }
    printf("kmeans freeing pointers");
    free_array_of_pointers(vectors, N);
    printf("kmeans returning");

    return centroids;
}

static PyObject* convert_to_python_list(double **array, int rows, int cols) {
    PyObject *outer_list, *inner_list, *value;
    outer_list = PyList_New(rows);
    if (outer_list == NULL) {
      return NULL;
    }

    for (int i = 0; i < rows; i++) {
        inner_list = PyList_New(cols);
        if (inner_list == NULL) {
            Py_DECREF(outer_list);
            return NULL;
        }
        for (int j = 0; j < cols; j++) {
            value = PyFloat_FromDouble(array[i][j]);
            if (value == NULL) {
                Py_DECREF(inner_list);
                Py_DECREF(outer_list);
                return NULL;
            }
            PyList_SET_ITEM(inner_list, j, value);
        }
        PyList_SET_ITEM(outer_list, i, inner_list);
    }

    return outer_list;
}

static PyObject* fit(PyObject *self, PyObject *args) {
    int K, N, d, iter, i, j;
    double eps, **vectors, **centroids;
    PyObject *vectors_obj, *centroids_obj, *vector, *centroid, *new_centroids_obj;

    printf("fit starting");
    if(!PyArg_ParseTuple(args, "iiiidOO", &K, &N, &d, &iter, &eps, &vectors_obj, &centroids_obj)) {
        printf("fit error");
        return NULL;
    }
    printf("fit ending");

    printf("checking length");
    if (PyObject_Length(vectors_obj) < 0 || PyObject_Length(centroids_obj) < 0) {
      return NULL;
    }

    printf("memory vectors");
    vectors = (double **)calloc(N, sizeof(double *));
    if (vectors == NULL) {
      return NULL;
    }

    printf("memory centroids");
    centroids = (double **)calloc(K, sizeof(double *));
    if (centroids == NULL) {
      Py_DECREF(vectors_obj);
      Py_DECREF(centroids_obj);
      free(vectors);
      return NULL;
    }

    printf("Building vectors");
    for (i=0; i<N; i++) {
      vector = PyList_GetItem(vectors_obj, i);
      vectors[i] = calloc(d, sizeof(double));
      if (vectors[i] == NULL) {
        Py_DECREF(vector);
        Py_DECREF(vectors_obj);
        Py_DECREF(centroids_obj);
        free_array_of_pointers(vectors, i);
        free(centroids);
        return NULL;
      }
      for (j=0; j<d; j++) {
        vectors[i][j] = PyFloat_AsDouble(PyList_GetItem(vector, j));
      }
    }

    printf("Building centroids");
    for (i=0; i<K; i++) {
      centroid = PyList_GetItem(centroids_obj, i);
      centroids[i] = calloc(d, sizeof(double));
      if (vectors[i] == NULL) {
        Py_DECREF(vector);
        Py_DECREF(centroid);
        Py_DECREF(vectors_obj);
        Py_DECREF(centroids_obj);
        free_array_of_pointers(vectors, N);
        free_array_of_pointers(centroids, i);
        return NULL;
      }
      for (j=0; j<d; j++) {
        centroids[i][j] = PyFloat_AsDouble(PyList_GetItem(centroid, j));
      }
    }
    
    printf("Before kmeans");
    centroids = kmeans(K, N, d, iter, eps, vectors, centroids);
    printf("After kmeans");
    
    Py_DECREF(vector);
    Py_DECREF(centroid);
    Py_DECREF(vectors_obj);
    Py_DECREF(centroids_obj);

    printf("convert_to_python_list before");
    new_centroids_obj = convert_to_python_list(centroids, K, d);
    printf("convert_to_python_list after");

    if (new_centroids_obj == NULL) {
      return NULL;
    }
    printf("Returning final value");
    return Py_BuildValue("O", new_centroids_obj);
}

static PyMethodDef mykmeansspMethods[] = {
    {
      "fit",
      (PyCFunction) fit,
      METH_VARARGS,
      PyDoc_STR(
        "Perform some operation using K-means algorithm.\n \
        Input:\n \
        K : int - The number of centroids.\n \
        N : int - The number of data points.\n \
        d : int - The dimensionality of the data points.\n \
        iter : int - The maximum number of iterations for the algorithm.\n \
        eps : float - The threshold for convergence.\n \
        vectors : double** - A 2D list representing the input data points. Each inner list represents a data point and should have 'd' elements.\n \
        centroids : double** - A 2D list representing the initial centroids for the clusters. It should have 'K' inner lists, each representing a centroid point and should have 'd' elements.\n \
        "
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
    printf("PyInit_mykmeanssp starting");
    PyObject *m;
    m = PyModule_Create(&mykmeansspModule);
    if (!m) {
        printf("PyInit_mykmeanssp error");
        return NULL;
    }
    printf("PyInit_mykmeanssp ending");
    return m;
}

