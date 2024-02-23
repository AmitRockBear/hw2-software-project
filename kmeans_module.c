
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define eps 0.001

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

double** stdin_to_matrix(int rows, int columns) {
  double **vectors;
  int line_count, d_counter;
  char *line, *start_iterator, *end_iterator;
  size_t len;
  
  vectors = (double**)calloc(rows, sizeof(double *));
  if (vectors == NULL) {
      return NULL;
  }

  len = 0;
  line = NULL;
  line_count = 0;
  while (getline(&line, &len, stdin) != -1 && line_count < rows) {
      vectors[line_count] = calloc(columns, sizeof(double));
      if (vectors[line_count] == NULL) {
        free_array_of_pointers(vectors, line_count);
        free(line);
        return NULL;
      }

      start_iterator = line;
      end_iterator = line;
      d_counter = 0;
      while(*end_iterator != '\n' && *end_iterator != '\0') {
        while(*end_iterator != ',' && *end_iterator != '\n' && *end_iterator != '\0') {
          end_iterator++;
        }
        vectors[line_count][d_counter] = strtod(start_iterator, NULL);
        if (*end_iterator != '\n' && *end_iterator != '\0') {
          end_iterator++;
          start_iterator = end_iterator;
        }
        d_counter++;
      }

      line_count++;
  }

  free(line);

  return vectors;
}

double* deep_copy_vector(double* vector, int dimension) {
  double* copied_vector;

  copied_vector = calloc(dimension, sizeof(double));
  if (copied_vector == NULL) {
    return NULL;
  }

  for (j=0; j<dimension; j++) {
    copied_vector[j] = vector[j];
  }

  return copied_vector;
}

double** init_centroids(double** vectors, int vectors_num, int centroids_num, int dimension) {
  double distance, min_distance, sum;
  int k;
  int i, j, chosen, random_number;
  double **centroids, *probabilities;
  
  srand(time(NULL));
  
  chosen = rand() % vectors_num;

  centroids = calloc((size_t)centroids_num, sizeof(double *));
  if (centroids == NULL) {
      return NULL;
  }

  centroids[0] = deep_copy_vector(vectors[chosen], dimension);
  if (centroids[0] == NULL) {
    free(centroids);
    return NULL;
  }

  for (i=1; i<centroids_num; i++) {
    probabilities = calloc((size_t)vectors, sizeof(double *));
    if (probabilities == NULL) {
        free_array_of_pointers(centroids, i);
        return NULL;
    }

    // Calculate min_distance
    for (j=0; j<vectors_num; j++) {
      min_distance = distance_from_centroid(vectors[j], centroids[0], dimension);
      for (k=1; k<i; k++) {
        distance =  distance_from_centroid(vectors[j], centroids[k], dimension);
        if (distance < min_distance) {
          min_distance = distance;
        }
      }
      probabilities[j] = min_distance;
    }

    sum = 0;
    for (j=0; j<vectors_num; j++) {
      sum = sum + probabilities[j];
    }

    for (j=0; j<vectors_num; j++) {
      probabilities[j] = probabilities[j] / sum;
    }

    sum = 0;
    for (j=0; j<vectors_num; j++) {
      sum = sum + probabilities[j]
      probabilities[j] = sum;
    }

    random_number = rand() % 1;
    for (j=0; j<vectors_num; j++) {
      if (random_number < probabilities[j]) {
        chosen = j;
        break;
      }
    }

    centroids[i] = deep_copy_vector(vectors[chosen], dimension);
    if (centroids[i] == NULL) {
      free_array_of_pointers(centroids, i);
      return NULL;
    }
  }

  return centroids;
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

int calculate_centroids_convergence(double** centroids, double** vectors, int centroids_num, int centroid_size, int vectors_num, int max_iterations) {
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

void print_output(double** centroids, int centroids_num, int centroid_size) {
  int i, j;

  for (i=0; i<centroids_num; i++) {
    for (j=0; j<centroid_size; j++) {
      printf("%.4f",centroids[i][j]);
      if (j != centroid_size-1)
        printf(",");
      else
        printf("\n");
    }
  }
}

int isInteger(const char *str) {
    if (*str == '\0')
        return 0;

    while (*str != '\0') {
        if (!('0' <= *str && *str <= '9'))
            return 0;
        str++;
    }
    return 1;
}

int main(int argc, char* argv[]) {
    int K, N, d, iter, res;
    double **vectors, **centroids;

    if (argc != 4 && argc != 5) {
      printf("An Error Has Occurred");
      return 1;
    }

    iter = 200;

    if (isInteger(argv[2]) == 0) {
      printf("Invalid number of points!");
      return 1;
    }
    N = atoi(argv[2]);
    if (N <= 1) {
      printf("Invalid number of points!");
      return 1;
    }

    if (isInteger(argv[1]) == 0) {
      printf("Invalid number of clusters!");
      return 1;
    }
    K = atoi(argv[1]);
    if (K <= 1 || K >= N) {
      printf("Invalid number of clusters!");
      return 1;
    }

    if (isInteger(argv[3]) == 0) {
      printf("Invalid dimension of point!");
      return 1;
    }
    d = atoi(argv[3]);
    if (d < 1) {
      printf("Invalid dimension of point!");
      return 1;
    }

    if (argc == 5) {
      if (isInteger(argv[4]) == 0) {
        printf("Invalid maximum iteration!");
        return 1;
      }
      iter = atoi(argv[4]);
      if (iter <= 1 || iter >= 1000) {
        printf("Invalid maximum iteration!");
        return 1;
      }
    }

    vectors = stdin_to_matrix(N, d);
    if (vectors == NULL) {
      printf("An Error Has Occurred");
      return 1;
    }

    centroids = deep_copy_matrix(vectors, N, K, d);
    if (centroids == NULL) {
      free_array_of_pointers(vectors, N);
      printf("An Error Has Occurred");
      return 1;
    }

    res = calculate_centroids_convergence(centroids, vectors, K, d, N, iter);
    if (res == 1) {
      free_array_of_pointers(vectors, N);
      free_array_of_pointers(centroids, K);
      printf("An Error Has Occurred");
      return 1;
    }

    print_output(centroids, K, d);

    free_array_of_pointers(vectors, N);
    free_array_of_pointers(centroids, K);

    return 0;
}

