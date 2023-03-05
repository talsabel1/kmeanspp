#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* function declarations */
double euclidean_distance(double *vector1, double *vector2, int max_d);

void kmeans(int num_of_clusters, int num_of_iterations, int vector_dimension, int count, double vector_list[][vector_dimension], double eps, double init_centroids[][vector_dimension]);

double euclidean_distance(double *vector1, double *vector2, int max_d) {
    double total_sum = 0;
    int i;
    for (i = 0; i < max_d; i++) {
        total_sum += pow(vector1[i] - vector2[i], 2);
    }
    return sqrt(total_sum);

}

void kmeans(int num_of_clusters, int num_of_iterations, int vector_dimension, int count, double vector_list[][vector_dimension], double eps, double init_centroids[][vector_dimension]) {

    double **centroids = calloc(num_of_clusters, sizeof(double *));
    int *cluster_sizes_copy = calloc(num_of_clusters, sizeof(int));
    double *temp_centroid = calloc(vector_dimension, sizeof(double));
    int *cluster_sizes = calloc(num_of_clusters, sizeof(int));
    double ***clusters = calloc(num_of_clusters, sizeof(double **));
    double max_distance = -1;
    int min_index1;
    int min_index2;
    int i;
    int j;
    int m;
    int n;
    int l;
    int k;
    int p;
    int q;
    double min_distance1;
    double min_distance2;
    double distance1;
    double distance2;
    double sum;


    for (i = 0; i < num_of_clusters; i++) {
        centroids[i] = calloc(vector_dimension, sizeof(double));
    }
    for (i = 0; i < num_of_clusters; i++) {
        for (j = 0; j < vector_dimension; j++) {
            centroids[i][j] = init_centroids[i][j];
        }
    }
    for (i = 0; i < num_of_iterations; i++) {

        for (j = 0; j < count; j++) {
            min_distance1 = euclidean_distance(vector_list[j], centroids[0], vector_dimension);
            min_index1 = 0;
            for (k = 1; k < num_of_clusters; k++) {
                distance1 = euclidean_distance(vector_list[j], centroids[k], vector_dimension);
                if (distance1 < min_distance1) {
                    min_distance1 = distance1;
                    min_index1 = k;
                }
            }
            cluster_sizes[min_index1] += 1;
        }
        for (j = 0; j < num_of_clusters; j++) {
            cluster_sizes_copy[j] = cluster_sizes[j];
        }
        for (j = 0; j < num_of_clusters; j++) {
            clusters[j] = calloc(cluster_sizes[j], sizeof(double *));
            for (m = 0; m < cluster_sizes[j]; m++) {
                clusters[j][m] = calloc(vector_dimension, sizeof(double));
            }
        }
        for (j = 0; j < count; j++) {
            min_distance2 = euclidean_distance(vector_list[j], centroids[0], vector_dimension);
            min_index2 = 0;
            for (p = 1; p < num_of_clusters; p++) {
                distance2 = euclidean_distance(vector_list[j], centroids[p], vector_dimension);
                if (distance2 < min_distance2) {
                    min_distance2 = distance2;
                    min_index2 = p;
                }
            }
            for (q = 0; q < vector_dimension; q++) {
                clusters[min_index2][cluster_sizes[min_index2] - 1][q] = vector_list[j][q];
            }
            cluster_sizes[min_index2] -= 1;
        }
        max_distance = -1;
        for (j = 0; j < num_of_clusters; j++) {
            for (l = 0; l < vector_dimension; l++) {
                temp_centroid[l] = centroids[j][l];
            }
            for (m = 0; m < vector_dimension; m++) {
                sum = 0;
                for (n = 0; n < cluster_sizes_copy[j]; n++) {
                    sum += clusters[j][n][m];
                }
                centroids[j][m] = sum / cluster_sizes_copy[j];

            }
            if (euclidean_distance(temp_centroid, centroids[j], vector_dimension) > max_distance) {
                max_distance = euclidean_distance(temp_centroid, centroids[j], vector_dimension);
            }
        }
        if (max_distance <= eps) {
            break;
        }
        for (i = 0; i < num_of_clusters; i++) {
            for (j = 0; j < cluster_sizes_copy[i]; j++) {
                free(clusters[i][j]);
            }
            free(clusters[i]);
        }
    }
    for (i = 0; i < num_of_clusters; i++) {
        for (j = 0; j < vector_dimension; j++) {
            if (j == vector_dimension - 1)
                printf("%.4f\n", centroids[i][j]);
            else
                printf("%.4f%c", centroids[i][j], ',');
        }
    }

    /* free memory */
    for (i = 0; i < num_of_clusters; i++) {
        free(centroids[i]);
    }
    free(centroids);
    free(temp_centroid);
    free(cluster_sizes);
    for (i = 0; i < num_of_clusters; i++) {
        for (j = 0; j < cluster_sizes_copy[i]; j++) {
            free(clusters[i][j]);
        }
        free(clusters[i]);
    }
    free(clusters);
    free(cluster_sizes_copy);

}

static PyObject* fit(PyObject *self, PyObject *args) {
    int num_of_clusters;
    int num_of_iterations;
    int vector_dimension;
    int count;
    PyObject *pyObjVector_list;
    double eps;
    PyObject *pyObjInit_centroids;
    PyObject *item_row;
    PyObject *item_col;
    double num;
    if (!PyArg_ParseTuple(args, "iiiiOdO", &num_of_clusters, &num_of_iterations, &vector_dimension, &count, &pyObjVector_list, &eps, &pyObjInit_centroids)) {
        return NULL;
    }

    if (vector_dimension < 0 || count < 0){
        return NULL;
    }

    double vector_list[count][vector_dimension];
    int i;
    int j;
    for (i = 0; i < count; i++) {
        item_row = PyList_GetItem(pyObjVector_list, i);
        for (j = 0; j < vector_dimension; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            vector_list[i][j] = num;
        }
    }

    double init_centroids[count][vector_dimension];
    for (i = 0; i < num_of_clusters; i++) {
        item_row = PyList_GetItem(pyObjInit_centroids, i);
        for (j = 0; j < vector_dimension; j++) {
            item_col = PyList_GetItem(item_row, j);
            num = PyFloat_AsDouble(item_col);
            init_centroids[i][j] = num;
        }
    }

    kmeans(num_of_clusters, num_of_iterations, vector_dimension, count, vector_list, eps, init_centroids);

    Py_RETURN_NONE;
}

static PyMethodDef kmeansMethods[] = {
        {
                "fit", // name exposed to Python
                fit, // C wrapper function
                     METH_VARARGS, // received variable args (but really just 1)
                "Performs the kmeans algorithm according to the given centroid initialization" // documentation
        },
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mykmeanssp = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",     // name of module exposed to Python
        "A module that performs the kmeans algorithm", // module documentation
        -1,
        kmeansMethods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&mykmeanssp);
    if (!m) {
        return NULL;
    }
    return m;
}