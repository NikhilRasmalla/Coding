#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <omp.h>

// Include coordReader.c functions
int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void writeTourToFile(int *tour, int tourLength, char *filename);

// Function to calculate the Euclidean distance between two points
double calculateDistance(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

void farthestInsertion(double **coords, int numOfCoords, char *outputFileName) {
    int i, j, k;

    // Set the number of threads
    omp_set_num_threads(2);

    // Initialize a dynamic array to store the tour
    int *tour = (int *)malloc((numOfCoords + 1) * sizeof(int)); // +1 for the return to starting point

    // Initialize a dynamic array to keep track of visited vertices
    bool *visited = (bool *)malloc(numOfCoords * sizeof(bool));

    // Initialize the tour with the first vertex
    tour[0] = 0;
    visited[0] = true;


    // Declare distances array
    double **distances = (double **)malloc(numOfCoords * sizeof(double *));
    for (i = 0; i < numOfCoords; i++) {
        distances[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    // Record the start time
    double startTime = omp_get_wtime();

    // Calculate distances in parallel
    #pragma omp parallel for private(j, k)
    for (i = 0; i < numOfCoords; i++) {
        for (j = 0; j < numOfCoords; j++) {
            distances[i][j] = calculateDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
        }
    }

    // Find the farthest vertex from the starting point and create the initial partial tour
    int farthestVertex = -1;
    double maxDistance = -DBL_MAX;
    for (j = 1; j < numOfCoords; j++) {
        double distance = distances[0][j];
        if (distance > maxDistance) {
            maxDistance = distance;
            farthestVertex = j;
        }
    }
    tour[1] = farthestVertex;
    tour[2] = 0;
    visited[farthestVertex] = true;

    // Build the tour by inserting farthest unvisited vertices
    for (i = 3; i <= numOfCoords; i++) {
        farthestVertex = -1;
        maxDistance = -DBL_MAX;

        // Find the farthest unvisited vertex from any vertex in the partial tour
        #pragma omp parallel for private(j, k) reduction(max:maxDistance)
        for (j = 0; j < i - 1; j++) {
            for (k = 1; k < numOfCoords; k++) {
                if (!visited[k]) {
                    double distance = distances[tour[j]][k];
                    if (distance > maxDistance) {
                        #pragma omp critical
                        {
                            if (distance > maxDistance) {
                                maxDistance = distance;
                                farthestVertex = k;
                            }
                        }
                    }
                }
            }
        }

        // Find the position to insert the farthest vertex with minimal cost
        int position = -1;
        double minCost = DBL_MAX;
        for (j = 0; j < i - 1; j++) {
            double cost = distances[tour[j]][farthestVertex]
                          + distances[farthestVertex][tour[j + 1]]
                          - distances[tour[j]][tour[j + 1]];

            if (cost < minCost) {
                minCost = cost;
                position = j + 1;
            }
        }

        // Insert the farthest vertex at the optimal position
        for (j = i; j > position; j--) {
            tour[j] = tour[j - 1];
        }
        tour[position] = farthestVertex;
        visited[farthestVertex] = true;
    }

    // Record the end time
    double endTime = omp_get_wtime();

    // Print the tour
    printf("Farthest Insertion Tour: ");
    for (k = 0; k <= numOfCoords; k++) {
        printf("%d ", tour[k]);
    }
    printf("\n");

    // Write the tour to the specified output file
    writeTourToFile(tour, numOfCoords + 1, outputFileName); // +1 for the return to starting point

    // Print the execution time
    printf("Execution Time: %f seconds\n", endTime - startTime);

    // Clean up
    free(tour);
    free(visited);

    // Clean up distances array
    for (i = 0; i < numOfCoords; i++) {
        free(distances[i]);
    }
    free(distances);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    char *inputFileName = argv[1];
    char *outputFileName = argv[2];

    int numOfCoords = readNumOfCoords(inputFileName);
    double **coords = readCoords(inputFileName, numOfCoords);

    // Call the farthest insertion algorithm function
    farthestInsertion(coords, numOfCoords, outputFileName);

    // Free dynamically allocated memory for coords
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
    }
    free(coords);

    return EXIT_SUCCESS;
}
