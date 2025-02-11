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
void *writeTourToFile(int *tour, int tourLength, char *filename);

// Function to calculate the Euclidean distance between two points
double calculateDistance(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

void cheapestInsertion(double **coords, int numOfCoords, char *outputFileName) {
    int i, j, k;

    // Initialize a dynamic array to store the tour
    int *tour = (int *)malloc((numOfCoords + 1) * sizeof(int)); // +1 for the return to starting point

    // Initialize a dynamic array to keep track of visited vertices
    bool *visited = (bool *)malloc(numOfCoords * sizeof(bool));

    // Set the number of threads to 2
    omp_set_num_threads(2);

    // Initialize the tour with the first vertex
    tour[0] = 0;
    visited[0] = true;

    // Variables to store the best results
    double *minCostArray = (double *)malloc(numOfCoords * sizeof(double));
    int *bestVertexArray = (int *)malloc(numOfCoords * sizeof(int));
    int *bestPositionArray = (int *)malloc(numOfCoords * sizeof(int));

    // Record the start time
    double startTime = omp_get_wtime();

    // Find the closest unvisited vertex to the current partial tour
    for (i = 1; i < numOfCoords; i++) {
        // Initialize the best results for each iteration
        #pragma omp parallel for private(j, k)
        for (j = 0; j < numOfCoords; j++) {
            minCostArray[j] = DBL_MAX;
            bestVertexArray[j] = -1;
            bestPositionArray[j] = -1;
        }

        #pragma omp parallel for private(j, k)
        for (j = 0; j < numOfCoords; j++) {
            if (!visited[j]) {
                // Find the position to insert the vertex to minimize the cost
                for (k = 0; k < i; k++) {
                    double currentCost = calculateDistance(coords[tour[k % i]][0], coords[tour[k % i]][1],
                                                            coords[j][0], coords[j][1]) +
                                         calculateDistance(coords[tour[(k + 1) % i]][0], coords[tour[(k + 1) % i]][1],
                                                            coords[j][0], coords[j][1]) -
                                         calculateDistance(coords[tour[k % i]][0], coords[tour[k % i]][1],
                                                            coords[tour[(k + 1) % i]][0], coords[tour[(k + 1) % i]][1]);

                    if (currentCost < minCostArray[omp_get_thread_num()]) {
                        minCostArray[omp_get_thread_num()] = currentCost;
                        bestVertexArray[omp_get_thread_num()] = j;
                        bestPositionArray[omp_get_thread_num()] = k + 1;
                    }
                }
            }
        }

        // Combine the results from all threads to find the global minimum
        double minCost = DBL_MAX;
        int bestVertex = -1;
        int bestPosition = -1;
        for (j = 0; j < numOfCoords; j++) {
            if (minCostArray[j] < minCost) {
                minCost = minCostArray[j];
                bestVertex = bestVertexArray[j];
                bestPosition = bestPositionArray[j];
            }
        }

        // Insert the vertex at the best position in the tour
        for (k = i; k > bestPosition; k--) {
            tour[k] = tour[k - 1];
        }
        tour[bestPosition] = bestVertex;

        // Mark the selected vertex as visited
        visited[bestVertex] = true;
    }

    // Add the starting point to complete the tour
    tour[numOfCoords] = 0;

    // Record the end time
    double endTime = omp_get_wtime();

    // Print the tour
    printf("Cheapest Insertion Tour: ");
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
    free(minCostArray);
    free(bestVertexArray);
    free(bestPositionArray);
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

    // Call the cheapest insertion algorithm function
    cheapestInsertion(coords, numOfCoords, outputFileName);

    // Free dynamically allocated memory for coords
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
    }
    free(coords);

    return EXIT_SUCCESS;
}
