#include <stdio.h>        // Standard I/O functions
#include <stdlib.h>       // Standard library functions, including dynamic memory allocation
#include <math.h>         // Mathematical functions, particularly for sqrt()
#include <float.h>        // Defines constants for floating point types like DBL_MAX
#include <stdbool.h>      // Boolean data type (true/false)
#include <omp.h>          // OpenMP library for parallel processing

// Function prototypes
int readNumOfCoords(char *fileName);                       // Reads the number of coordinates from the input file
double **readCoords(char *filename, int numOfCoords);      // Reads coordinates from the input file
void writeTourToFile(int *tour, int tourLength, char *filename);  // Writes the calculated tour to an output file

// Calculates the Euclidean distance between two points (x1, y1) and (x2, y2)
double calculateDistance(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)); // Using the distance formula
}

// Generates a matrix where distMatrix[i][j] contains the distance between coordinate i and coordinate j
void calculateDistanceMatrix(double **coords, double **distMatrix, int numOfCoords) {
    for (int i = 0; i < numOfCoords; i++) { // Loop over each coordinate i
        for (int j = 0; j < numOfCoords; j++) { // Loop over each coordinate j
            distMatrix[i][j] = calculateDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]); // Store distance in matrix
        }
    }
}

// Implements the Smallest Sum Insertion algorithm to find a tour
void smallestSumInsertion(double **distMatrix, int numOfCoords, int *tour) {
    bool *visited = (bool *)malloc(numOfCoords * sizeof(bool)); // Array to track visited nodes
    for (int i = 0; i < numOfCoords; i++) visited[i] = false;  // Initialize all nodes as unvisited

    tour[0] = 0;          // Start the tour from the first node (index 0)
    visited[0] = true;    // Mark the first node as visited
    int tourLength = 1;   // Current length of the tour

    while (tourLength < numOfCoords) { // Continue until all nodes are visited
        int nextVertex = -1; // Variable to store the next vertex to add
        double minSum = DBL_MAX; // Initialize minimum sum to maximum double value

        #pragma omp parallel // Start a parallel region
        {
            int localNextVertex = -1; // Thread-local variable for the next vertex
            double localMinSum = DBL_MAX; // Thread-local variable for the minimum sum

            #pragma omp for // Parallelize the for loop
            for (int i = 0; i < numOfCoords; i++) { // Loop over all nodes
                if (!visited[i]) { // Consider only unvisited nodes
                    double sum = 0.0; // Sum of distances to current tour
                    for (int j = 0; j < tourLength; j++) {
                        sum += distMatrix[tour[j]][i]; // Calculate sum of distances
                    }
                    if (sum < localMinSum) { // Check if this sum is the minimum found
                        localMinSum = sum;
                        localNextVertex = i; // Update next vertex for this thread
                    }
                }
            }

            #pragma omp critical // Ensure only one thread updates the global minimum at a time
            {
                if (localMinSum < minSum) { // Check if this thread found a new global minimum
                    minSum = localMinSum;
                    nextVertex = localNextVertex; // Update global next vertex
                }
            }
        }

        // Find the best position to insert the next vertex in the tour
        int bestPos = -1; // Position to insert next vertex
        double minCost = DBL_MAX; // Initialize minimum insertion cost to maximum double value
        for (int i = 0; i < tourLength; i++) {
            int j = (i + 1) % tourLength; // Wrap around for circular tour
            double cost = distMatrix[tour[i]][nextVertex] + distMatrix[tour[j]][nextVertex] - distMatrix[tour[i]][tour[j]]; // Calculate insertion cost
            if (cost < minCost) { // Check if this cost is the minimum found
                minCost = cost;
                bestPos = i; // Update best position
            }
        }

        // Insert the next vertex in the tour
        for (int i = tourLength; i > bestPos + 1; i--) {
            tour[i] = tour[i - 1]; // Shift elements to make space for the new vertex
        }
        tour[bestPos + 1] = nextVertex; // Insert new vertex
        visited[nextVertex] = true; // Mark the vertex as visited
        tourLength++; // Increase tour length
    }
    tour[tourLength] = 0; // Closing the tour by returning to the starting point
    free(visited); // Free memory allocated for visited array
}

// Writes the tour to a file and also prints it to the console
void writeTourToFileAndConsole(int *tour, int tourLength, char *filename) {
    FILE *file = fopen(filename, "w"); // Open the file for writing
    if (file == NULL) { // Check if the file was opened successfully
        printf("Unable to open file: %s\n", filename);
        return;
    }

    fprintf(file, "%d\n", tourLength); // Write the length of the tour to the file
    printf("Tour: "); // Print the tour to the console
    for (int i = 0; i < tourLength; i++) { // Loop through the tour
        fprintf(file, "%d ", tour[i]); // Write each vertex of the tour to the file
        printf("%d ", tour[i]); // Print each vertex to the console
    }
    printf("\n");

    fclose(file); // Close the file
}

// Main function: entry point of the program
int main(int argc, char *argv[]) {
    if (argc != 3) { // Check if the program has received exactly 2 arguments
        printf("Usage: %s <input file> <output file>\n", argv[0]);
        return 1; // Exit with an error code
    }

    char *inputFile = argv[1]; // Get input file name from arguments
    char *outputFile = argv[2]; // Get output file name from arguments

    // Record the start time
    double startTime = omp_get_wtime(); // Get the current wall clock time

    int numOfCoords = readNumOfCoords(inputFile); // Read the number of coordinates from the input file
    if (numOfCoords < 0) { // Check if reading the number of coordinates failed
        printf("Error reading number of coordinates from file.\n");
        return 1; // Exit with an error code
    }

    double **coords = readCoords(inputFile, numOfCoords); // Read the coordinates from the input file
    if (coords == NULL) { // Check if reading coordinates failed
        printf("Error reading coordinates from file.\n");
        return 1; // Exit with an error code
    }

    double **distMatrix = (double **)malloc(numOfCoords * sizeof(double *)); // Allocate memory for the distance matrix
    for (int i = 0; i < numOfCoords; i++) {
        distMatrix[i] = (double *)malloc(numOfCoords * sizeof(double)); // Allocate memory for each row of the distance matrix
    }

    calculateDistanceMatrix(coords, distMatrix, numOfCoords); // Calculate the distance matrix using the coordinates

    int *tour = (int *)malloc((numOfCoords + 1) * sizeof(int)); // Allocate memory for the tour array (+1 for return to origin)
    smallestSumInsertion(distMatrix, numOfCoords, tour); // Generate the tour using the Smallest Sum Insertion algorithm

    // Record the end time before writing to the output file
    double endTime = omp_get_wtime(); // Get the current wall clock time

    writeTourToFileAndConsole(tour, numOfCoords + 1, outputFile); // Write the tour to the file and console (+1 for return to origin)

    // Print the execution time
    printf("Execution Time: %f seconds\n", endTime - startTime); // Calculate and print the total execution time

    // Free allocated memory for coordinates, distance matrix, and tour
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
        free(distMatrix[i]);
    }
    free(coords);
    free(distMatrix);
    free(tour);

    return 0; // Exit the program successfully
}
