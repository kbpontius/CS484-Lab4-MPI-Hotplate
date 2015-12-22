//
//  main.c
//  Lab 4 - MPI Hotplate
//
//  Created by Kyle on 12/15/15.
//  Copyright © 2015 Kyle Pontius. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

// MARK: HOTPLATE DECLARARTIONS

#define MAX_ARRAY_SIZE 16384
#define EPSILON  0.1

double When()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

// MARK: HOTPLATE CODE

float getDifference(int i, int j, float** array) {
    float middle = array[i][j];
    
    float up = array[i - 1][j];
    float down = array[i + 1][j];
    float left = array[i][j - 1];
    float right = array[i][j + 1];
    
    float result = fabs(middle - ((down + up + right + left) / 4));
    
    return result;
}

float calculateNewCell(int i, int j, float** array) {
    float middle = array[i][j];
    
    float up = array[i - 1][j];
    float down = array[i + 1][j];
    float left = array[i][j - 1];
    float right = array[i][j + 1];
    
    return (down + up + right + left + (middle * 4)) / 8;
}

void writeCSV(float** arr) {
    FILE *fp = fopen("hotplateOutput.csv", "w+");
    int i, j;
    
    for (i = 0; i < MAX_ARRAY_SIZE; i++) {
        fprintf(fp, "%f", arr[i][0]);
        
        for (j = 1; j < MAX_ARRAY_SIZE; j++) {
            fprintf(fp, ", %f",arr[i][j]);
        }
        
        fprintf(fp, "\n");
    }
    
    fclose(fp);
}

void setupArrays(int theSize, int iproc, int nproc, float** newArray, float** oldArray) {
    int i, j;
    
    for (i = 0; i < theSize + 2; i++) {
        for (j = 0; j < MAX_ARRAY_SIZE; j++) {
            newArray[i][j] = 50;                // All other cells
            oldArray[i][j] = 50;
        }
        
        newArray[i][0] = 0;                     // Left side
        oldArray[i][0] = 0;
        
        newArray[i][MAX_ARRAY_SIZE - 1] = 0;    // Right side
        oldArray[i][MAX_ARRAY_SIZE - 1] = 0;
    }
    
    if (iproc == 0) {
        // Fill in the top row.
        int i;
        
        // 'newArray[1]...' here because the first proc has one extra row at the beginning.
        for (i = 0; i < MAX_ARRAY_SIZE; i++) {
            newArray[1][i] = 0;
            oldArray[1][i] = 0;
        }
    } else if (iproc == nproc - 1) {
        // Fill in the bottom row.
        int i;
        
        // 'theSize' here because the last proc has one extra row on the end.
        for (i = 0; i < MAX_ARRAY_SIZE; i++) {
            newArray[theSize][i] = 100;
            oldArray[theSize][i] = 100;
        }
    }
}

void allocateArray(float*** array, int theSize) {
    *array = (float **) malloc((theSize + 2) * (sizeof(float *)));
    int i;
    
    for (i = 0; i < theSize + 2; i++) {
        (*array)[i] = (float *) malloc((MAX_ARRAY_SIZE) * (sizeof(float)));
    }
}

void calculateNextState(int start, int end, float** newArray, float** oldArray) {
    int i, j;
    
    for (i = start; i < end; i++) {
        for (j = 1; j < MAX_ARRAY_SIZE - 1; j++) {
            newArray[i][j] = calculateNewCell(i, j, oldArray);
        }
    }
}

void swapArrays(float** newArray, float** oldArray) {
    float** tempArray = newArray;
    newArray = oldArray;
    oldArray = tempArray;
}

int isDone(int theSize, float** array) {
    int i, j;
    
//    fprintf(stderr, "SIZE: %i\n", theSize);
    
    for (i = 1; i < theSize + 1; i++) {
        for (j = 1; j < MAX_ARRAY_SIZE - 1; j++) {
            float result = getDifference(i, j, array);
            if (result > EPSILON) {
                fprintf(stderr, "%f > EPSILON.\n", result);
                return 0;
            }
        }
    }
    
    return 1;
}

// MARK: MAIN()

int main(int argc, char *argv[])
{
    int done, reallyDone;
    int count;
    int start, end;
    int theSize;
    char hostName[255];
    
    double startTime;
    
    int nproc, iproc;
    
    start = 1;
    
    MPI_Status status;
    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    
    // # of Rows
    theSize = MAX_ARRAY_SIZE / nproc;
    
    float** newArray;
    float** oldArray;
    
    gethostname(hostName, 253);
    startTime = When();
    
    fprintf(stderr,"%s: Hello from %d of %d. START: %i\n", hostName, iproc, nproc, theSize);
    
    // Allocate arrays.
    allocateArray(&newArray, theSize);
    allocateArray(&oldArray, theSize);
    
    // Setup arrays.
    setupArrays(theSize, iproc, nproc, newArray, oldArray);
    
    /* Now run the relaxation */
    reallyDone = 0;
    
    start = 1;
    end = theSize + 1;
    
    if (iproc == 0) {
        start = 2;
    } else if (iproc == nproc - 1) {
        end = theSize;
    }
    
    for(count = 0; !reallyDone; count++)
    {
        if (iproc == 0) {
            fprintf(stderr, "%i being calculateNextState() on count %i\n", iproc, count);
        }
        
        calculateNextState(start, end, newArray, oldArray);

        /* SEND ROWS DOWN */
        if (iproc > 0)
        {
            MPI_Send(&newArray[1], 1, MPI_FLOAT, iproc - 1, 0, MPI_COMM_WORLD);
        }
        
        // Last row shouldn't receive anything new.
        if (iproc < nproc - 1) {
            MPI_Recv(&newArray[theSize + 1], 1, MPI_FLOAT, iproc + 1, 0, MPI_COMM_WORLD, &status);
        }
        
        /* SEND ROWS UP */
        if (iproc < nproc - 1)
        {
            MPI_Send(&newArray[theSize], 1, MPI_FLOAT, iproc + 1, 0, MPI_COMM_WORLD);\
        }
        
        // First row shouldn't receive anything.
        if (iproc > 0) {
            MPI_Recv(&newArray[0], 1, MPI_FLOAT, iproc - 1, 0, MPI_COMM_WORLD, &status);
        }
        
//        fprintf(stderr, "%d: Checking if isDone.\n", iproc);
        
        // Calculate new done value
        done = isDone(theSize, newArray);
        
        // Reduce all done values.
        MPI_Allreduce(&done, &reallyDone, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        
        // TODO: Get code example for MPI_Allreduce.
        if (reallyDone) { break; }
        
        swapArrays(newArray, oldArray);
    }
    
    /* print out the number of iterations to relax */
    fprintf(stderr, "%d:It took %d iterations and %lf seconds to relax the system\n", iproc, count, When() - startTime);
    MPI_Finalize();
}
