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

// MARK: MPI DECLARARTIONS

#define VECSIZE 8
#define ITERATIONS 10000

typedef struct {
    double val;
    int rank;
} element;

// MARK: HOTPLATE DECLARARTIONS

#define MAX_ARRAY_SIZE 16384
#define EPSILON  0.1

float** newArray;
float** oldArray;
float** tempArray;

double When()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

// MARK: HOTPLATE CODE

void swapArrays() {
    /* Swap the pointers */
    tempArray = newArray;
    newArray = oldArray;
    oldArray = tempArray;
}

float getDifference(int i, int j) {
    float middle = newArray[i][j];
    
    float up = newArray[i - 1][j];
    float down = newArray[i + 1][j];
    float left = newArray[i][j - 1];
    float right = newArray[i][j + 1];
    
    return fabs(middle - ((down + up + right + left) / 4));
}

float calculateNewCell(int i, int j) {
    float middle = oldArray[i][j];
    
    float up = oldArray[i - 1][j];
    float down = oldArray[i + 1][j];
    float left = oldArray[i][j - 1];
    float right = oldArray[i][j + 1];
    
    return (down + up + right + left + (middle * 4)) / 8;
}

void writeCSV(float** arr) {
    FILE *fp = fopen("hotplateOutput.csv", "w+");
    
    for (int i = 0; i < MAX_ARRAY_SIZE; i++) {
        fprintf(fp, "%f", arr[i][0]);
        
        for (int j = 1; j < MAX_ARRAY_SIZE; j++) {
            fprintf(fp, ", %f",arr[i][j]);
        }
        
        fprintf(fp, "\n");
    }
    
    fclose(fp);
}

void setupArrays() {
    for (int i = 0; i < MAX_ARRAY_SIZE; i++) {
        for (int j = 0; j < MAX_ARRAY_SIZE; j++) {
            newArray[i][j] = 50;                // All other cells
            oldArray[i][j] = 50;
        }
        
        newArray[0][i] = 0;                     // Top row
        oldArray[0][i] = 0;
        
        newArray[i][0] = 0;                     // Left side
        oldArray[i][0] = 0;
        
        newArray[i][MAX_ARRAY_SIZE - 1] = 0;    // Right side
        oldArray[i][MAX_ARRAY_SIZE - 1] = 0;
        
        newArray[MAX_ARRAY_SIZE - 1][i] = 100;  // Bottom row
        oldArray[MAX_ARRAY_SIZE - 1][i] = 100;
    }
}

void allocateArray(float*** array, int arraySize) {
    *array = (float **) malloc((arraySize + 2) * (sizeof(float *)));
    
    for (int i = 0; i < arraySize + 2; i++) {
        (*array)[i] = (float *) malloc((arraySize + 2) * (sizeof(float)));
    }
}

int arrayIsFinished() {
    float difference = 0;
    float isFinished = 1;
    
    // Avoid iterating over borders.
    for (int i = 1; i < MAX_ARRAY_SIZE - 1; i++) {
        if (isFinished == 1) {
            for (int j = 1; j < MAX_ARRAY_SIZE - 1; j++) {
                difference = getDifference(i, j);
                
                if (difference >= 0.1) {
                    //                        printf("\nDifference: %f", difference);
                    //                        printf("\nRow: %d || Col: %d", i, j);
                    isFinished = 0;
                    break;
                }
            }
        }
    }
    
    return isFinished;
}

int calculateSteadyState() {
    int isFinished = 0;
    int iterations = 0;
    
    while (isFinished == 0) {
        iterations++;
        swapArrays();
        
        // Avoid iterating over borders.
        for (int i = 1; i < MAX_ARRAY_SIZE - 1; i++) {
            for (int j = 1; j < MAX_ARRAY_SIZE - 1; j++) {
                newArray[i][j] = calculateNewCell(i, j);
            }
        }
        
        isFinished = arrayIsFinished();
    }
    
    return iterations;
}

// MARK: MAIN()

int main(int argc, char *argv[])
{
    int i, done, reallydone;
    int cnt;
    int start, end;
    int theSize;
    
    double startTime;
    
    int nproc, iproc;
    MPI_Status status;
    
    MPI_Init(&argc, &argv);
    startTime = When();
    
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    fprintf(stderr,"%d: Hello from %d of %d\n", iproc, iproc, nproc);
    
    /* Determine how much I should be doing and allocate the arrays */
    theSize = MAX_ARRAY_SIZE / nproc;
    allocateArray(&newArray, theSize);
    allocateArray(&oldArray, theSize);
    
    start = 1;
    end = theSize + 1;
    
    /* Initialize the cells */
    for (i = 0; i < theSize + 2; i++)
    {
        newArray[i] = oldArray[i] = 50;
    }
    
    /* Initialize the Boundaries */
    if (iproc == 0)
    {
        start = 2;
        newArray[1] = oldArray[1] = 100;
    }
    if (iproc == nproc - 1)
    {
        end = theSize;
        newArray[theSize] = oldArray[theSize] = 0;
    }
    
    /* Now run the relaxation */
    reallydone = 0;
    for(cnt = 0; !reallydone; cnt++)
    {
        /* First, I must get my neighbors boundary values */
        if (iproc != 0)
        {
            MPI_Send(&oldArray[1], 1, MPI_FLOAT, iproc - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&oldArray[0], 1, MPI_FLOAT, iproc - 1, 0, MPI_COMM_WORLD, &status);
        }
        
        if (iproc != nproc - 1)
        {
            MPI_Send(&oldArray[theSize], 1, MPI_FLOAT, iproc + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&oldArray[theSize + 1], 1, MPI_FLOAT, iproc + 1, 0, MPI_COMM_WORLD, &status);
        }
        
        /* Do the calculations */
        for (i = start; i < end; i++)
        {
            newArray[i] = (oldArray[i-1] + oldArray[i+1] + 2 * oldArray[i]) / 4.0;
        }
        
        /* Check to see if we are done */
        done = 1;
        for (i = start; i < end; i++)
        {
            if (fabs((newArray[i-1] + newArray[i+1])/ 2.0 - newArray[i]) > EPSILON)
            {
                done = 0;
                break;
            }
        }
        
        /* Do a reduce to see if everybody is done */
        MPI_Allreduce(&done, &reallydone, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        
        swapArrays();
    }
    
    /* print out the number of iterations to relax */
    fprintf(stderr, "%d:It took %d iterations and %lf seconds to relax the system\n", iproc, cnt, When() - startTime);
    MPI_Finalize();
}
