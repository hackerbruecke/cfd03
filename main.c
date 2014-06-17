#ifndef _MAIN_C_
#define _MAIN_C_

#include "LBDefinitions.h"
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "initLB.h"

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

	/*Declaring Variables required for LBM*/
	char *problem = NULL;
	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;
	int xlength[3];
	double tau;
	double velocityWall[D];
	int timesteps;
	int timestepsPerPlotting;
	double denref = 1.0;
	double dense = 1.0;

	problem = argv[2];

    printf("LBM simulation by Krivokapic, Malcher, Mody - CFD Lab SS2014\n");
    printf("============================================================\n");
    printf("Reading in parameters...\n");
	if (readParameters(xlength, &tau, velocityWall, &timesteps,
			&timestepsPerPlotting, problem, denref, &dense, argc, argv)) {
		printf("Reading in parameters failed. Aborting program!\n");
		exit(-1);
	}
    printf("...done\n");
    
    printf("Starting LBM...\n");
    /* (xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2) elements including boundaries */
    const int xl_to3 = (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);

    /*Allocate memory for Pointers*/
    collideField = malloc(sizeof *collideField * Q * xl_to3);
    streamField = malloc(sizeof *streamField * Q * xl_to3);
    flagField = malloc(sizeof *flagField * xl_to3);

    /* Initialize pointers */
	initialiseFields(collideField, streamField, flagField, xlength, problem);

	for (int t = 0; t < timesteps; ++t) {
		/*Perform the LBM Step*/
		double *swap = NULL;
		doStreaming(collideField, streamField, flagField, xlength);

		/*Exchange Streaming and Collision Fields*/
		swap = collideField;
		collideField = streamField;
		streamField = swap;

		/*Perform Collision Step*/
		doCollision(collideField, flagField, &tau, xlength);
		treatBoundary(collideField, flagField, velocityWall, xlength, problem, denref, &dense);

        /* Write output to vtk file for Post-processing */
		if (t % timestepsPerPlotting == 0) {
			writeVtkOutput(collideField, flagField, "lbm_", t, xlength);
		}
	}
    printf("============================================================\n");
    printf("LBM simulation completed for %d cells and %d timesteps\n", xl_to3, timesteps);
    printf("Freeing allocated memory...\n");

    free(collideField);
    free(streamField);
    free(flagField);
	return 0;
}

#endif

