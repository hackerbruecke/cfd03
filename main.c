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
#include <stdlib.h>

int main(int argc, char *argv[]) {
	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;
	int xlength[3];
	double tau;
	double velocityWall[D];
	int timesteps;
	int timestepsPerPlotting;

    printf("LBM simulation by Krivokapic, Mody, Malcher - CFD Lab SS2014\n");
    printf("============================================================\n");
    printf("Reading in parameters...\n");
	if (readParameters(xlength, &tau, velocityWall, &timesteps,
			&timestepsPerPlotting, argc, argv)) {
		printf("Reading in parameters failed. Aborting program!\n");
		exit(-1);
	}
    printf("...done\n");
    
    printf("Starting LBM...\n");
    /* (xlength+2)^D elements must be stored for all lattices including boundaries */
    const int xl_to3 = (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
    collideField = malloc(sizeof *collideField * Q * xl_to3);
    streamField = malloc(sizeof *streamField * Q * xl_to3);
    flagField = malloc(sizeof *flagField * xl_to3);
    
    /* Initialize pointers */
	initialiseFields(collideField, streamField, flagField, xlength);

    for (int t = 0; t < timesteps; ++t) {
		double *swap = NULL;
		doStreaming(collideField, streamField, flagField, xlength);
		swap = collideField;
		collideField = streamField;
		streamField = swap;

		doCollision(collideField, flagField, &tau, xlength);
		treatBoundary(collideField, flagField, velocityWall, xlength);

        /* Write output to vtk file for postprocessing */
		if (t % timestepsPerPlotting == 0) {
			writeVtkOutput(collideField, flagField, "lbm_out", t, xlength);
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

