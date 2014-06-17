#include "collision.h"
#include "LBDefinitions.h"
#include <stdio.h>
#include <math.h>

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){

    /* Updating all(19) post collision distribution fields using BGK approx rule in Current Cell */
    for (int i = 0; i < Q; ++i) {
        currentCell[i] = currentCell[i] - (currentCell[i] - feq[i]) / *tau;
    }
}



void doCollision(double *collideField, int *flagField,const double * const tau, int *xlength){
    double density;
    double velocity[D];
    double feq[Q];
    for (int z = 0; z < xlength[2] + 2; ++z) {
        for (int y = 0; y < xlength[1] + 2; ++y) {
            for (int x = 0; x < xlength[0] + 2; ++x) {
                double *currentCell = &collideField[idx(xlength, x, y, z, 0)];
	
                /*Updating values for velocity, density and Feq for Current cell*/                
                computeDensity(currentCell, &density);
                computeVelocity(currentCell, &density, velocity);
                computeFeq(&density, velocity, feq);

                computePostCollisionDistributions(currentCell, tau, feq);
            }
        }
    }
}

