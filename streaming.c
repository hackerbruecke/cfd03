#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField, int *flagField,
		int *xlength) {

	int dx, dy, dz;
	double fi;
	/*Setting distribution function for each moving direction/lattice velocity of every particle*/
	for (int z = 0; z < xlength[2] + 2; ++z) {
		for (int y = 0; y < xlength[1] + 2; ++y) {
			for (int x = 0; x < xlength[0] + 2; ++x) {

				/*Check for the cell to be FLUID*/
				if (flagField[fidx(xlength, x, y, z)] == FLUID) {
					for (int i = 0; i < Q; ++i) {

						/*dx = c_i_x*dt, dt = 1*/
						dx = LATTICEVELOCITIES[i][0];
						dy = LATTICEVELOCITIES[i][1];
						dz = LATTICEVELOCITIES[i][2];

						/*New value for our distribution function (DF) of the index 'i'

						 (We set it to DF(i) of the next particle, whose i-th lattice velocity
						 points towards considered particle (x,y,z))

						 Position of that next particle is given by (x-dx, y-dy, z-dz)*/

						fi = collideField[idx(xlength, x - dx, y - dy, z - dz,
								i)];
						streamField[idx(xlength, x, y, z, i)] = fi;
					}
				}
			}
		}
	}
}
