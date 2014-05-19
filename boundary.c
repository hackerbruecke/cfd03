#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>

#if 1

static inline int inb(int x, int y, int z, int xmax, int ymax, int zmax) {
	return  x > 0 && x < xmax &&
			y > 0 && y < ymax &&
			z > 0 && z < zmax;
}

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, const int *xlength) {
	/* XY plane */
	int xmax = xlength[0] + 1;
	int ymax = xlength[1] + 1;
	int zmax = xlength[2] + 1;

	double density = 0;
	double *currentCell;
	int dx, dy, dz;

	/* XY plane */
	for (int x = 0; x < xlength[0] + 2; ++x) {
		for (int y = 0; y < xlength[1] + 2; ++y) {
			for (int z = 0; z < xlength[2] + 2; ++z) {
				for (int i = 0; i < Q; ++i) {
					/* Bottom, z = 0: Invert 0, 1, 2, 3, 4 */
					/*index = plane_fi[XY0][i];*/
					dx = LATTICEVELOCITIES[i][0];
					dy = LATTICEVELOCITIES[i][1];
					dz = LATTICEVELOCITIES[i][2];
					if (inb(x+dx, y+dy, z+dz, xmax, ymax, zmax) && flagField[fidx(xlength, x+dx, y+dy, z+dz)] == FLUID) {
/*					if (inb(x, y, z, xmax, ymax, zmax) && flagField[fidx(xlength, x+dx, y+dy, z+dz)] == FLUID) {*/
						switch (flagField[fidx(xlength, x, y, z)]) {
						case NO_SLIP:
						{
							collideField[idx(xlength, x, y, z, i)] = collideField[idx(xlength, x+dx, y+dy, z+dz, inv(i))];
						}
						break;
						case MOVING_WALL:
						{
							/* Moving wall */
							double cdotu = 0.0;
							for (int d = 0; d < D; ++d) {
								cdotu += LATTICEVELOCITIES[i][d] * wallVelocity[d];
							}
							/* Take first distribution function of current cell */
							currentCell = &collideField[idx(xlength, x, y, z, 0)];
							computeDensity(currentCell, &density);
							/* End moving wall */
							double finv = collideField[idx(xlength, x+dx, y+dy, z+dz, inv(i))];
							collideField[idx(xlength, x, y, z, i)] = finv + 2*LATTICEWEIGHTS[i]*density*cdotu/(C_S*C_S);
						}
						break;
						case FLUID:
							break;
						default:
							printf("None\n");
						}
					}
				}
			}
		}
	}
}
#else
void treatBoundary(double *collideField, int* flagField,
		const double * const inputParameters, const int *xlength) {
	int i, x, y, z, dx, dy, dz;
	int nx = xlength[0]; /*substitution for the sake of simplicity*/
	int ny = xlength[1];
	int nz = xlength[2];
	double finv;
	double cu = 0;
	double density;
	double *currentCell;

	/*Going through the hole boundary domain*/

	/*For each boundary cell/particle we update all the distribution functions,
	that point towards fluid, according to no-slip and moving-wall boundary conditions*/

	/*Two opposite/parallel boundary planes*/
	for (int k = 0; k < 2; ++k) {
		x = k * (nx + 1);
		for (y = 0; y < ny + 2; ++y) {
			for (z = 0; z < nz + 1; ++z) {
				for (i = 0; i < Q; ++i) {
					/*dx = c_i_x*dt, dt = 1*/
					dx = LATTICEVELOCITIES[i][0];
					dy = LATTICEVELOCITIES[i][1];
					dz = LATTICEVELOCITIES[i][2];

					/*Checking if the i-th directions is facing the fluid particles*/
					/*y-z plane*/
					if (x+dx > 0 && x+dx < nx + 1 &&
						y+dy > 0 && y+dy < ny + 1 &&
						z+dz > 0 && z+dz < nz + 1) {

						/*i-th distribution fun of boundary cell becomes distribution fun of
						the inverse lattice velocity of the pointed particle/lattice */

						/*Inverse lattice velocity to c_i is accessed
						at Q-i-1-th position of (x+dx, y+dy, z+dz) cell*/
						finv = collideField[idx(xlength, x+dx, y+dy, z+dz, inv(i))];
						collideField[idx(xlength, x, y, z, i)] = finv;
					}
					/*x and y swapped, so that another pair of boundary planes may be accessed within the same iteration*/
					/*x-z plane*/
					if (y+dx > 0 && y+dx < nx + 1 &&
						x+dy > 0 && x+dy < ny + 1 &&
						z+dz > 0 && z+dz < nz + 1) {

						finv = collideField[idx(xlength, y+dx, x+dy, z+dz, inv(i))];
						collideField[idx(xlength, y, x, z, i)] = finv;
					}
					/*x and z swapped, out of same reason*/
					/*x-y plane*/
					if (z+dx > 0 && z+dx < nx + 1 &&
						y+dy > 0 && y+dy < ny + 1 &&
						x+dz > 0 && x+dz < nz + 1) {

						finv = collideField[idx(xlength, z+dx, y+dy, x+dz, inv(i))];
					#if 1 /* No wall velocity*/
						/*Additional term when we are considering moving-wall (top plane - k=1)*/
						/*Velocity of the moving-wall is taken into account*/
						if(k == 1) {
							currentCell = collideField + idx(xlength, z+dx, y+dy, x+dz, 0);
							computeDensity(currentCell, &density);
							cu = 0; /*dot product of c_i and velocity of wall*/
							for (int d = 0; d < D; ++d) {
								cu += LATTICEVELOCITIES[i][d] * inputParameters[d];
							}
							finv += 2*LATTICEWEIGHTS[i]*density*cu/(C_S*C_S);
						}
					#endif
						collideField[idx(xlength, z, y, x, i)] = finv;
					}
				}
			}
		}
	}
}

#endif
