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

	double rho = 0;
	double *currentCell;
	int dx, dy, dz;

	/* Inflow: TODO: Read values from file!!! */
	const double rho_ref = 1.0; /* rho_ref */
	double v_in[3] = { 0, 0, 0.1 }; /* V_in */
	double feq_in[Q];
	/* Pressure in: TODO: Read value from file!!! */
	const double rho_in = 1.005;
	/* Outflow */
	double v[3];
	double feq_out[Q];

	/* XY plane */
	for (int x = 0; x < xlength[0] + 2; ++x) {
		for (int y = 0; y < xlength[1] + 2; ++y) {
			for (int z = 0; z < xlength[2] + 2; ++z) {
				/* Ignore fluid cells */
				if (flagField[fidx(xlength, x, y, z)] == FLUID)
					continue;
				/* Take first distribution function of current cell */
				currentCell = &collideField[idx(xlength, x, y, z, 0)];
				/* Outflow/Moving wall conditions */
				computeDensity(currentCell, &rho);
				/* Outflow conditions */
				computeVelocity(currentCell, &rho_ref, v);
				computeFeq(&rho_ref, v, feq_out);

				for (int i = 0; i < Q; ++i) {
					dx = LATTICEVELOCITIES[i][0];
					dy = LATTICEVELOCITIES[i][1];
					dz = LATTICEVELOCITIES[i][2];

					if (inb(x+dx, y+dy, z+dz, xmax, ymax, zmax) && flagField[fidx(xlength, x+dx, y+dy, z+dz)] == FLUID) {
						double finv = collideField[idx(xlength, x+dx, y+dy, z+dz, inv(i))];

						switch (flagField[fidx(xlength, x, y, z)]) {
						case NO_SLIP:
						{
							collideField[idx(xlength, x, y, z, i)] = finv;
						}
						break;
						case MOVING_WALL:
						{
							/* Moving wall */
							double cdotu = 0.0;
							for (int d = 0; d < D; ++d) {
								cdotu += LATTICEVELOCITIES[i][d] * wallVelocity[d];
							}
							/* End moving wall */
							collideField[idx(xlength, x, y, z, i)] = finv + 2*LATTICEWEIGHTS[i]*rho*cdotu/(C_S*C_S);
						}
						break;
						case FREE_SLIP:
						{
							/* Normal inverting for 2, 6, 8, 10, 12, 16 */
							switch (i) {
							case 2: /* DOWN */
								collideField[idx(xlength, x, y, z, 0)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 14)];
								collideField[idx(xlength, x, y, z, 1)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 15)];
								collideField[idx(xlength, x, y, z, 3)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 17)];
								collideField[idx(xlength, x, y, z, 4)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 18)];
								break;
							case 6: /* FRONT */
								collideField[idx(xlength, x, y, z, 0)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 4)];
								collideField[idx(xlength, x, y, z, 5)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 11)];
								collideField[idx(xlength, x, y, z, 7)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 13)];
								collideField[idx(xlength, x, y, z, 14)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 18)];
								break;
							case 8: /* LEFT */
								collideField[idx(xlength, x, y, z, 1)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 3)];
								collideField[idx(xlength, x, y, z, 5)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 7)];
								collideField[idx(xlength, x, y, z, 11)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 13)];
								collideField[idx(xlength, x, y, z, 15)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 17)];
								break;
							case 10: /* RIGHT */
								collideField[idx(xlength, x, y, z, 3)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 1)];
								collideField[idx(xlength, x, y, z, 7)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 5)];
								collideField[idx(xlength, x, y, z, 13)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 11)];
								collideField[idx(xlength, x, y, z, 17)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 15)];
								break;
							case 12: /* BACK */
								collideField[idx(xlength, x, y, z, 4)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 0)];
								collideField[idx(xlength, x, y, z, 11)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 5)];
								collideField[idx(xlength, x, y, z, 13)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 7)];
								collideField[idx(xlength, x, y, z, 18)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 14)];
								break;
							case 16: /* TOP */
								collideField[idx(xlength, x, y, z, 14)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 0)];
								collideField[idx(xlength, x, y, z, 15)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 1)];
								collideField[idx(xlength, x, y, z, 17)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 3)];
								collideField[idx(xlength, x, y, z, 18)] = collideField[idx(xlength, x+dx, y+dy, z+dz, 4)];
								break;
							}
						}
						break;
						case INFLOW:
						{
							/* Inflow conditions */
							computeFeq(&rho_ref, v_in, feq_in);
							collideField[idx(xlength, x, y, z, i)] = feq_in[i];
						}
						break;
						case OUTFLOW:
						{
							collideField[idx(xlength, x, y, z, i)] = feq_out[inv(i)] + feq_out[i] - finv;
						}
						break;
						case PRESSURE_IN:
							/* Pressure in conditions TODO: Crosscheck! */
							computeVelocity(currentCell, &rho_in, v_in);
							computeFeq(&rho_in, v_in, feq_in);
							collideField[idx(xlength, x, y, z, i)] = feq_in[inv(i)] + feq_in[i] - finv;
							break;
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
