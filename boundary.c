#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>

int velocityIndex(int u, int v, int w) {
	return w == 0 ? (6 + (v + 1) * 3 + u) : ((w + 1) * 7 + (v + 1) * 2 + u);
}

void treatBoundary(double *collideField, int* flagField,
		const double * const wallVelocity, int* xlength, char *problem,
		double denref, double *dense) {
	int i, dx, dy, dz;
	int Nx = xlength[0];
	int Ny = xlength[1];
	int Nz = xlength[2]; /*substitution for the sake of simplicity*/
	double *currentCell;

	double finv = 0;

	/*Going through the whole boundary domain*/

	/*For each boundary cell/particle we update all the distribution functions,
	 that point towards fluid, according to boundary conditions*/

	for (int z = 0; z < xlength[2] + 2; ++z) {
		for (int y = 0; y < xlength[1] + 2; ++y) {
			for (int x = 0; x < xlength[0] + 2; ++x) {

				/* Conditions for NO_SLIP */
				if (flagField[fidx(xlength, x, y, z)] == NO_SLIP) {
					for (i = 0; i < Q; ++i) {

						dx = LATTICEVELOCITIES[i][0];
						dy = LATTICEVELOCITIES[i][1];
						dz = LATTICEVELOCITIES[i][2];

						if (x + dx > 0 && x + dx < Nx + 1 && y + dy > 0
								&& y + dy < Ny + 1 && z + dz > 0
								&& z + dz < Nz + 1) {

							/*Checking if the i-th direction is facing the fluid particles*/

							if (flagField[fidx(xlength, x + dx, y + dy, z + dz)]
									== FLUID) {

								finv = collideField[idx(xlength, x + dx, y + dy,
										z + dz, Q - i - 1)];
								collideField[idx(xlength, x, y, z, i)] = finv;
							}
						}
					}/* Conditions for OUTFLOW */
				} else if (flagField[fidx(xlength, x, y, z)] == OUTFLOW) {

					denref = 1.0;
					double velocity[D];
					double feq[Q];

					for (i = 0; i < Q; ++i) {

						dx = LATTICEVELOCITIES[i][0];
						dy = LATTICEVELOCITIES[i][1];
						dz = LATTICEVELOCITIES[i][2];

						if (x + dx > 0 && x + dx < Nx + 1 && y + dy > 0
								&& y + dy < Ny + 1 && z + dz > 0
								&& z + dz < Nz + 1) {

							/*Checking if the i-th direction is facing the fluid particles*/

							if (flagField[fidx(xlength, x + dx, y + dy, z + dz)]
									== FLUID) {

								currentCell = &collideField[idx(xlength, x + dx,
										y + dy, z + dz, 0)];

								computeVelocity(currentCell, &denref, velocity);
								computeFeq(&denref, velocity, feq);

								finv = feq[i] + feq[Q - i - 1]
										- collideField[idx(xlength, x + dx,
												y + dy, z + dz, Q - i - 1)];

								collideField[idx(xlength, x, y, z, i)] = finv;
							}
						}
					}/* Conditions for INFLOW */
				} else if (flagField[fidx(xlength, x, y, z)] == INFLOW) {

					double feq[Q];

					for (i = 0; i < Q; ++i) {

						dx = LATTICEVELOCITIES[i][0];
						dy = LATTICEVELOCITIES[i][1];
						dz = LATTICEVELOCITIES[i][2];

						if (x + dx > 0 && x + dx < Nx + 1 && y + dy > 0
								&& y + dy < Ny + 1 && z + dz > 0
								&& z + dz < Nz + 1) {

							/*Checking if the i-th direction is facing the fluid particles*/

							if (flagField[fidx(xlength, x + dx, y + dy, z + dz)]
									== FLUID) {

								computeFeq(&denref, wallVelocity, feq);

								finv = feq[i];
								collideField[idx(xlength, x, y, z, i)] = finv;
							}
						}
					}
				}
				/* Conditions for PRESSURE_IN */
				else if (flagField[fidx(xlength, x, y, z)] == PRESSURE_IN) {

					/*dense = 1.005;*/
					double velocity[D];
					double feq[Q];

					for (i = 0; i < Q; ++i) {

						dx = LATTICEVELOCITIES[i][0];
						dy = LATTICEVELOCITIES[i][1];
						dz = LATTICEVELOCITIES[i][2];

						if (x + dx > 0 && x + dx < Nx + 1 && y + dy > 0
								&& y + dy < Ny + 1 && z + dz > 0
								&& z + dz < Nz + 1) {

							/*Checking if the i-th direction is facing the fluid particles*/

							if (flagField[fidx(xlength, x + dx, y + dy, z + dz)]
									== FLUID) {

								currentCell = &collideField[idx(xlength, x + dx,
										y + dy, z + dz, 0)];
								computeVelocity(currentCell, dense, velocity);
								computeFeq(dense, velocity, feq);

								finv = feq[i] + feq[Q - i - 1]
										- collideField[idx(xlength, x + dx,
												y + dy, z + dz, Q - i - 1)];

								collideField[idx(xlength, x, y, z, i)] = finv;
							}
						}
					}
				}
				/* Conditions for FREE_SLIP */
				else if (flagField[fidx(xlength, x, y, z)] == FREE_SLIP) {

					for (i = 0; i < Q; ++i) {

						dx = LATTICEVELOCITIES[i][0];
						dy = LATTICEVELOCITIES[i][1];
						dz = LATTICEVELOCITIES[i][2];

						if (x + dx > 0 && x + dx < Nx + 1 && y + dy > 0
								&& y + dy < Ny + 1 && z + dz > 0
								&& z + dz < Nz + 1) {

							/*Checking if the i-th direction is facing the fluid particles*/

							if (flagField[fidx(xlength, x + dx, y + dy, z + dz)]
									== FLUID) {
								if (flagField[fidx(xlength, x + dx, y, z)]
										== FLUID) {
									collideField[idx(xlength, x, y, z, i)] =
											collideField[idx(xlength, x + dx, y,
													z,
													velocityIndex(-dx, dy, dz))];
								} else if (flagField[fidx(xlength, x, y + dy, z)]
										== FLUID) {
									collideField[idx(xlength, x, y, z, i)] =
											collideField[idx(xlength, x, y + dy,
													z,
													velocityIndex(dx, -dy, dz))];
								} else if (flagField[fidx(xlength, x, y, z + dz)]
										== FLUID) {
									collideField[idx(xlength, x, y, z, i)] =
											collideField[idx(xlength, x, y,
													z + dz,
													velocityIndex(dx, dy, -dz))];
								} else {
									collideField[idx(xlength, x, y, z, i)] =
											collideField[idx(xlength, x + dx,
													y + dy, z + dz, Q - i - 1)];
								}
							}
						}
					}
				}
			}
		}
	}
}
