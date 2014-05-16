#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, int* flagField,
		const double * const wallVelocity, int* xlength, char *problem) {
	int i, k, x, y, z, dx, dy, dz;
	int Nx = xlength[0];
	int Ny = xlength[1];
	int Nz = xlength[2]; /*substitution for the sake of simplicity*/
	double *currentCell;

	double finv, cu = 0;
	double density;
	double velocity[D];
	double feq[Q];

	/*Going through the whole boundary domain*/

	/*For each boundary cell/particle we update all the distribution functions,
	 that point towards fluid, according to no-slip and moving-wall boundary conditions*/

	/*Two opposite/parallel boundary planes*/

if(*problem == 'a')
{
	for (k = 0; k < 2; ++k) {
		x = k * (Nx + 1);
		for (y = 0; y < Ny + 2; ++y) {
			for (z = 0; z < Nz + 1; ++z) {
				for (i = 0; i < Q; ++i) {
					/*dx = c_i_x*dt, dt = 1*/
					dx = LATTICEVELOCITIES[i][0];
					dy = LATTICEVELOCITIES[i][1];
					dz = LATTICEVELOCITIES[i][2];

					/*Checking if the i-th directions is facing the fluid particles*/
					/*y-z plane*/
					if (x + dx > 0 && x + dx < Nx + 1 && y + dy > 0
							&& y + dy < Ny + 1 && z + dz > 0
							&& z + dz < Nz + 1) {

						/* PLEASE CHANGE description ----- */

						/*i-th distribution fun of boundary cell becomes distribution fun of
						 the inverse lattice velocity of the pointed particle/lattice */

						/*Inverse lattice velocity to c_i is accessed
						 at Q-i-1-th position of (x+dx, y+dy, z+dz) cell*/

						finv = collideField[Q
						* ((z + dz) * (Nz + 2) * (Nz + 2)
								+ (y + dy) * (Ny + 2) + x + dx) + Q
						- i - 1];
						collideField[Q
						* (z * (Nz + 2) * (Nz + 2) + y * (Ny + 2)
								+ x) + i] = finv;
					}

					/* APPLY FREE SLIP HERE */
					//*x and y swapped, so that another pair of boundary planes may be accessed within the same iteration*/
					/*x-z plane*/
					if (y + dx > 0 && y + dx < Ny + 1 && x + dy > 0
							&& x + dy < Nx + 1 && z + dz > 0
							&& z + dz < Nz + 1) {

						/*finv = collideField[Q * ((z+dz)*(Nz+2)*(Nz+2) + (x+dy)*(Ny+2) + y+dx) + Q-i-1];
						 collideField[Q * (z*(Nz+2)*(Nz+2) + x*(Ny+2) + y) + i] = finv; */
					}
					/*x and z swapped, out of same reason*/
					/*x-y plane*/
					if (z + dx > 0 && z + dx < Nz + 1 && y + dy > 0
							&& y + dy < Ny + 1 && x + dz > 0
							&& x + dz < Nx + 1) {

						currentCell =
						&collideField[Q
						* ((x + dz) * (Nz + 2) * (Nz + 2)
								+ (y + dy) * (Ny + 2) + z
								+ dx) + i];

						/*Updating values for velocity, density and Feq for Current cell*/
						computeDensity(currentCell, &density);
						computeVelocity(currentCell, &density, velocity);
						computeFeq(&density, velocity, feq);

						if (k == 0) {
							finv = feq[i];

							collideField[Q
							* (x * (Nz + 2) * (Nz + 2)
									+ y * (Ny + 2) + z) + i] = finv;
						}

						if (k == 1) {
							finv = feq[i] + feq[Q - i - 1]
							- collideField[Q
							* ((x + dz) * (Nz + 2)
									* (Nz + 2)
									+ (y + dy) * (Ny + 2)
									+ z + dx) + Q - i - 1];

							collideField[Q
							* (x * (Nz + 2) * (Nz + 2)
									+ y * (Ny + 2) + z) + i] = finv;
						}
					}
				}
			}
		}
	}
	if(*problem == 'b')
	{
		for (k = 0; k < 2; ++k) {
			x = k * (Nx + 1);
			for (y = 0; y < Ny + 2; ++y) {
				for (z = 0; z < Nz + 1; ++z) {
					for (i = 0; i < Q; ++i) {
						/*dx = c_i_x*dt, dt = 1*/
						dx = LATTICEVELOCITIES[i][0];
						dy = LATTICEVELOCITIES[i][1];
						dz = LATTICEVELOCITIES[i][2];

						/*Checking if the i-th directions is facing the fluid particles*/
						/*y-z plane*/
						if (x + dx > 0 && x + dx < Nx + 1 && y + dy > 0
								&& y + dy < Ny + 1 && z + dz > 0
								&& z + dz < Nz + 1) {

							/* PLEASE CHANGE description ----- */

							/*i-th distribution fun of boundary cell becomes distribution fun of
							 the inverse lattice velocity of the pointed particle/lattice */

							/*Inverse lattice velocity to c_i is accessed
							 at Q-i-1-th position of (x+dx, y+dy, z+dz) cell*/

							finv = collideField[Q
							* ((z + dz) * (Nz + 2) * (Nz + 2)
									+ (y + dy) * (Ny + 2) + x + dx)
							+ Q - i - 1];
							collideField[Q
							* (z * (Nz + 2) * (Nz + 2)
									+ y * (Ny + 2) + x) + i] = finv;
						}

						/* APPLY FREE SLIP HERE */
						//*x and y swapped, so that another pair of boundary planes may be accessed within the same iteration*/
						/*x-z plane*/
						if (y + dx > 0 && y + dx < Ny + 1 && x + dy > 0
								&& x + dy < Nx + 1 && z + dz > 0
								&& z + dz < Nz + 1) {

							/*finv = collideField[Q * ((z+dz)*(Nz+2)*(Nz+2) + (x+dy)*(Ny+2) + y+dx) + Q-i-1];
							 collideField[Q * (z*(Nz+2)*(Nz+2) + x*(Ny+2) + y) + i] = finv; */
						}
						/*x and z swapped, out of same reason*/
						/*x-y plane*/
						if (z + dx > 0 && z + dx < Nz + 1 && y + dy > 0
								&& y + dy < Ny + 1 && x + dz > 0
								&& x + dz < Nx + 1) {

							currentCell = &collideField[Q
							* ((x + dz) * (Nz + 2) * (Nz + 2)
									+ (y + dy) * (Ny + 2) + z + dx)
							+ i];

							/*Updating values for velocity, density and Feq for Current cell*/
							computeDensity(currentCell, &density);
							computeVelocity(currentCell, &density,
									velocity);
							computeFeq(&density, velocity, feq);

							if (k == 0) {
								/* PRESSURE_IN CONDITIONS TO BE APPLIED */

								finv = feq[i];

								collideField[Q
								* (x * (Nz + 2) * (Nz + 2)
										+ y * (Ny + 2) + z) + i] =
								finv;
							}

							if (k == 1) {
								finv = feq[i] + feq[Q - i - 1]
								- collideField[Q
								* ((x + dz) * (Nz + 2)
										* (Nz + 2)
										+ (y + dy)
										* (Ny + 2)
										+ z + dx) + Q - i
								- 1];

								collideField[Q
								* (x * (Nz + 2) * (Nz + 2)
										+ y * (Ny + 2) + z) + i] =
								finv;
							}
						}
					}
				}
			}
		}

		if(*problem == 'c')
		{
			for (k = 0; k < 2; ++k) {
				x = k * (Nx + 1);
				for (y = 0; y < Ny + 2; ++y) {
					for (z = 0; z < Nz + 1; ++z) {
						for (i = 0; i < Q; ++i) {
							/*dx = c_i_x*dt, dt = 1*/
							dx = LATTICEVELOCITIES[i][0];
							dy = LATTICEVELOCITIES[i][1];
							dz = LATTICEVELOCITIES[i][2];

							/*Checking if the i-th directions is facing the fluid particles*/
							/*y-z plane*/
							if (x + dx > 0 && x + dx < Nx + 1 && y + dy > 0
									&& y + dy < Ny + 1 && z + dz > 0
									&& z + dz < Nz + 1) {

								/* PLEASE CHANGE description ----- */

								/*i-th distribution fun of boundary cell becomes distribution fun of
								 the inverse lattice velocity of the pointed particle/lattice */

								/*Inverse lattice velocity to c_i is accessed
								 at Q-i-1-th position of (x+dx, y+dy, z+dz) cell*/

								finv = collideField[Q
								* ((z + dz) * (Nz + 2) * (Nz + 2)
										+ (y + dy) * (Ny + 2) + x
										+ dx) + Q - i - 1];
								collideField[Q
								* (z * (Nz + 2) * (Nz + 2)
										+ y * (Ny + 2) + x) + i] =
								finv;
							}

							/* APPLY FREE SLIP HERE */
							//*x and y swapped, so that another pair of boundary planes may be accessed within the same iteration*/
							/*x-z plane*/
							if (y + dx > 0 && y + dx < Ny + 1 && x + dy > 0
									&& x + dy < Nx + 1 && z + dz > 0
									&& z + dz < Nz + 1) {

								/*finv = collideField[Q * ((z+dz)*(Nz+2)*(Nz+2) + (x+dy)*(Ny+2) + y+dx) + Q-i-1];
								 collideField[Q * (z*(Nz+2)*(Nz+2) + x*(Ny+2) + y) + i] = finv; */
							}
							/*x and z swapped, out of same reason*/
							/*x-y plane*/
							if (z + dx > 0 && z + dx < Nz + 1 && y + dy > 0
									&& y + dy < Ny + 1 && x + dz > 0
									&& x + dz < Nx + 1) {

								currentCell = &collideField[Q
								* ((x + dz) * (Nz + 2) * (Nz + 2)
										+ (y + dy) * (Ny + 2) + z
										+ dx) + i];

								/*Updating values for velocity, density and Feq for Current cell*/
								computeDensity(currentCell, &density);
								computeVelocity(currentCell, &density,
										velocity);
								computeFeq(&density, velocity, feq);

								if (k == 0) {
									finv = feq[i];

									collideField[Q
									* (x * (Nz + 2) * (Nz + 2)
											+ y * (Ny + 2) + z) + i] =
									finv;
								}

								if (k == 1) {
									finv =
									feq[i] + feq[Q - i - 1]
									- collideField[Q* ((x + dz)	* (Nz+ 2)* (Nz+ 2)+ (y+ dy)	* (Ny+ 2)+ z + dx)+ Q - i - 1];

									collideField[Q* (x * (Nz + 2) * (Nz + 2)+ y * (Ny + 2) + z) + i] = finv;
								}
							}
						}
					}
				}
			}
