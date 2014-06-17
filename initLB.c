#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall,
		int *timesteps, int *timestepsPerPlotting, char *problem, double denref,
		double *dense, int argc, char *argv[]) {
	/* Checking if there are exactly 3 parameters in the input. Else give Error message. */
	if (argc != 3) {
		printf("Please input Data File and Problem character :\n"
				"1 - TiltedPlate\n"
				"2 - PlaneShear\n"
				"3 - FlowOverStep\n"
				"Eg: ./lbsim problem.dat 1\n");
		return -1;
	}

	const char *filename = argv[1];

	/* Reading all the input parameters */

	if (*problem == '1') {
		read_int(filename, "xlength1", &xlength[0]);
		read_int(filename, "ylength1", &xlength[1]);
		read_int(filename, "zlength1", &xlength[2]);
		read_int(filename, "timesteps1", timesteps);

		read_double(filename, "velocityWallx1", &velocityWall[0]);
		read_double(filename, "velocityWally1", &velocityWall[1]);
		read_double(filename, "velocityWallz1", &velocityWall[2]);
		read_double(filename, "tau1", tau);
		READ_INT(filename, *timestepsPerPlotting);
		return 0;
	} else if (*problem == '2') {
		read_int(filename, "xlength2", &xlength[0]);
		read_int(filename, "ylength2", &xlength[1]);
		read_int(filename, "zlength2", &xlength[2]);
		read_int(filename, "timesteps2", timesteps);

		read_double(filename, "denref", &denref);
		read_double(filename, "dense", dense);
		read_double(filename, "tau2", tau);
		READ_INT(filename, *timestepsPerPlotting);
		return 0;
	} else if (*problem == '3') {
		read_int(filename, "xlength3", &xlength[0]);
		read_int(filename, "ylength3", &xlength[1]);
		read_int(filename, "zlength3", &xlength[2]);
		read_int(filename, "timesteps3", timesteps);

		read_double(filename, "velocityWallx3", &velocityWall[0]);
		read_double(filename, "velocityWally3", &velocityWall[1]);
		read_double(filename, "velocityWallz3", &velocityWall[2]);
		read_double(filename, "tau3", tau);
		READ_INT(filename, *timestepsPerPlotting);
	}

	return 0;
}

void initialiseFields(double *collideField, double *streamField, int *flagField,
		int *xlength, char *problem) {
	/* Set all values of flagField to 0, i.e. FLUID state. We will apply boundary conditions
	 in the following */
	memset(flagField, FLUID,
			(xlength[2] + 2) * (xlength[1] + 2) * (xlength[0] + 2)
					* sizeof(*flagField));

	/*Scenario for Problem 1 (As given in Worksheet) */
	if (*problem == '1') {

		/* The values for Boundary Z = 0 set to Inflow and Zmax plane set to Outflow */
		for (int x = 0; x < xlength[0] + 2; ++x) {
			for (int y = 0; y < xlength[1] + 2; ++y) {
				flagField[fidx(xlength, x, y, 0)] = INFLOW;
				flagField[fidx(xlength, x, y, xlength[2] + 1)] = OUTFLOW;
			}
		}
		/* The values for Boundary on X = 0 and Xmax plane set to NO_SLIP */
		for (int y = 0; y < xlength[1] + 2; ++y) {
			for (int z = 0; z < xlength[2] + 2; ++z) {
				flagField[fidx(xlength, 0, y, z)] = NO_SLIP;
				flagField[fidx(xlength, xlength[0] + 1, y, z)] = NO_SLIP;
			}
		}
		/* The values for Boundary on Y = 0 and Ymax plane set to FREE_SLIP */
		for (int x = 0; x < xlength[0] + 2; ++x) {
			for (int z = 0; z < xlength[2] + 2; ++z) {
				flagField[fidx(xlength, x, 0, z)] = FREE_SLIP;
				flagField[fidx(xlength, x, xlength[1] + 1, z)] = FREE_SLIP;
			}
		}

		/* Titled Plate geometry read from the .vtk file */
		int **tilt;
		tilt = read_pgm("lbm_tilted_plate.vtk");

		for (int z = 1; z < xlength[2] + 1; ++z) {
			for (int y = 0; y < xlength[1] + 1; ++y) {
				for (int x = 1; x < xlength[0] + 1; ++x) {
					flagField[fidx(xlength, x, y, z)] =
							tilt[xlength[0] + 1 - x][xlength[2] + 1 - z];
				}
			}
		}
		free(tilt);
	}
	/*Scenario for Problem 2 (As given in Worksheet) */
	else if (*problem == '2') {
		/* The values for Boundary Z = 0 set to PressureIn and Zmax plane set to Outflow */
		for (int x = 0; x < xlength[0] + 2; ++x) {
			for (int y = 0; y < xlength[1] + 2; ++y) {
				flagField[fidx(xlength, x, y, 0)] = PRESSURE_IN;
				flagField[fidx(xlength, x, y, xlength[2] + 1)] = OUTFLOW;
			}
		}
		/* The values for Boundary on X = 0 and Xmax plane set to NO_SLIP */
		for (int y = 0; y < xlength[1] + 2; ++y) {
			for (int z = 0; z < xlength[2] + 2; ++z) {
				flagField[fidx(xlength, 0, y, z)] = FREE_SLIP;
				flagField[fidx(xlength, xlength[0] + 1, y, z)] = FREE_SLIP;
			}
		}
		/* The values for Boundary on Y = 0 and Ymax plane set to FREE_SLIP */
		for (int x = 0; x < xlength[0] + 2; ++x) {
			for (int z = 0; z < xlength[2] + 2; ++z) {
				flagField[fidx(xlength, x, 0, z)] = NO_SLIP;
				flagField[fidx(xlength, x, xlength[1] + 1, z)] = NO_SLIP;
			}
		}

	}
	/*Scenario for Problem 3 (As given in Worksheet) */
	else if (*problem == '3') {
		/* The values for Boundary Z = 0 on INFLOW and Zmax plane set to OUTFLOW */
		for (int y = 0; y < xlength[1] + 2; ++y) {
			for (int x = 0; x < xlength[0] + 2; ++x) {
				flagField[fidx(xlength, x, y, 0)] = INFLOW;
				flagField[fidx(xlength, x, y, xlength[2] + 1)] = OUTFLOW;
			}
		}
		/* The values for Boundary on X = 0 and Xmax plane set to NO_SLIP */
		for (int z = 0; z < xlength[2] + 2; ++z) {
			for (int y = 0; y < xlength[1] + 2; ++y) {
				flagField[fidx(xlength, 0, y, z)] = NO_SLIP;
				flagField[fidx(xlength, xlength[0] + 1, y, z)] = NO_SLIP;
			}
		}
		/* The values for Boundary on Y = 0 and Ymax plane set to NO_SLIP */
		for (int z = 0; z < xlength[2] + 2; ++z) {
			for (int x = 0; x < xlength[0] + 2; ++x) {
				flagField[fidx(xlength, x, 0, z)] = NO_SLIP;
				flagField[fidx(xlength, x, xlength[1] + 1, z)] = NO_SLIP;
			}
		}

		/* Step Block geometry defined */
		for (int z = 0; z < (xlength[2] + 2) / 2; ++z) {
			for (int y = 0; y < xlength[1] + 2; ++y) {
				for (int x = 0; x < (xlength[0] + 2) / 2; ++x) {
					flagField[fidx(xlength, x, y, z)] = NO_SLIP;
				}
			}
		}
	}

	/* Stream and Collide Fields are initialized to the respective Lattice-Weights of the Cell */
	for (int z = 0; z < xlength[2] + 2; ++z) {
		for (int y = 0; y < xlength[1] + 2; ++y) {
			for (int x = 0; x < xlength[0] + 2; ++x) {
				for (int i = 0; i < Q; ++i) {
					collideField[idx(xlength, x, y, z, i)] = LATTICEWEIGHTS[i];
					streamField[idx(xlength, x, y, z, i)] = LATTICEWEIGHTS[i];
				}
			}
		}
	}
}

