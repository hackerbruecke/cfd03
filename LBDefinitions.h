#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#define Q 19
#define D 3

static const int LATTICEVELOCITIES[19][3] =
{
    {  0, -1, -1 }, { -1,  0, -1 }, {  0,  0, -1 }, {  1,  0, -1 }, {  0,  1, -1 },
    { -1, -1,  0 }, {  0, -1,  0 }, {  1, -1,  0 }, { -1,  0,  0 }, {  0,  0,  0 },
    {  1,  0,  0 }, { -1,  1,  0 }, {  0,  1,  0 }, {  1,  1,  0 }, {  0, -1,  1 },
    { -1,  0,  1 }, {  0,  0,  1 }, {  1,  0,  1 }, {  0,  1,  1 }
};
static const double LATTICEWEIGHTS[19] = 
{
    1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36,
    1.0/36, 2.0/36, 1.0/36, 2.0/36, 12.0/36,
    2.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36,
    1.0/36, 2.0/36, 1.0/36, 1.0/36
};
static const double C_S = 0.57735026919l; 

typedef enum {
	FLUID = 0,
	NO_SLIP = 1,
	MOVING_WALL = 2,
	FREE_SLIP = 3,
	INFLOW = 4,
	OUTFLOW = 5,
	PRESSURE_IN = 6
} STATE;
#if 0 /* TODO: Remove when enum works */
static const int FLUID = 0;
static const int NO_SLIP = 1;
static const int MOVING_WALL = 2;
#endif

static inline int idx(const int *xlength, int x, int y, int z, int i) {
    return Q * (z * (xlength[0]+2) * (xlength[1]+2) + y * (xlength[0]+2) + x) + i;
}

static inline int fidx(const int *xlength, int x, int y, int z) {
    return z * (xlength[0]+2) * (xlength[1]+2) + y * (xlength[0]+2) + x;
}

static inline int inv(int idx) {
    return 18 - idx;
}
#endif

