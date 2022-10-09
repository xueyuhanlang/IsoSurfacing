/*
	File: marching_cubes_33.h
	Programmed by: David Vega - dvega@uc.edu.ve
	           Javier Abache - jabache@uc.edu.ve
	March 2012
	updated by David Vega
	December 2014,
	June 2015,
	November 2017.
	February 2018.
*/

/********************************CUSTOMIZING**********************************/
/*
GRD_data_type is the variable type of the grid data, by default it is float,
you can use unsigned char, unsigned short int, int, etc ...
*/
//Do not change the following lines after compiling the library:
#ifndef GRD_data_type
#define GRD_data_type float
#endif
#define _ORTHO_GRD
//#define MC_Normal_neg //If defined, the front and back surface are exchanged.
/*****************************************************************************/
#ifndef marching_cubes_33_h
#define marching_cubes_33_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
The structure _GRD contains a function F[][][] evaluated at a grid of regularly
spaced points. N[] is the number of intervals in each dimension, so the number
of points is N[] + 1, and the grid contain (N[2] + 1)*(N[1] + 1)*(N[0] + 1)
points. L[] is the grid size in each dimension. r0[] are the coordinates of
the first grid point. d[] is the distance between adjacent points in each
dimension (can be different for each dimension), nonortho has a value of 1 when
the grid is inclined else 0. _A and A_ are the matrices that transform from
inclined to orthogonal coordinates and vice versa, respectively. If the grid is
periodic (is infinitely repeated along each dimension) the flag periodic must
be 1 (for x) + 2 (for y) + 4 (for z), else 0.

In this program, if _ORTHO_GRD is defined, then nonortho, _A and A_ can be
removed from this structure, and only work with orthogonal grids. periodic
and L[] are not used here and also can be removed.
*/

#ifndef LAGRANGE3D4GRD_H
typedef struct
{
	GRD_data_type ***F;
	int N[3];
	float L[3], Ang[3];//not necessary
	double r0[3], d[3];
#ifndef _ORTHO_GRD
	int nonortho; //necessary if _ORTHO_GRD is not defined
	double _A[3][3], A_[3][3]; //necessary if _ORTHO_GRD is not defined
#endif
	int periodic;//not necessary
	char title[160];//not necessary
} _GRD;
#endif

/*The struct surface (pointed by S) contains the data of an isosurface. The
isovalue is S->ISO, the number of points is S->nV + 1, and the number of
triangles is S->nT + 1. The coordinates of the points are stored in the array
V[][dim2][3]. The coordinates of the point number M are obtained from the
vector V[M/S->dim2][M%S->dim2][]. N[][dim2][3] is the array that contains the
normal vectors calculated at the points of the surface. color[][dim2] contains
the color index of each point. T[][dim2][3] is the array that contains the
sets of three point (vertex) indices that form each triangle.*/
typedef struct
{
	float ISO;
	int dim2, flag;
	int nV, nT;
	int (**T)[3];
	float (**V)[3];
	float (**N)[3];
	int **color;
} surface;

#ifndef UTILMAT_H
//from https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Reciprocal_of_the_square_root
//http://www.lomont.org/Math/Papers/2003/InvSqrt.pdf (0x5f375a86)
float invSqrt(float x);
//{
//	union
//	{
//		float x;
//		int i;
//	} u;
//	u.x = x;
//	x *= 0.5f;
//	u.i = 0x5f375a86 - (u.i>>1);
//	/* The next line can be repeated any number of times to increase accuracy */
//	u.x = u.x * (1.5f - x * u.x * u.x);
//	return u.x;
//}
#endif

#if !defined(UTILMAT_H) && !defined(chgcar_cube_h) && !defined(_ORTHO_GRD) 
//c = Ab, A is a 3x3 upper triangular matrix. If t != 0, A is transposed.
void _multTSA_bf(double A[3][3], float *b, float *c, int t)
{
	if(t)
	{
		c[2] = A[0][2]*b[0] + A[1][2]*b[1] + A[2][2]*b[2];
		c[1] = A[0][1]*b[0] + A[1][1]*b[1];
		c[0] = A[0][0]*b[0];
	}
	else
	{
		c[0] = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
		c[1] = A[1][1]*b[1] + A[1][2]*b[2];
		c[2] = A[2][2]*b[2];
	}
}
//Performs the multiplication of the matrix A and the vector b: c = Ab. If t != 0, A is transposed.
void _multA_bf(double A[3][3], float* b, float* c, int t)
{
	double u,v;
	if(t)
	{
		u = A[0][0]*b[0] + A[1][0]*b[1] + A[2][0]*b[2];
		v = A[0][1]*b[0] + A[1][1]*b[1] + A[2][1]*b[2];
		c[2] = A[0][2]*b[0] + A[1][2]*b[1] + A[2][2]*b[2];
	}
	else
	{
		u = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
		v = A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2];
		c[2] = A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2];
	}
	c[0] = u;
	c[1] = v;
}
void (*mult_Abf)(double (*)[3], float *, float *, int) = _multA_bf;
#endif

/*
Sets the size of the second dimension of the arrays of the structure surface.
The value is a power of 2. The default value is 4096. This value is stored in
the global variable _MC_DVE and in the member dim2 of the structure.
*/
void set_MC_DVE(unsigned int size);
extern int _MC_DVE, _MC_A, _MC_N;
/*
DefaultColorMC = 0x00FF0080
										 B G R
Red 128, green 0, blue 255.
*/
extern int DefaultColorMC;


/******************************************************************
Saves all the surface *S data (in binary format) to a "filename" file. The
return value is 0 if the call succeeds, else -1.
*/
int write_bin_s(char *filename, surface *S);
/******************************************************************
Reads (from a "filename" file) the surface data stored in binary format.
The return value is a surface pointer if the call succeeds, else NULL.
*/
surface* read_bin_s(char *filename);
/******************************************************************
Saves all the surface *S data (in plain text format) to a "filename" file.
The return value is 0 if the call succeeds, else -1.
*/
int write_txt_s(char* filename, surface *S);
/******************************************************************
Calculates the isosurface with isovalue equal to iso, using a grid grd. The
return value is a pointer to surface. The pointer is NULL if there is not
enough memory.*/
surface* calculate_isosurface(_GRD* grd, float iso);

/*The above function calls the following ones, in this order:*/
//1. init_temp_isosurface. It initializes variables and allocates temporary structures.
int init_temp_isosurface(_GRD* grd);
//2. calc_isosurface. Use this function to calculate several isosurfaces with the same grid.
surface* calc_isosurface(float iso);
//3. clear_temp_isosurface. It frees temporary structures.
void clear_temp_isosurface();
/******************************************************************
Releases the allocated memory pointed to by S.
*/
void free_surface_memory(surface **S);

#endif //marching_cubes_33_h

