/*
	File: marching_cubes_33.c
	Programmed by: David Vega - dvega@uc.edu.ve
	           Javier Abache - jabache@uc.edu.ve
	March 2012
	updated by David Vega
	December 2014,
	June 2015,
	November 2017.
	February 2018.
	June 2018
*/

#include "MC33_LookUpTable.h"
#include "marching_cubes_33.h"

#ifndef marching_cubes_33_c
#define marching_cubes_33_c


#ifdef __cplusplus
#include <algorithm>
#endif

#ifndef default_MC_N
#define default_MC_N 12
#endif

float invSqrt(float x)
{
	union
	{
		float x;
		int i;
	} u;
	u.x = x;
	x *= 0.5f;
	u.i = 0x5f375a86 - (u.i >> 1);
	/* The next line can be repeated any number of times to increase accuracy */
	u.x = u.x * (1.5f - x * u.x * u.x);
	return u.x;
}

/* _MC_DVE and _MC_A are calculated using _MC_N.
_MC_A is a mask (n & _MC_A is equivalent to n % _MC_DVE).
n >> _MC_N is equivalent to n / _MC_DVE. */
int _MC_N = default_MC_N, _MC_DVE = 1<<default_MC_N, _MC_A = (1<<default_MC_N) - 1;

void set_MC_DVE(unsigned int size)
{
	_MC_N = 0;
	--size;
	while(size)
	{
		++_MC_N;
		size>>=1;
	}
	if(_MC_N < 6)
		_MC_N = 6;
	_MC_DVE = 1<<_MC_N;
	_MC_A = _MC_DVE - 1;
}

int DefaultColorMC = 0xff5c5c5c;//gray RGB color as unsigned char [3]
/*
_MC_memoryfault takes a non-zero value if the system memory is insufficient
to store the surface.
*/
int _MC_memoryfault;

#ifndef incDimTVN
#define incDimTVN (1<<9)
#endif

/*_MCnT and _MCnV are the value of the first dimension of arrays T and V of the
structure surface. They are used in store_point_normal and store_triangle.
functions*/
int _MCnT, _MCnV;
surface* _MC_S;

GRD_data_type ***_MC_F = 0;
float _MC_O[3], _MC_D[3];
int _MCnx, _MCny, _MCnz;
#ifndef _ORTHO_GRD
int _MC_NORT;
double (*_MC__A)[3], (*_MC_A_)[3];
#endif

void _MC_memory_error(int i)
{
	_MC_memoryfault = 1;
	if(i)
	{
		_MC_S->nT -= _MC_DVE;
		--_MCnT;
	}
	else
	{
		_MC_S->nV -= _MC_DVE;
		--_MCnV;
	}
}

// temporary structures that store the indices of triangle vertices:
int **_Ox, **_Oy, **_Nx, **_Ny, *_OL, *_NL;
/******************************************************************
Vertices:           Faces:
    3 __________2        ___________
   /|          /|      /|          /|
  / |         / |     / |   2     / |
7/__________6/  |    /  |     4  /  |
|   |       |   |   |¯¯¯¯¯¯¯¯¯¯¯| 1 |     z
|   0_______|___1   | 3 |_______|___|     |
|  /        |  /    |  /  5     |  /      |____y
| /         | /     | /     0   | /      /
4/__________5/      |/__________|/      x


This function return a vector with all six test face results (face[6]). Each
result value is 1 if the positive face vertices are joined, -1 if the negative
vertices are joined, and 0 (unchanged) if the test must no be applied. The
return value of this function is the the sum of all six results.*/
int face_tests(int *face, int ind, int sw, float *v)
{
	if(ind&0x80)//vertex 0
	{
		face[0] = ((ind&0xCC) == 0x84? (v[0]*v[5] < v[1]*v[4]? -sw: sw): 0);//0x84 = 10000100, vertices 0 and 5
		face[3] = ((ind&0x99) == 0x81? (v[0]*v[7] < v[3]*v[4]? -sw: sw): 0);//0x81 = 10000001, vertices 0 and 7
		face[4] = ((ind&0xF0) == 0xA0? (v[0]*v[2] < v[1]*v[3]? -sw: sw): 0);//0xA0 = 10100000, vertices 0 and 2
	}
	else
	{
		face[0] = ((ind&0xCC) == 0x48? (v[0]*v[5] < v[1]*v[4]? sw: -sw): 0);//0x48 = 01001000, vertices 1 and 4
		face[3] = ((ind&0x99) == 0x18? (v[0]*v[7] < v[3]*v[4]? sw: -sw): 0);//0x18 = 00011000, vertices 3 and 4
		face[4] = ((ind&0xF0) == 0x50? (v[0]*v[2] < v[1]*v[3]? sw: -sw): 0);//0x50 = 01010000, vertices 1 and 3
	}
	if(ind&0x02)//vertex 6
	{
		face[1] = ((ind&0x66) == 0x42? (v[1]*v[6] < v[2]*v[5]? -sw: sw): 0);//0x42 = 01000010, vertices 1 and 6
		face[2] = ((ind&0x33) == 0x12? (v[3]*v[6] < v[2]*v[7]? -sw: sw): 0);//0x12 = 00010010, vertices 3 and 6
		face[5] = ((ind&0x0F) == 0x0A? (v[4]*v[6] < v[5]*v[7]? -sw: sw): 0);//0x0A = 00001010, vertices 4 and 6
	}
	else
	{
		face[1] = ((ind&0x66) == 0x24? (v[1]*v[6] < v[2]*v[5]? sw: -sw): 0);//0x24 = 00100100, vertices 2 and 5
		face[2] = ((ind&0x33) == 0x21? (v[3]*v[6] < v[2]*v[7]? sw: -sw): 0);//0x21 = 00100001, vertices 2 and 7
		face[5] = ((ind&0x0F) == 0x05? (v[4]*v[6] < v[5]*v[7]? sw: -sw): 0);//0x05 = 00000101, vertices 5 and 7
	}
	return face[0] + face[1] + face[2] + face[3] + face[4] + face[5];
}
/* Faster function for the face test, the test is applied to only one face
(int face). This function is only used for the cases 3 and 6 of MC33*/
int face_test1(int face, float *v)
{
	switch(face)
	{
	case 0:
		return (v[0]*v[5] < v[1]*v[4]? 0x48: 0x84);
	case 1:
		return (v[1]*v[6] < v[2]*v[5]? 0x24: 0x42);
	case 2:
		return (v[3]*v[6] < v[2]*v[7]? 0x21: 0x12);
	case 3:
		return (v[0]*v[7] < v[3]*v[4]? 0x18: 0x81);
	case 4:
		return (v[0]*v[2] < v[1]*v[3]? 0x50: 0xA0);
	case 5:
		return (v[4]*v[6] < v[5]*v[7]? 0x05: 0x0A);
	}
	return 0;
}

#if defined(__MINGW32__)//for TDM-GCC signbit(float) bug.
#define signbf __signbitf
#else
#define signbf signbit
// #define signbf(x) (*(reinterpret_cast<unsigned int*>(&x))>>31)
#endif

/******************************************************************
Interior test function. If the test is positive, the function returns a value
different fom 0. The integer i must be 0 to test if the vertices 0 and 6 are
joined. 1 for vertices 1 and 7, 2 for vertices 2 and 4, and 3 for 3 and 5.
For case 13, the integer flagtplane must be 1, and the function returns 2 if
one of the vertices 0, 1, 2 or 3 is joined to the center point of the cube
(case 13.5.2), returns 1 if one of the vertices 4, 5, 6 or 7 is joined to the
center point of the cube (case 13.5.2 too), and it returns 0 if the vertices
are no joined (case 13.5.1)*/
int interior_test(int i, int flagtplane, float *v)
{
	//Signs of cube vertices were changed to use signbit function in calc_isosurface
	//A0 = -v[0], B0 = -v[1], C0 = -v[2], D0 = -v[3]
	//A1 = -v[4], B1 = -v[5], C1 = -v[6], D1 = -v[7]
	//But the function still works
	float At = v[4] - v[0], Bt = v[5] - v[1],
				Ct = v[6] - v[2], Dt = v[7] - v[3];
	float t = At*Ct - Bt*Dt;//the "a" value.
	if(i&0x01)//false for i = 0 and 2, and true for i = 1 and 3
	{
		if(t <= 0.0f) return 0;
	}
	else
	{
		if(t >= 0.0f) return 0;
	}
	t = 0.5f*(v[3]*Bt + v[1]*Dt - v[2]*At - v[0]*Ct)/t;//t = -b/2a
	if(t <= 0.0f || t >= 1.0f)
		return 0;

	At = v[0] + At*t;
	Bt = v[1] + Bt*t;
	Ct = v[2] + Ct*t;
	Dt = v[3] + Dt*t;
	if(i&0x01)
	{
		if(At*Ct < Bt*Dt && signbf(Bt) == signbf(Dt))
			return (signbf(Bt) == signbf(v[i])) + flagtplane;
	}
	else
	{
		if(At*Ct > Bt*Dt && signbf(At) == signbf(Ct))
			return (signbf(At) == signbf(v[i])) + flagtplane;
	}
	return 0;
}

/******************************************************************
Assign memory for the vertex r[3], normal n[3]. The return value is the new
vertex label.
*/
int store_point_normal(float *r, float *n)
{
	int i, nv = (++_MC_S->nV)&_MC_A;
	float t = 0.0f, *p;
	void *pv, *pn, *pc;
	if(!nv)//create a new _MC_S->V[*][][], and _MC_S->N[*][][], and _MC_S->color[*][]
	{
		if(!((++_MCnV)&(incDimTVN - 1)))//expand the _MC_S->V[] dimension in incDimTVN
		{
			pc = realloc(_MC_S->color,(_MCnV + incDimTVN)*sizeof(void*));
			pn = realloc(_MC_S->N,(_MCnV + incDimTVN)*sizeof(void*));
			pv = realloc(_MC_S->V,(_MCnV + incDimTVN)*sizeof(void*));
			if(pv)
			{
				_MC_S->V = (float (**)[3])pv;
				_MC_S->N = (float (**)[3])pn;
				_MC_S->color = (int **)pc;
			}
			else
			{
				if(pn)
					_MC_S->N = (float (**)[3])pn;
				if(pc)
					_MC_S->color = (int **)pc;
				_MC_memory_error(0);
				return 0;
			}
		}
		_MC_S->color[_MCnV] = (int *)malloc(_MC_DVE*sizeof(int));
		_MC_S->N[_MCnV] = (float (*)[3])malloc(3*_MC_DVE*sizeof(float));
		_MC_S->V[_MCnV] = (float (*)[3])malloc(3*_MC_DVE*sizeof(float));
		if(!_MC_S->V[_MCnV])
		{
			free(_MC_S->N[_MCnV]);
			free(_MC_S->color[_MCnV]);
			_MC_memory_error(0);
			return 0;
		}
	}
	_MC_S->color[_MCnV][nv] = DefaultColorMC;
	p = _MC_S->V[_MCnV][nv];
	for(i=0;i<3;++i)
		p[i] = r[i]*_MC_D[i] + _MC_O[i];
#ifndef _ORTHO_GRD
	if(_MC_NORT)
	{
		mult_Abf(_MC__A,p,p,0);
		mult_Abf(_MC_A_,n,n,1);
	}
#endif
	for(i=0;i<3;++i)
		t += n[i]*n[i];
#ifndef MC_Normal_neg
//	t = 1.0f/sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	t = invSqrt(t);
#else //MC_Normal_neg reverse the direction of the normal
//t = -1.0f/sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	t = -invSqrt(t);
#endif
	p = _MC_S->N[_MCnV][nv];
	for(i=0;i<3;++i)
		p[i] = t*n[i];
	return _MC_S->nV;
}

/******************************************************************
Store the triangle, an array of three vertex indices (integers).
*/
void store_triangle(int *t)
{
	int nt = (++_MC_S->nT)&_MC_A;
	void *pt;
	if(!nt)//create a new _MC_S->T[*][][]
	{
		if(!((++_MCnT)&(incDimTVN - 1)))//expand the _MC_S->T[] dimension in incDimTVN
		{
			pt = realloc(_MC_S->T,(_MCnT + incDimTVN)*sizeof(void*));
			if(pt)
				_MC_S->T = (int (**)[3])pt;
			else
			{
				_MC_memory_error(1);
				return;
			}
		}
		_MC_S->T[_MCnT] = (int (*)[3])malloc(3*_MC_DVE*sizeof(int));
		if(!_MC_S->T[_MCnT])
		{
			_MC_memory_error(1);
			return;
		}
	}
	memcpy(_MC_S->T[_MCnT][nt],t,3*sizeof(int));
}

/******************************************************************
Auxiliary function that calculates the normal if a vertex
of the cube lies on the isosurface.
*/
int surfint(int x, int y, int z, float *r, float *n)
{
	r[0] = x; r[1] = y; r[2] = z;
	if(x == 0)
		n[0] = _MC_F[z][y][0] - _MC_F[z][y][1];
	else if(x == _MCnx)
		n[0] = _MC_F[z][y][x - 1] - _MC_F[z][y][x];
	else
		n[0] = 0.5f*(_MC_F[z][y][x - 1] - _MC_F[z][y][x + 1]);
	if(y == 0)
		n[1] = _MC_F[z][0][x] - _MC_F[z][1][x];
	else if(y == _MCny)
		n[1] = _MC_F[z][y - 1][x] - _MC_F[z][y][x];
	else
		n[1] = 0.5f*(_MC_F[z][y - 1][x] - _MC_F[z][y + 1][x]);
	if(z == 0)
		n[2] = _MC_F[0][y][x] - _MC_F[1][y][x];
	else if(z == _MCnz)
		n[2] = _MC_F[z - 1][y][x] - _MC_F[z][y][x];
	else
		n[2] = 0.5f*(_MC_F[z - 1][y][x] - _MC_F[z + 1][y][x]);
	return store_point_normal(r,n);
}

/******************************************************************
This function find the MC33 case (using the index i, and the face and interior
tests). The correct triangle pattern is obtained from the arrays contained in
the MC33_LookUpTable.h file. The necessary vertices (intersection points) are
also calculated here (using trilinear interpolation).
       _____2_____
     /|          /|
   11 |<-3     10 |
   /____6_____ /  1     z
  |   |       |   |     |
  |   |_____0_|___|     |____y
  7  /        5  /     /
  | 8         | 9     x
  |/____4_____|/

The temporary matrices: _Ox, _Oy, _Nx and _Ny, and vectors: _OL and _NL are filled
and used here.*/

void find_case(int x, int y, int z, int i, float *v)
{
	const unsigned short int *pcase = 0;
	float t, r[3], n[3];
	int k, m, c;
	int f[6];//for the face tests
	int p[13] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	m = i&0x80;
	c = Case_Index[m? i^0xFF: i];
//	static int count4bits[] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
//	h = count4bits[i&0x0F] + count4bits[i>>4];
//	m = (h&0xFC? (h == 4? (i&0x80? 1: 0): 0): 1);//count bit
	switch(c>>8)//find the MC33 case
	{
	case 1://********************************************
		if(c&0x0080) m ^= 0x80;
		pcase = Case_1 + (c&0x7F);
		break;
	case 2://********************************************
		if(c&0x0080) m ^= 0x80;
		pcase = Case_2 + (c&0x7F);
		break;
	case 3://********************************************
		if(c&0x0080) m ^= 0x80;
		if((m? i: i^0xFF)&face_test1((c&0x7F)>>1,v))
			pcase = Case_3_2 + 4*(c&0x7F);
		else
			pcase = Case_3_1 + 2*(c&0x7F);
		break;
	case 4://********************************************
		if(c&0x0080) m ^= 0x80;
		if(interior_test((c&0x7F),0,v))
			pcase = Case_4_2 + 6*(c&0x7F);
		else
			pcase = Case_4_1 + 2*(c&0x7F);
		break;
	case 5://********************************************
		if(c&0x0080) m ^= 0x80;
		pcase = Case_5 + (c&0x7F);
		break;
	case 6://********************************************
		if(c&0x0080) m ^= 0x80;
		if((m? i: i^0xFF)&face_test1((c&0x7F)%6,v))
			pcase = Case_6_2 + 5*(c&0x7F);
		else if(interior_test((c&0x7F)/6,0,v))
			pcase = Case_6_1_2 + 7*(c&0x7F);
		else
			pcase = Case_6_1_1 + 3*(c&0x7F);
		break;
	case 7://********************************************
		if(c&0x0080) m ^= 0x80;
		switch(face_tests(f,i,(m? 1: -1),v))
		{
		case -3:
			pcase = Case_7_1 + 3*(c&0x7F);
			break;
		case -1:
			if(f[4] + f[5] == 1)
				pcase = Case_7_2_1 + 5*(c&0x7F);
			else
				pcase = (f[(33825>>((c&0x7F)<<1))&3] == 1? Case_7_2_3: Case_7_2_2) + 5*(c&0x7F);
			break;
		case 1:
			if(f[4] + f[5] == -1)
				pcase = Case_7_3_3 + 9*(c&0x7F);
			else
				pcase = (f[(33825>>((c&0x7F)<<1))&3] == 1? Case_7_3_2: Case_7_3_1) + 9*(c&0x7F);
			break;
		case 3:
			if(interior_test((c&0x7F)>>1,0,v))
				pcase = Case_7_4_2 + 9*(c&0x7F);
			else
				pcase = Case_7_4_1 + 5*(c&0x7F);
		}
		break;
	case 8://********************************************
		pcase = Case_8 + (c&0x7F);
		break;
	case 9://********************************************
		pcase = Case_9 + (c&0x7F);
		break;
	case 10://********************************************
		switch(face_tests(f,i,(m? 1: -1),v))
		{
		case -2:
			if(c&0x7F? interior_test(0,0,v)||interior_test(c&0x01? 1: 3,0,v): interior_test(0,0,v))
				pcase = Case_10_1_2_1 + 8*(c&0x7F);
			else
				pcase = Case_10_1_1_1 + 4*(c&0x7F);
			break;
		case 2:
			if(c&0x7F? interior_test(2,0,v)||interior_test(c&0x01? 3: 1,0,v): interior_test(1,0,v))
				pcase = Case_10_1_2_2 + 8*(c&0x7F);
			else
				pcase = Case_10_1_1_2 + 4*(c&0x7F);
			break;
		case 0:
			pcase = (f[4>>((c&0x7F)<<1)] == 1? Case_10_2_2: Case_10_2_1) + 8*(c&0x7F);
		}
		break;
	case 11://********************************************
		pcase = Case_11 + (c&0x7F);
		break;
	case 12://********************************************
		switch(face_tests(f,i,(m? 1: -1),v))
		{
		case -2:
			if(interior_test((int)_12_test_index[0][(int)(c&0x7F)],0,v))
				pcase = Case_12_1_2_1 + 8*(c&0x7F);
			else
				pcase = Case_12_1_1_1 + 4*(c&0x7F);
			break;
		case 2:
			if(interior_test((int)_12_test_index[1][(int)(c&0x7F)],0,v))
				pcase = Case_12_1_2_2 + 8*(c&0x7F);
			else
				pcase = Case_12_1_1_2 + 4*(c&0x7F);
			break;
		case 0:
			pcase = (f[(int)_12_test_index[2][(int)(c&0x7F)]] == 1? Case_12_2_2: Case_12_2_1) + 8*(c&0x7F);
		}
		break;
	case 13://********************************************
		c = face_tests(f,i,(m? 1: -1),v);
		switch(abs(c))
		{
		case 6:
			pcase = Case_13_1 + 4*(c > 0);
			break;
		case 4:
			c >>= 2;
			i = 0;
			while(f[i] != -c)
				++i;
			pcase = Case_13_2 + 6*(3*c + 3 + i);
			i = 1;
			break;
		case 2:
			c = ((((((((f[0] < 0)<<1)|(f[1] < 0))<<1)|(f[2] < 0))<<1)|
					(f[3] < 0))<<1)|(f[4] < 0);
			pcase = Case_13_3 + 10*(25 - c + (((c > 10) + (c > 20))<<1));
			break;
		case 0:
			c = ((f[1] < 0)<<1)|(f[5] < 0);
			if(f[0]*f[1]*f[5] == 1)
				pcase = Case_13_4 + 12*c;
			else
			{
				i = interior_test(c,1,v);
				if(i)
					pcase = Case_13_5_2 + 10*(c|((i&1)<<2));
				else
				{
					pcase = Case_13_5_1 + 6*c;
					i = 1;
				}
			}
		}
		break;
	case 14://********************************************
		pcase = Case_14 + (c&0x7F);
		break;
	}
	while(i)
	{
		i = *(pcase++);
		for(k=0;k<3;++k)
		{
			c = i&0x0F;
			if(p[c] < 0)
			{
				switch(c)//the vertices r[3] and normals n[3] are calculated here
				{
				case 0:
					if(z || x)
						p[0] = _Oy[y][x];
					else
					{
						if(v[0] == 0.0f)
						{
							p[0] = surfint(0,y,0,r,n);
							if(signbf(v[3])) p[3] = p[0];
							if(signbf(v[4])) p[8] = p[0];
						}
						else if(v[1] == 0.0f)
						{
							p[0] = surfint(0,y + 1,0,r,n);
							if(signbf(v[2])) _NL[0] = p[1] = p[0];
							if(signbf(v[5])) _Ox[y + 1][0] = p[9] = p[0];
						}
						else
						{
							t = v[0]/(v[0] - v[1]);
							r[0] = r[2] = 0.0f;
							r[1] = y + t;
							n[0] = (v[4] - v[0])*(1.0f - t) + (v[5] - v[1])*t;
							n[1] = v[1] - v[0];
							n[2] = (v[3] - v[0])*(1.0f - t) + (v[2] - v[1])*t;
							p[0] = store_point_normal(r,n);
						}
					}
					break;
				case 1:
					if(x)
						p[1] = _NL[x];
					else
					{
						if(v[1] == 0.0f)
						{
							_NL[0] = p[1] = surfint(0,y + 1,z,r,n);
							if(signbf(v[0])) p[0] = p[1];
							if(signbf(v[5]))
							{
								p[9] = p[1];
								if(z == 0) _Ox[y + 1][0] = p[9];
							}
						}
						else if(v[2] == 0.0f)
						{
							_NL[0] = p[1] = surfint(0,y + 1,z + 1,r,n);
							if(signbf(v[3])) _Ny[y][0] = p[2] = p[1];
							if(signbf(v[6])) _Nx[y + 1][0] = p[10] = p[1];
						}
						else
						{
							t = v[1]/(v[1] - v[2]);
							r[0] = 0.0f; r[1] = y + 1;
							r[2] = z + t;
							n[0] = (v[5] - v[1])*(1.0f - t) + (v[6] - v[2])*t;
							n[1] = (y + 1 < _MCny? 0.5f*((_MC_F[z][y][0] - _MC_F[z][y + 2][0])*(1.0f - t)
										+ (_MC_F[z + 1][y][0] - _MC_F[z + 1][y + 2][0])*t):
										(v[1] - v[0])*(1.0f - t) + (v[2] - v[3])*t);
							n[2] = v[2] - v[1];
							_NL[0] = p[1] = store_point_normal(r,n);
						}
					}
					break;
				case 2:
					if(x)
						p[2] = _Ny[y][x];
					else
					{
						if(v[3] == 0.0f)
						{
							_Ny[y][0] = p[2] = surfint(0,y,z + 1,r,n);
							if(signbf(v[0])) p[3] = p[2];
							if(signbf(v[7]))
							{
								p[11] = p[2];
								if(y == 0) _Nx[0][0] = p[11];
							}
						}
						else if(v[2] == 0.0f)
						{
							_Ny[y][0] = p[2] = surfint(0,y + 1,z + 1,r,n);
							if(signbf(v[1])) _NL[0] = p[1] = p[2];
							if(signbf(v[6])) _Nx[y + 1][0] = p[10] = p[2];
						}
						else
						{
							t = v[3]/(v[3] - v[2]);
							r[0] = 0.0f; r[2] = z + 1;
							r[1] = y + t;
							n[0] = (v[7] - v[3])*(1.0f - t) + (v[6] - v[2])*t;
							n[1] = v[2] - v[3];
							n[2] = (z + 1 < _MCnz? 0.5f*((_MC_F[z][y][0] - _MC_F[z + 2][y][0])*(1.0f - t)
										+ (_MC_F[z][y + 1][0] - _MC_F[z + 2][y + 1][0])*t):
										(v[3] - v[0])*(1.0f - t) + (v[2] - v[1])*t);
							_Ny[y][0] = p[2] = store_point_normal(r,n);
						}
					}
					break;
				case 3:
					if(y || x)
						p[3] = _OL[x];
					else
					{
						if(v[0] == 0.0f)
						{
							p[3] = surfint(0,0,z,r,n);
							if(signbf(v[1])) p[0] = p[3];
							if(signbf(v[4])) p[8] = p[3];
						}
						else if(v[3] == 0.0f)
						{
							p[3] = surfint(0,0,z + 1,r,n);
							if(signbf(v[2])) _Ny[0][0] = p[2] = p[3];
							if(signbf(v[7])) _Nx[0][0] = p[11] = p[3];
						}
						else
						{
							t = v[0]/(v[0] - v[3]);
							r[0] = r[1] = 0.0f;
							r[2] = z + t;
							n[0] = (v[4] - v[0])*(1.0f - t) + (v[7] - v[3])*t;
							n[1] = (v[1] - v[0])*(1.0f - t) + (v[2] - v[3])*t;
							n[2] = v[3] - v[0];
							p[3] = store_point_normal(r,n);
						}
					}
					break;
				case 4:
					if(z)
						p[4] = _Oy[y][x + 1];
					else
					{
						if(v[4] == 0.0f)
						{
							_Oy[y][x + 1] = p[4] = surfint(x + 1,y,0,r,n);
							if(signbf(v[7])) p[7] = p[4];
							if(signbf(v[0])) p[8] = p[4];
							if(y == 0)
								_OL[x + 1] = p[7];
						}
						else if(v[5] == 0.0f)
						{
							_Oy[y][x + 1] = p[4] = surfint(x + 1,y + 1,0,r,n);
							if(signbf(v[6])) _NL[x + 1] = p[5] = p[4];
							if(signbf(v[1])) _Ox[y + 1][x] = p[9] = p[4];
						}
						else
						{
							t = v[4]/(v[4] - v[5]);
							r[0] = x + 1; r[2] = 0.0f;
							r[1] = y + t;
							n[0] = (x + 1 < _MCnx? 0.5f*((_MC_F[0][y][x] - _MC_F[0][y][x + 2])*(1.0f - t)
										+ (_MC_F[0][y + 1][x] - _MC_F[0][y + 1][x + 2])*t):
										(v[4] - v[0])*(1.0f - t) + (v[5] - v[1])*t);
							n[1] = v[5] - v[4];
							n[2] = (v[7] - v[4])*(1.0f - t) + (v[6] - v[5])*t;
							_Oy[y][x + 1] = p[4] = store_point_normal(r,n);
						}
					}
					break;
				case 5:
					if(v[5] == 0.0f)
					{
						if(signbf(v[4]))
						{
							if(z)
							{
								_NL[x + 1] = p[5] = p[4] = _Oy[y][x + 1];
								if(signbf(v[1])) p[9] = p[5];
							}
							else
							{
								_NL[x + 1] = p[5] = _Oy[y][x + 1] = p[4] = surfint(x + 1,y + 1,0,r,n);
								if(signbf(v[1])) _Ox[y + 1][x] = p[9] = p[5];
							}
						}
						else if(signbf(v[1]))
						{
							if(z)
								_NL[x + 1] = p[5] = p[9] = _Ox[y + 1][x];
							else
								_NL[x + 1] = p[5] = _Ox[y + 1][x] = p[9] = surfint(x + 1,y + 1,0,r,n);
						}
						else
							_NL[x + 1] = p[5] = surfint(x + 1,y + 1,z,r,n);
					}
					else if(v[6] == 0.0f)
					{
						_NL[x + 1] = p[5] = surfint(x + 1,y + 1,z + 1,r,n);
						if(signbf(v[2])) _Nx[y + 1][x] = p[10] = p[5];
						if(signbf(v[7])) _Ny[y][x + 1] = p[6] = p[5];
					}
					else
					{
						t = v[5]/(v[5] - v[6]);
						r[0] = x + 1; r[1] = y + 1;
						r[2] = z + t;
						n[0] = (x + 1 < _MCnx? 0.5f*((_MC_F[z][y + 1][x] - _MC_F[z][y + 1][x + 2])*(1.0f - t)
									+ (_MC_F[z + 1][y + 1][x] - _MC_F[z + 1][y + 1][x + 2])*t):
									(v[5] - v[1])*(1.0f - t) + (v[6] - v[2])*t);
						n[1] = (y + 1 < _MCny? 0.5f*((_MC_F[z][y][x + 1] - _MC_F[z][y + 2][x + 1])*(1.0f - t)
									+ (_MC_F[z + 1][y][x + 1] - _MC_F[z + 1][y + 2][x + 1])*t):
									(v[5] - v[4])*(1.0f - t) + (v[6] - v[7])*t);
						n[2] = v[6] - v[5];
						_NL[x + 1] = p[5] = store_point_normal(r,n);
					}
					break;
				case 6:
					if(v[7] == 0.0f)
					{
						if(signbf(v[3]))
						{
							if(y)
							{
								_Ny[y][x + 1] = p[6] = p[11] = _Nx[y][x];
								if(signbf(v[4])) p[7] = p[6];
							}
							else
							{
								_Ny[y][x + 1] = p[6] = _Nx[0][x] = p[11] = surfint(x + 1,0,z + 1,r,n);
								if(signbf(v[4])) _OL[x + 1] = p[7] = p[6];
							}
						}
						else if(signbf(v[4]))
						{
							if(y)
								_Ny[y][x + 1] = p[6] = p[7] = _OL[x + 1];
							else
								_Ny[y][x + 1] = p[6] = _OL[x + 1] = p[7] = surfint(x + 1,0,z + 1,r,n);
						}
						else
							_Ny[y][x + 1] = p[6] = surfint(x + 1,y,z + 1,r,n);
					}
					else if(v[6] == 0.0f)
					{
						_Ny[y][x + 1] = p[6] = surfint(x + 1,y + 1,z + 1,r,n);
						if(signbf(v[5])) _NL[x + 1] = p[5] = p[6];
						if(signbf(v[2])) _Nx[y + 1][x] = p[10] = p[6];
					}
					else
					{
						t = v[7]/(v[7] - v[6]);
						r[0] = x + 1; r[2] = z + 1;
						r[1] = y + t;
						n[0] = (x + 1 < _MCnx? 0.5f*((_MC_F[z + 1][y][x] - _MC_F[z + 1][y][x + 2])*(1.0f - t)
									+ (_MC_F[z + 1][y + 1][x] - _MC_F[z + 1][y + 1][x + 2])*t):
									(v[7] - v[3])*(1.0f - t) + (v[6] - v[2])*t);
						n[1] = v[6] - v[7];
						n[2] = (z + 1 < _MCnz? 0.5f*((_MC_F[z][y][x + 1] - _MC_F[z + 2][y][x + 1])*(1.0f - t)
										+ (_MC_F[z][y + 1][x + 1] - _MC_F[z + 2][y + 1][x + 1])*t):
									(v[7] - v[4])*(1.0f - t) + (v[6] - v[5])*t);
						_Ny[y][x + 1] = p[6] = store_point_normal(r,n);
					}
					break;
				case 7:
					if(y)
						p[7] = _OL[x + 1];
					else
					{
						if(v[4] == 0.0f)
						{
							_OL[x + 1] = p[7] = surfint(x + 1,0,z,r,n);
							if(signbf(v[0])) p[8] = p[7];
							if(signbf(v[5]))
							{
								p[4] = p[7];
								if(z == 0)
									_Oy[0][x + 1] = p[4];
							}
						}
						else if(v[7] == 0.0f)
						{
							_OL[x + 1] = p[7] = surfint(x + 1,0,z + 1,r,n);
							if(signbf(v[6])) _Ny[0][x + 1] = p[6] = p[7];
							if(signbf(v[3])) _Nx[0][x] = p[11] = p[7];
						}
						else
						{
							t = v[4]/(v[4] - v[7]);
							r[0] = x + 1; r[1] = 0.0f;
							r[2] = z + t;
							n[0] = (x + 1 < _MCnx? 0.5f*((_MC_F[z][0][x] - _MC_F[z][0][x + 2])*(1.0f - t)
										+ (_MC_F[z + 1][0][x] - _MC_F[z + 1][0][x + 2])*t):
										(v[4] - v[0])*(1.0f - t) + (v[7] - v[3])*t);
							n[1] = (v[5] - v[4])*(1.0f - t) + (v[6] - v[7])*t;
							n[2] = v[7] - v[4];
							_OL[x + 1] = p[7] = store_point_normal(r,n);
						}
					}
					break;
				case 8:
					if(z || y)
						p[8] = _Ox[y][x];
					else
					{
						if(v[0] == 0.0f)
						{
							p[8] = surfint(x,0,0,r,n);
							if(signbf(v[1])) p[0] = p[8];
							if(signbf(v[3])) p[3] = p[8];
						}
						else if(v[4] == 0.0f)
						{
							p[8] = surfint(x + 1,0,0,r,n);
							if(signbf(v[5]))
								_Oy[0][x + 1] = p[4] = p[8];
							if(signbf(v[7]))
								_OL[x + 1] = p[7] = p[8];
						}
						else
						{
							t = v[0]/(v[0] - v[4]);
							r[1] = r[2] = 0.0f;
							r[0] = x + t;
							n[0] = v[4] - v[0];
							n[1] = (v[1] - v[0])*(1.0f - t) + (v[5] - v[4])*t;
							n[2] = (v[3] - v[0])*(1.0f - t) + (v[7] - v[4])*t;
							p[8] = store_point_normal(r,n);
						}
					}
					break;
				case 9:
					if(z)
						p[9] = _Ox[y + 1][x];
					else
					{
						if(v[1] == 0.0f)
						{
							_Ox[y + 1][x] = p[9] = surfint(x,y + 1,0,r,n);
							if(signbf(v[2]))
							{
								p[1] = p[9];
								if(x == 0) _NL[0] = p[1];
							}
							if(signbf(v[0])) p[0] = p[9];
						}
						else if(v[5] == 0.0f)
						{
							_Ox[y + 1][x] = p[9] = surfint(x + 1,y + 1,0,r,n);
							if(signbf(v[6])) _NL[x + 1] = p[5] = p[9];
							if(signbf(v[4])) _Oy[y][x + 1] = p[4] = p[9];
						}
						else
						{
							t = v[1]/(v[1] - v[5]);
							r[1] = y + 1; r[2] = 0.0f;
							r[0] = x + t;
							n[0] = v[5] - v[1];
							n[1] = (y + 1 < _MCny? 0.5f*((_MC_F[0][y][x] - _MC_F[0][y + 2][x])*(1.0f - t)
										+ (_MC_F[0][y][x + 1] - _MC_F[0][y + 2][x + 1])*t):
										(v[1] - v[0])*(1.0f - t) + (v[5] - v[4])*t);
							n[2] = (v[2] - v[1])*(1.0f - t) + (v[6] - v[5])*t;
							_Ox[y + 1][x] = p[9] = store_point_normal(r,n);
						}
					}
					break;
				case 10:
					if(v[2] == 0.0f)
					{
						if(signbf(v[1]))
						{
							if(x)
							{
								_Nx[y + 1][x] = p[10] = p[1] = _NL[x];
								if(signbf(v[3])) p[2] = p[10];
							}
							else
							{
								_Nx[y + 1][0] = p[10] = _NL[0] = p[1] = surfint(0,y + 1,z + 1,r,n);
								if(signbf(v[3])) _Ny[y][0] = p[2] = p[10];
							}
						}
						else if(signbf(v[3]))
						{
							if(x)
								_Nx[y + 1][x] = p[10] = p[2] = _Ny[y][x];
							else
								_Nx[y + 1][0] = p[10] = _Ny[y][0] = p[2] = surfint(0,y + 1,z + 1,r,n);
						}
						else
							_Nx[y + 1][x] = p[10] = surfint(x,y + 1,z + 1,r,n);
					}
					else if(v[6] == 0.0f)
					{
						_Nx[y + 1][x] = p[10] = surfint(x + 1,y + 1,z + 1,r,n);
						if(signbf(v[5])) _NL[x + 1] = p[5] = p[10];
						if(signbf(v[7])) _Ny[y][x + 1] = p[6] = p[10];
					}
					else
					{
						t = v[2]/(v[2] - v[6]);
						r[1] = y + 1; r[2] = z + 1;
						r[0] = x + t;
						n[0] = v[6] - v[2];
						n[1] = (y + 1 < _MCny? 0.5f*((_MC_F[z + 1][y][x] - _MC_F[z + 1][y + 2][x])*(1.0f - t)
									+ (_MC_F[z + 1][y][x + 1] - _MC_F[z + 1][y + 2][x + 1])*t):
									(v[2] - v[3])*(1.0f - t) + (v[6] - v[7])*t);
						n[2] = (z + 1 < _MCnz? 0.5f*((_MC_F[z][y + 1][x] - _MC_F[z + 2][y + 1][x])*(1.0f - t)
									+ (_MC_F[z][y + 1][x + 1] - _MC_F[z + 2][y + 1][x + 1])*t):
									(v[2] - v[1])*(1.0f - t) + (v[6] - v[5])*t);
						_Nx[y + 1][x] = p[10] = store_point_normal(r,n);
					}
					break;
				case 11:
					if(y)
						p[11] = _Nx[y][x];
					else
					{
						if(v[3] == 0.0f)
						{
							_Nx[0][x] = p[11] = surfint(x,0,z + 1,r,n);
							if(signbf(v[0])) p[3] = p[11];
							if(signbf(v[2]))
							{
								p[2] = p[11];
								if(x == 0)
									_Ny[0][0] = p[2];
							}
						}
						else if(v[7] == 0.0f)
						{
							_Nx[0][x] = p[11] = surfint(x + 1,0,z + 1,r,n);
							if(signbf(v[4])) _OL[x + 1] = p[7] = p[11];
							if(signbf(v[6])) _Ny[0][x + 1] = p[6] = p[11];
						}
						else
						{
							t = v[3]/(v[3] - v[7]);
							r[1] = 0.0f; r[2] = z + 1;
							r[0] = x + t;
							n[0] = v[7] - v[3];
							n[1] = (v[2] - v[3])*(1.0f - t) + (v[6] - v[7])*t;
							n[2] = (z + 1 < _MCnz? 0.5f*((_MC_F[z][0][x] - _MC_F[z + 2][0][x])*(1.0f - t)
										+ (_MC_F[z][0][x + 1] - _MC_F[z + 2][0][x + 1])*t):
										(v[3] - v[0])*(1.0f - t) + (v[7] - v[4])*t);
							_Nx[0][x] = p[11] = store_point_normal(r,n);
						}
					}
				break;
				case 12:
					r[0] = x + 0.5f; r[1] = y + 0.5f; r[2] = z + 0.5f;
					n[0] = v[4] + v[5] + v[6] + v[7] - v[0] - v[1] - v[2] - v[3];
					n[1] = v[1] + v[2] + v[5] + v[6] - v[0] - v[3] - v[4] - v[7];
					n[2] = v[2] + v[3] + v[6] + v[7] - v[0] - v[1] - v[4] - v[5];
					p[12] = store_point_normal(r,n);
					break;
				}
			}
			f[k] = p[c];//now f contains the vertex indices of the triangle
			i >>= 4;
		}
		if(f[0] != f[1] && f[0] != f[2] && f[1] != f[2])//to avoid zero area triangles
		{
#ifndef MC_Normal_neg
			if(m)//The order of triangle vertices is reversed if m is not zero
#else //it is also reversed if MC_Normal_neg was defined
			if(!m)
#endif
				{f[2] = f[0]; f[0] = p[c];}
			store_triangle(f);
		}
	}
}

int write_bin_s(char *filename, surface *S)
{
	FILE *out;
	int i, nl = S->nV%S->dim2 + 1;
	int n = S->nT/S->dim2;
	out = fopen(filename,"wb");
	if(!out)
		return -1;
	fputs(".sup",out);

	fwrite(&S->ISO,sizeof(float),1,out);
	fwrite(&S->nV,sizeof(int),1,out);
	fwrite(&S->nT,sizeof(int),1,out);

	for(i=0;i<n;++i)
		fwrite(S->T[i],3*S->dim2*sizeof(int),1,out);
	fwrite(S->T[n],3*(S->nT%S->dim2 + 1)*sizeof(int),1,out);

	n = S->nV/S->dim2;
	for(i=0;i<n;++i)
		fwrite(S->V[i],3*S->dim2*sizeof(float),1,out);
	fwrite(S->V[n],3*nl*sizeof(float),1,out);
	for(i=0;i<n;++i)
		fwrite(S->N[i],3*S->dim2*sizeof(float),1,out);
	fwrite(S->N[n],3*nl*sizeof(float),1,out);
	for(i=0;i<n;++i)
		fwrite(S->color[i],S->dim2*sizeof(int),1,out);
	i = fwrite(S->color[n],nl*sizeof(int),1,out);
	fclose(out);
	return -(i == 1);
}

void free_surface_memory(surface **S)
{
	int n, m;
	if(!(*S))
		return;
	if((*S)->T)
	{
		n = (*S)->nT/(*S)->dim2 + 1;
		while(n)
			free((*S)->T[--n]);
		free((*S)->T);
	}
	n = m = (*S)->nV/(*S)->dim2 + 1;
	if((*S)->V)
	{
		while(n)
			free((*S)->V[--n]);
		free((*S)->V);
	}
	if((*S)->N)
	{
		n = m;
		while(n)
			free((*S)->N[--n]);
		free((*S)->N);
	}
	if((*S)->color)
	{
		n = m;
		while(n)
			free((*S)->color[--n]);
		free((*S)->color);
	}
	free(*S);
	*S = 0;
}

surface* read_bin_s(char *filename)
{
	FILE *in;
	surface *S;
	int i, n, nl, k = 0;

	in = fopen(filename,"rb");
	if(!in)
		return 0;
	fread(&n,sizeof(int),1,in);
	S = (surface*)malloc(sizeof(surface));
	if(n != 0x7075732e || !S)
	{
		fclose(in);
		free(S);
		return 0;
	}
	S->dim2 = _MC_DVE;
	fread(&S->ISO,sizeof(float),1,in);
	fread(&S->nV,sizeof(int),1,in);
	fread(&S->nT,sizeof(int),1,in);
	n = S->nT>>_MC_N;
	nl = (S->nT&_MC_A) + 1;

	S->T = (int (**)[3])malloc((n + 1)*sizeof(void*));
	if(!S->T)
	{
		free(S);
		fclose(in);
		return 0;
	}
	S->T[n] = (int (*)[3])malloc(3*nl*sizeof(int));
	i = n;
	while(i)
		S->T[--i] = (int (*)[3])malloc(3*_MC_DVE*sizeof(int));
	if(S->T[0])
	{
		for(i=0;i<n;++i)
			if(!fread(S->T[i],3*_MC_DVE*sizeof(int),1,in))
				k = 1;
		if(!fread(S->T[i],3*nl*sizeof(int),1,in))
			k = 1;
	}
	else
		k = 1;

	n = S->nV>>_MC_N;
	nl = (S->nV&_MC_A) + 1;
	S->V = (float (**)[3])malloc((n + 1)*sizeof(void*));
	S->N = (float (**)[3])malloc((n + 1)*sizeof(void*));
	S->color = (int **)malloc((n + 1)*sizeof(void*));
	if(!(S->V && S->N && S->color))
		k = 1;
	if(k)
	{
		for(i=0;i<n+1;++i)
			free(S->T[i]);
		free(S->T);
		free(S->V);
		free(S->N);
		free(S->color);
		free(S);
		fclose(in);
		return 0;
	}
	S->V[n] = (float (*)[3])malloc(3*nl*sizeof(float));
	S->N[n] = (float (*)[3])malloc(3*nl*sizeof(float));
	S->color[n] = (int *)malloc(nl*sizeof(int));
	i = n;
	while(i)
	{
		S->V[--i] = (float (*)[3])malloc(3*_MC_DVE*sizeof(float));
		S->N[i] = (float (*)[3])malloc(3*_MC_DVE*sizeof(float));
		S->color[i] = (int *)malloc(_MC_DVE*sizeof(int));
	}
	if(S->V[i] && S->N[i] && S->color[i])
	{
		for(i=0;i<n;++i)
			if(!fread(S->V[i],3*_MC_DVE*sizeof(float),1,in))
				k = 1;
		if(!fread(S->V[n],3*nl*sizeof(float),1,in))
			k = 1;
		for(i=0;i<n;++i)
			if(!fread(S->N[i],3*_MC_DVE*sizeof(float),1,in))
				k = 1;
		if(!fread(S->N[n],3*nl*sizeof(float),1,in))
			k = 1;
		for(i=0;i<n;++i)
			if(!fread(S->color[i],_MC_DVE*sizeof(int),1,in))
				k = 1;
		if(!fread(S->color[n],nl*sizeof(int),1,in))
			k = 1;
	}
	else k = 1;
	if(k)
		free_surface_memory(&S);
	fclose(in);
	return S;
}

int write_txt_s(char* filename, surface *S)
{
	FILE *out;
	int i, *t, nt = S->nT + 1;
	int nv = S->nV + 1;
	float *r;

	out = fopen(filename,"w");
	if(!out)
		return -1;

	fprintf(out,"isovalue: %10.5E\n\nVERTICES:\n",S->ISO);
	fprintf(out,"%d\n\n",nv);
	for(i=0;i<nv;++i)
	{
		r = S->V[i/S->dim2][i%S->dim2];
		fprintf(out,"%9.6f %9.6f %9.6f\n",r[0],r[1],r[2]);
	}

	fprintf(out,"\n\nTRIANGLES:\n");
	fprintf(out,"%d\n\n",nt);
	for(i=0;i<nt;++i)
	{
		t = S->T[i/S->dim2][i%S->dim2];
		fprintf(out,"%8d %8d %8d\n",t[0],t[1],t[2]);
	}

	fprintf(out,"\n\nNORMALS:\n");
	for(i=0;i<nv;++i)
	{
		r = S->N[i/S->dim2][i%S->dim2];
		fprintf(out,"%9.6f %9.6f %9.6f\n",r[0],r[1],r[2]);
	}

	fprintf(out,"\n\nCOLORS:\n");
	for(i=0;i<nv;++i)
		fprintf(out,"%d\n",S->color[i/S->dim2][i%S->dim2]);
	i = fprintf(out,"\nEND\n");
	fclose(out);
	return -(i < 0);
}

void _MC_free_temp_O_N()
{
	free(_Ox); free(_Nx); free(_Oy); free(_Ny);
	free(_OL); free(_NL);
}

void clear_temp_isosurface()
{
	int y;
	for(y=0;y!=_MCny;++y)
	{
		free(_Ox[y]); free(_Nx[y]); free(_Oy[y]); free(_Ny[y]);
	}
	free(_Ox[_MCny]); free(_Nx[_MCny]);
	_MC_free_temp_O_N();
	_MC_F = 0;
}

int init_temp_isosurface(_GRD* grd)
{
	int x, y;
	if(!grd)
		return -1;
	_MCnx = grd->N[0];
	_MCny = grd->N[1];
	_MCnz = grd->N[2];
	for(y=0;y<3;++y)
	{
		_MC_O[y] = grd->r0[y];
		_MC_D[y] = grd->d[y];
	}
#ifndef _ORTHO_GRD
	_MC_NORT = grd->nonortho;
	if(_MC_NORT)
	{
		_MC__A = grd->_A;
		_MC_A_ = grd->A_;
		mult_Abf(grd->A_,_MC_O,_MC_O,0);
	}
#endif
	x = _MCnx*sizeof(int);
	_OL = (int*)malloc(x + sizeof(int));//edges 3 (only read) and 7
	_NL = (int*)malloc(x + sizeof(int));//edges 1 and 5 (only write)
	_Oy = (int**)malloc(_MCny*sizeof(int*));//edges 0 (only read) and 4
	_Ny = (int**)malloc(_MCny*sizeof(int*));//edges 2 and 6 (only write)
	_Ox = (int**)malloc((_MCny + 1)*sizeof(int*));//edges 8 (only read) and 9
	_Nx = (int**)malloc((_MCny + 1)*sizeof(int*));//edges 10 (only write) and 11
	if(!_Nx || !_NL)
	{
		_MC_free_temp_O_N();
		return -1;
	}
	_MC_F = grd->F;
	for(y=0;y<_MCny;++y)
	{
		_Ox[y] = (int*)malloc(x);
		_Nx[y] = (int*)malloc(x);
		_Oy[y] = (int*)malloc(x + sizeof(int));
		_Ny[y] = (int*)malloc(x + sizeof(int));
	}
	if(!_Ny[y - 1])
	{
		_Ox[y] = _Nx[y] = 0;
		clear_temp_isosurface();
		return -1;
	}
	_Ox[y] = (int*)malloc(x);
	_Nx[y] = (int*)malloc(x);
	if(!_Nx[y])
	{
		clear_temp_isosurface();
		return -1;
	}
	return 0;
}

surface* calc_isosurface(float iso)
{
	int x, y, z, nx = _MCnx, ny = _MCny, nz = _MCnz, i;
	GRD_data_type ***F = _MC_F, **F0, **F1, *V00, *V01, *V11, *V10;
	float V[12];
	float *v1 = V, *v2 = V + 4;
	if(!F)//The init_temp_isosurface function was not executed
		return 0;
	_MC_S = (surface*)malloc(sizeof(surface));
	if(!_MC_S)
		return 0;
	_MC_S->dim2 = _MC_DVE;
	_MCnT = _MCnV = _MC_S->nV = _MC_S->nT = -1;
	_MC_S->ISO = iso;
	_MC_S->T = 0;
	_MC_S->V = _MC_S->N = 0;
	_MC_S->color = 0;
	_MC_memoryfault = 0;
	for(z=0;z!=nz;++z)
	{
		F0 = *F;
		F1 = *(++F);
		for(y=0;y!=ny;++y)
		{
			V00 = *F0;
			V01 = *(++F0);
			V10 = *F1;
			V11 = *(++F1);
			v2[0] = iso - *V00;//the difference was inverted to use signbit function
			v2[1] = iso - *V01;
			v2[2] = iso - *V11;
			v2[3] = iso - *V10;
			//the eight least significant bits of i correspond to vertex indices. (x...x01234567)
			//If the bit is 1 then the vertex value is greater than zero.
			i = (((((signbf(v2[0])<<1)|signbf(v2[1]))<<1)|signbf(v2[2]))<<1)|signbf(v2[3]);
			for(x=0;x!=nx;++x)
			{
#ifdef __cplusplus
				std::swap(v1,v2);
#else
				{register float *P = v1; v1 = v2; v2 = P;}//v1 and v2 are exchanged
#endif
				v2[0] = iso - *(++V00);
				v2[1] = iso - *(++V01);
				v2[2] = iso - *(++V11);
				v2[3] = iso - *(++V10);
				i = ((((((((i&0x0F)<<1)|signbf(v2[0]))<<1)|signbf(v2[1]))<<1)|signbf(v2[2]))<<1)|signbf(v2[3]);
				if(i && i^0xFF)//i is different from 0 and 0xFF
				{
					if(v1 > v2) memcpy(v1 + 4,v2,4*sizeof(float));
					find_case(x,y,z,i,v1);
				}
			}
#ifdef __cplusplus
			std::swap(_OL,_NL);
#else
			{int* P = _OL; _OL = _NL; _NL = P;}//_OL and _NL are exchanged
#endif
		}
#ifdef __cplusplus
		std::swap(_Ox,_Nx);
		std::swap(_Oy,_Ny);
#else
		{int** P = _Ox; _Ox = _Nx; _Nx = P;}//_Ox and _Nx are exchanged
		{int** P = _Oy; _Oy = _Ny; _Ny = P;}//_Oy and _Ny are exchanged
#endif
	}
	if(_MC_memoryfault)
		free_surface_memory(&_MC_S);
	return _MC_S;
}

surface* calculate_isosurface(_GRD* grd, float iso)
{
	if(init_temp_isosurface(grd))
		return 0;
	calc_isosurface(iso);
	clear_temp_isosurface();
	return _MC_S;
}

#endif //marching_cubes_33_c

