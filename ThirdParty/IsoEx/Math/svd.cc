/*===========================================================================*\
 *                                                                           *
 *                                IsoEx                                      *
 *        Copyright (C) 2002 by Computer Graphics Group, RWTH Aachen         *
 *                         www.rwth-graphics.de                              *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/

//=============================================================================
//
//  SVD decomposition and back-substitution
//
//=============================================================================

#define SVD_C

//== INCLUDES =================================================================

#include "svd.hh"
#include <math.h>
#include <iostream>


//== NAMESPACES ===============================================================
namespace IsoEx {
namespace Math {
//=============================================================================

  
inline double dsqr(double a)
{
  double dsqrarg(a);
  return (dsqrarg == 0.0 ? 0.0 : dsqrarg*dsqrarg);
}


//-----------------------------------------------------------------------------


inline double dmax(double a, double b)
{
  double dmaxarg1(a), dmaxarg2(b);
  return (dmaxarg1 > dmaxarg2 ? dmaxarg1 : dmaxarg2);
}


//-----------------------------------------------------------------------------


inline int imin(int a, int b)
{
  int iminarg1(a), iminarg2(b);
  return (iminarg1 < iminarg2 ? iminarg1 : iminarg2);
}


//-----------------------------------------------------------------------------


inline double sign(double a, double b)
{
  return (b >= 0.0 ? fabs(a) : -fabs(a));
}


//-----------------------------------------------------------------------------


inline double dpythag(double a, double b)
{
  double absa(fabs(a)), absb(fabs(b));
  if (absa > absb) return absa*sqrt(1.0+dsqr(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+dsqr(absa/absb)));
}


//-----------------------------------------------------------------------------


template <typename MAT_MxN,
          typename VEC_N,
          typename MAT_NxN>
bool svd_decomp( MAT_MxN& A, VEC_N& S, MAT_NxN& V )
{
  int     m(A.rows()), n(A.cols());
  int     flag, i, its, j, jj, k, l(0), nm(0);
  double  anorm(0.0), c, f, g(0.0), h, s, scale(0.0), x, y, z;
  VEC_N   rv1(n);
  bool    convergence(true);
  
  
  
  for (i=0;i<n;i++)
  {
    l=i+1;
    rv1(i)=scale*g;
    g=s=scale=0.0;
    if (i < m)
    {
      for (k=i;k<m;k++) scale += fabs(A(k,i));
      if (scale)
      {
	for (k=i;k<m;k++)
	{
	  A(k,i) /= scale;
	  s += A(k,i)*A(k,i);
	}
	f=A(i,i);
	g = -sign(sqrt(s),f);
	h=f*g-s;
	A(i,i)=f-g;
	for (j=l;j<n;j++)
	{
	  for (s=0.0,k=i;k<m;k++) s += A(k,i)*A(k,j);
	  f=s/h;
	  for (k=i;k<m;k++) A(k,j) += f*A(k,i);
	}
	for (k=i;k<m;k++) A(k,i) *= scale;
      }
    }
    S(i)=scale *g;
    g=s=scale=0.0;
    if (i < m && i != (n - 1))
    {
      for (k=l;k<n;k++) scale += fabs(A(i,k));
      if (scale)
      {
	for (k=l;k<n;k++)
	{
	  A(i,k) /= scale;
	  s += A(i,k)*A(i,k);
	}
	f=A(i,l);
	g = -sign(sqrt(s),f);
	h=f*g-s;
	A(i,l)=f-g;
	for (k=l;k<n;k++) rv1(k)=A(i,k)/h;
	for (j=l;j<m;j++)
	{
	  for (s=0.0,k=l;k<n;k++) s += A(j,k)*A(i,k);
	  for (k=l;k<n;k++) A(j,k) += s*rv1(k);
	}
	for (k=l;k<n;k++) A(i,k) *= scale;
      }
    }
    anorm=dmax(anorm,(fabs(S(i))+fabs(rv1(i))));
  }

  for (i = (n-1); i >= 0;i--)
  {
    if (i < (n-1))
    {
      if (g)
      {
	for (j=l;j<n;j++) V(j,i)=(A(i,j)/A(i,l))/g;
	for (j=l;j<n;j++)
	{
	  for (s=0.0,k=l;k<n;k++) s += A(i,k)*V(k,j);
	  for (k=l;k<n;k++) V(k,j) += s*V(k,i);
	}
      }
      for (j=l;j<n;j++) V(i,j)=V(j,i)=0.0;
    }
    V(i,i)=1.0;
    g=rv1(i);
    l=i;
  }

  
  for (i=imin(m,n)-1;i>=0;i--)
  {
    l=i+1;
    g=S(i);
    for (j=l;j<n;j++) A(i,j)=0.0;
    if (g)
    {
      g=1.0/g;
      for (j=l;j<n;j++)
      {
	for (s=0.0,k=l;k<m;k++) s += A(k,i)*A(k,j);
	f=(s/A(i,i))*g;
	for (k=i;k<m;k++) A(k,j) += f*A(k,i);
      }
      for (j=i;j<m;j++) A(j,i) *= g;
    }
    else for (j=i;j<m;j++) A(j,i)=0.0;
    ++A(i,i);
  }

  
  for (k = (n-1); k>=0; k--)
  {
    for (its=1;its<=100;its++)
    {
      flag=1;
      for (l=k;l>=0;l--)
      {
	nm=l-1;
	if ((double)(fabs(rv1(l))+anorm) == anorm)
	{
	  flag=0;
	  break;
	}
	if ((double)(fabs(S(nm))+anorm) == anorm) break;
      }
      if (flag)
      {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++)
	{
	  f=s*rv1(i);
	  rv1(i)=c*rv1(i);
	  if ((double)(fabs(f)+anorm) == anorm) break;
	  g=S(i);
	  h=dpythag(f,g);
	  S(i)=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++)
	  {
	    y=A(j,nm);
	    z=A(j,i);
	    A(j,nm)=y*c+z*s;
	    A(j,i)=z*c-y*s;
	  }
	}
      }
      z=S(k);
      if (l == k)
      {
	if (z < 0.0)
	{
	  S(k) = -z;
	  for (j=0;j<n;j++) V(j,k) = -V(j,k);
	}
	break;
      }

      
      if (its == 30)
      {
	std::cerr << "SVD FAILURE (no convergence after 30 iters)\n";
	convergence = false;
      }
		
      x=S(l);
      nm=k-1;
      y=S(nm);
      g=rv1(nm);
      h=rv1(k);
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=dpythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
      c=s=1.0;

      for (j=l;j<=nm;j++)
      { 
	i=j+1;
	g=rv1(i);
	y=S(i);
	h=s*g;
	g=c*g;
	z=dpythag(f,h);
	rv1(j)=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++)
	{
	  x=V(jj,j);
	  z=V(jj,i);
	  V(jj,j)=x*c+z*s;
	  V(jj,i)=z*c-x*s;
	}
	z=dpythag(f,h);
	S(j)=z;
	if (z)
	{
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=0;jj<m;jj++)
	{
	  y=A(jj,j);
	  z=A(jj,i);
	  A(jj,j)=y*c+z*s;
	  A(jj,i)=z*c-y*s;
	}
      }
      rv1(l)=0.0;
      rv1(k)=f;
      S(k)=x;
    }
  }

  return convergence;
}


//----------------------------------------------------------------------------


template <typename MAT_MxN,
          typename VEC_N,
          typename MAT_NxN,
          typename VEC_M>
void
svd_backsub( const MAT_MxN& A, const VEC_M& S, const MAT_NxN& V,
	     const VEC_M& b, VEC_N& x )
{
  int     m(A.rows()), n(V.cols());
  int     jj, j, i;
  double  s;
  VEC_N   tmp(n);


  // Calculate U^T * B
  for ( j=0; j<n; j++ )
  {
    s = 0.0;    
    if (S(j))  // non-zero result only if S_j is nonzero
    { 
      for ( i=0; i<m; i++ ) s += A(i,j)*b(i);
      s /= S(j);
    }
    tmp(j) = s;
  }

  
  // Matrix multiply by V to get answer
  for ( j=0; j<n; j++ )
  { 
    s = 0.0;
    for ( jj=0; jj<n; jj++ ) s += V(j,jj)*tmp(jj);
    x(j) = s;
  }
}


//=============================================================================
} // namespace Math
} // namespace IsoEx
//=============================================================================

