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

#ifndef SVD_HH
#define SVD_HH


/** \file svd.hh
    This file provides singular value decompositon.
*/

//== NAMESPACES ===============================================================
namespace IsoEx {
namespace Math {
//=============================================================================
  
/**
   Computes the SVD of A into U*S*V^T. A will be destroyed!
   The diagonal matrix S is stored as a Nx1 vector.
   This is the implementation described in Numerical Recipies.
*/
template <typename MAT_MxN,
          typename VEC_M,
          typename MAT_NxN>
bool
svd_decomp( MAT_MxN& A, VEC_M& S, MAT_NxN& V );



/**
   SVD backsubstitution.
   This is the implementation described in Numerical Recipies.
*/
template <typename MAT_MxN,
          typename VEC_N,
          typename MAT_NxN,
          typename VEC_M>
void
svd_backsub( const MAT_MxN& A, const VEC_M& S, const MAT_NxN& V,
	     const VEC_M& b, VEC_N& x );



//=============================================================================
} // namespace Math
} // namespace IsoEx
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(SVD_C)
#define SVD_TEMPLATES
#include "svd.cc"
#endif
//=============================================================================
#endif // SVD_HH defined
//=============================================================================
