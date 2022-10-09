/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *      Copyright (C) 2001-2011 by Computer Graphics Group, RWTH Aachen      *
 *                           www.openflipper.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenFlipper.                                        *
 *                                                                           *
 *  OpenFlipper is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenFlipper is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenFlipper. If not,                                  *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 14430 $                                                       *
 *   $Author: ebke $                                                      *
 *   $Date: 2012-04-23 12:26:49 +0200 (Mo, 23 Apr 2012) $                   *
 *                                                                           *
\*===========================================================================*/



//=============================================================================
//
//  CLASS Matrix4x4T
//
//=============================================================================


#ifndef ACG_MATRIX4X4_HH
#define ACG_MATRIX4X4_HH


//== INCLUDES =================================================================

#include "VectorT.hh"
#include <math.h>
#include <iostream>


//== NAMESPACES  ==============================================================


namespace ACG {


//== MACROS / HELPERS =========================================================

	      
#define matrixOperator(src, dst, op) { \
    Scalar *_a = dst;\
    const Scalar *_b = src;\
    *_a++ op *_b++; *_a++ op *_b++; *_a++ op *_b++; *_a++ op *_b++; \
    *_a++ op *_b++; *_a++ op *_b++; *_a++ op *_b++; *_a++ op *_b++; \
    *_a++ op *_b++; *_a++ op *_b++; *_a++ op *_b++; *_a++ op *_b++; \
    *_a++ op *_b++; *_a++ op *_b++; *_a++ op *_b++; *_a   op *_b  ; \
}
  
#define MAT(m,r,c) ((m)[(r)+((c)<<2)])

template<typename Scalar> inline bool checkEpsilon(Scalar x) {
  return fabs(x) < (1e-6);
}

template<> inline bool checkEpsilon(float x) {
  return fabs(x) < (1e-4);
}



//== CLASS DEFINITION =========================================================


/** Simple 4x4 Matrix. Inherited by GLMatrix.
*/

template <class Scalar>
class Matrix4x4T
{
public:
   
  /// constructor: uninitialized values
  Matrix4x4T() {}
 
  /// construct from other matrix type
  template <class OtherScalar>
  inline Matrix4x4T(const Matrix4x4T<OtherScalar> _rhs) {
    operator=(_rhs);
  }
  
  /** setup matrix using an array of N*N scalar values.
      elements are ordered 'column first' (like OpenGL) */
  inline Matrix4x4T(const Scalar _array[16]) {
      /*
       * Why aren't we doing a memcpy here?
       */
    matrixOperator(_array, mat_, =);
  }

  /// destructor
  ~Matrix4x4T() {}


  /// assignment from other matrix type
  template<typename otherScalar>
  inline Matrix4x4T<Scalar>& operator=(const Matrix4x4T<otherScalar>& _rhs) {
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
	operator()(i,j) = _rhs(i,j);
    return *this;
  }



  /// access operator (read and write)
  inline Scalar& operator()(unsigned int row, unsigned int col) {
    return MAT(mat_,row,col);
  }

  /// access operator (read only)
  inline const Scalar& operator()(unsigned int row, unsigned int col) const {
    return MAT(mat_,row,col);
  }


  /// compare two matrices (up to some epsilon)
  inline bool operator== (const Matrix4x4T<Scalar>& _rhs) const {
    int i;
    const Scalar *a = mat_;
    const Scalar *b = _rhs.mat_;
    for(i=0;i< 16;i++,a++,b++)
      if(! checkEpsilon( *a - *b ))
	return false;
    return true;
  }

  /// compare two matrices
  inline bool operator!= (const Matrix4x4T<Scalar>& _rhs) const {
    return !( operator==(_rhs) );
  }


  /// self + _rhs
  inline Matrix4x4T operator+ (const Matrix4x4T<Scalar>& _rhs) const {
    Matrix4x4T<Scalar> m(_rhs);
    matrixOperator(mat_, m.mat_, += );
    return m;
  }

  /// self - _rhs
  inline Matrix4x4T operator- (const Matrix4x4T<Scalar>& _rhs) const {
    Matrix4x4T<Scalar> m(*this);
    matrixOperator(_rhs.mat_, m.mat_, -= );
    return m;
  }

  /// self * _rhs
  Matrix4x4T operator*(const Matrix4x4T<Scalar>& inst) const;

  /// self * scalar
  Matrix4x4T operator*(const Scalar& scalar);


  /// self += _rhs
  inline Matrix4x4T& operator+= ( const Matrix4x4T<Scalar>& _rhs) {
    matrixOperator(_rhs.mat_, mat_, += );
    return *this;
  }

  /// self -= _rhs
  inline Matrix4x4T& operator-= ( const Matrix4x4T<Scalar>& _rhs) {
    matrixOperator(_rhs.mat_, mat_, -= );
    return *this;
  }

  /// self *= _rhs
  Matrix4x4T& operator*= (const Matrix4x4T<Scalar>& _rhs);

  /// multiply from left: self = _rhs * self
  Matrix4x4T& leftMult(const Matrix4x4T<Scalar>& _rhs);


  /// matrix by vector multiplication
  template <typename T>
  inline VectorT<T,4> operator*(const VectorT<T,4>& _v) const;

  /// transform point (x',y',z',1) = M * (x,y,z,1)
  template <typename T>
  inline VectorT<T,3> transform_point(const VectorT<T,3>& _v) const;

  /// transform vector (x',y',z',0) = A * (x,y,z,0)
  template <typename T>
  inline VectorT<T,3> transform_vector(const VectorT<T,3>& _v) const;

  /// sets all elements to zero
  inline void clear();

  /// setup an identity matrix
  inline void identity();


  /// check if the matrix is the identity ( up to an epsilon )
  inline bool is_identity() const {
    int i;
    const Scalar *a = mat_;
    Scalar b = 0.0;
    for(i=0;i< 16;i++,a++,b++) {
      if ( ( i == 0) || ( i == 5 ) || ( i == 10 ) || ( i == 15 ) )
         b = 1.0;
      else
         b = 0.0;
      if(! checkEpsilon( *a - b ))
        return false;
    }
    return true;
  }


  /// transpose matrix
  inline void transpose();

  
  /// matrix inversion (returns true on success)
  bool invert();

  Scalar determinant() const {
      return  mat_[12] * mat_[9] * mat_[6] * mat_[3] - mat_[8] * mat_[13] * mat_[6] * mat_[3] -
              mat_[12] * mat_[5] * mat_[10] * mat_[3] + mat_[4] * mat_[13] * mat_[10] * mat_[3] +
              mat_[8] * mat_[5] * mat_[14] * mat_[3] - mat_[4] * mat_[9] * mat_[14] * mat_[3] -
              mat_[12] * mat_[9] * mat_[2] * mat_[7] + mat_[8] * mat_[13] * mat_[2] * mat_[7] +
              mat_[12] * mat_[1] * mat_[10] * mat_[7] - mat_[0] * mat_[13] * mat_[10] * mat_[7] -
              mat_[8] * mat_[1] * mat_[14] * mat_[7] + mat_[0] * mat_[9] * mat_[14] * mat_[7] +
              mat_[12] * mat_[5] * mat_[2] * mat_[11] - mat_[4] * mat_[13] * mat_[2] * mat_[11] -
              mat_[12] * mat_[1] * mat_[6] * mat_[11] + mat_[0] * mat_[13] * mat_[6] * mat_[11] +
              mat_[4] * mat_[1] * mat_[14] * mat_[11] - mat_[0] * mat_[5] * mat_[14] * mat_[11] -
              mat_[8] * mat_[5] * mat_[2] * mat_[15] + mat_[4] * mat_[9] * mat_[2] * mat_[15] +
              mat_[8] * mat_[1] * mat_[6] * mat_[15] - mat_[0] * mat_[9] * mat_[6] * mat_[15] -
              mat_[4] * mat_[1] * mat_[10] * mat_[15] + mat_[0] * mat_[5] * mat_[10] * mat_[15];
  }


  /** access to data array. not very nice, but in case of 4x4 matrices
      this member can be used to pass matrices to OpenGL
      e.g. glLoadMatrixf(m.get_raw_data()); */
  inline const Scalar* get_raw_data() const { return mat_; }
  inline const Scalar* raw() const { return mat_; }
  inline const Scalar* data() const { return mat_; }

protected:

    Scalar mat_[16];
};


/// typedef
typedef Matrix4x4T<float>  Matrix4x4f;
/// typedef
typedef Matrix4x4T<double> Matrix4x4d;




//== IO to/from streams =======================================================


/// output matrix to ostream os
template<typename Scalar>
inline std::ostream& 
operator<<(std::ostream& os, const Matrix4x4T<Scalar>& m)
{
  for(int i=0; i<4; i++)
  {
    for(int j=0; j<4; j++)
      os << m(i,j) << " ";
    os << "\n";
  }
  return os;
}
 
 
/// read the space-separated components of a vector from a stream */
template<typename Scalar>
inline std::istream& 
operator>>(std::istream& is, Matrix4x4T<Scalar>& m) 
{
  for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
      is >> m(i,j);
  return is;
}


//=============================================================================

#undef matrixOperator
#undef MAT

//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_MATRIX4X4_C)
#define ACG_MATRIX4X4_TEMPLATES
#include "Matrix4x4T.cc"
#endif
//=============================================================================
#endif // ACG_MATRIX4X4_HH defined
//=============================================================================

