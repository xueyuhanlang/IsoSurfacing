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
//  CLASS MatrixT and VectorT
//
//=============================================================================


#ifndef ISOEX_MATRIXT_HH
#define ISOEX_MATRIXT_HH


//== INCLUDES =================================================================

#include <vector>

//== NAMESPACES ===============================================================

namespace IsoEx {
namespace Math {

//== CLASS DEFINITION =========================================================


/** \class MatrixT MatrixT.hh <IsoEx/Math/MatrixT.hh>
    Simple matrix class whose scalar type is its template parameter.
*/	      
template <typename T>
class MatrixT
{
public:
   
  /// Construct with number of rows and columns.
  MatrixT(unsigned int _rows, unsigned int _cols) 
    : rows_(_rows), cols_(_cols)
  { data_.resize(_cols*_rows); }


  /// Read & write element access
  T& operator()(unsigned int _i, unsigned int _j)
  {
    assert (_i < rows_ && _j < cols_);
    return data_[_i * cols_ + _j];
  }

  /// Read only element access
  const T& operator()(unsigned int _i, unsigned int _j) const
  {
    assert (_i < rows_ && _j < cols_);
    return data_[_i * cols_ + _j];
  }
  
  /// Number of rows
  unsigned int rows() const { return rows_; }

  /// Number of columns
  unsigned int cols() const { return cols_; }
  
private:

  unsigned int    rows_, cols_;
  std::vector<T>  data_;
};


//== CLASS DEFINITION =========================================================


/** \class VectorT MatrixT.hh <IsoEx/Math/MatrixT.hh>
    Simple vector class whose scalar type is its template parameter.
*/	      
template <typename T> 
class VectorT
{
public:

  /// Construct with dimension
  VectorT(unsigned int _n) : n_(_n)
  { data_.resize(_n); }

  /// Read & write element access
  T& operator()(unsigned int _i)
  {
    assert(_i < n_);
    return data_[_i];
  }

  /// Read only element access
  const T& operator()(unsigned int _i) const
  {
    assert(_i < n_);
    return data_[_i];
  }

  /// Return vector's dimension
  unsigned int dim() const { return n_; }

private:
  unsigned int    n_;
  std::vector<T>  data_;
};


//=============================================================================
} // namespace Math
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_MATRIXT_HH defined
//=============================================================================

