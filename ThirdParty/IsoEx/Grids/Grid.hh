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
//  CLASS Grid
//
//=============================================================================


#ifndef ISOEX_GRIDBASE_HH
#define ISOEX_GRIDBASE_HH


//== INCLUDES =================================================================

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <vector>

#include <IsoEx/Config/IsoExDefines.hh>

//== NAMESPACES ===============================================================

namespace IsoEx {

//== CLASS DEFINITION =========================================================

	      
/** \class Grid Grid.hh <IsoEx/Grids/Grid.hh>
    
    This is the abstract base class for all grid objects.

    \ingroup grids
*/	      
template< class VectorType >
class ISOEXDLLEXPORT Grid
{
public:
   
  //-------------------------------------------------------------- public types

  typedef VectorType Vec3;
  typedef typename VectorType::value_type real;
  

  /// CubeIdx can be used to refer to cubes
  typedef unsigned int CubeIdx;

  /// The grid points can be referred to by PointIdx
  typedef unsigned int PointIdx;


  /** A CubeIterator is used to iterate over all cubes in the grid,
      from begin() to end(), just like for STL containers.  Its
      operator* returns a CubeIdx.
  */
  class CubeIterator
  {
  public:
    /** Constructor. 
	\internal
    */
    CubeIterator(unsigned int _idx) : idx_(_idx) {}
    /// Get cube index from iterator
    CubeIdx& operator*()  { return idx_; }
    /// Get cube index from iterator
    CubeIdx* operator->() { return &idx_; }
    /// Comparison
    bool operator==(const CubeIterator& _rhs) const { return idx_==_rhs.idx_;}
    /// Comparison
    bool operator!=(const CubeIterator& _rhs) const { return !(*this==_rhs); }
    /// Pre-increment
    CubeIterator& operator++() { ++idx_; return *this; }
  private:
    unsigned int idx_;
  };



  //------------------------------------------------------------ public methods


  /// Default constructor
  Grid() {}
  /// Destructor
  virtual ~Grid() {}



  /// \name Iterate over all cubes
  //@{

  /// Returns begin iterator for the grid's cubes
  CubeIterator begin() const { return CubeIterator(0); }

  /// Returns end iterator for the grid's cubes
  CubeIterator end()   const { return CubeIterator(n_cubes()); }

  //@}


  //------------------------------------------------ abstract virtual interface


  /// \name Abstract interface of grids
  //@{

  /// Number of cubes in the grid.
  virtual unsigned int n_cubes() const = 0;

  /// Number of cubes in the grid.
  virtual unsigned int n_points() const = 0;

  /// Return the PointIdx of the \b _corners'th corner of the cube \b _idx
  virtual PointIdx point_idx(CubeIdx _idx, unsigned char _corner) const = 0;

  /// Return the 3D point refered to by \b _idx.
  virtual Vec3  point(PointIdx _idx) const = 0;

  /// See IsoEx::Implicit::is_inside()
  virtual bool is_inside(PointIdx _pidx, const real _isovalue = 0) const = 0;

  /// See IsoEx::Implicit::scalar_distance()
  virtual real scalar_distance(PointIdx _pidx) const = 0;

  /// See IsoEx::Implicit::directed_distance()
  virtual bool directed_distance(const Vec3&  _p0,
				 const Vec3&  _p1,
				 Vec3&        _point,
				 Vec3&        _normal,
				 real&        _distance,
                const real _isovalue = 0
               ) const = 0;
  //@}
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_GRIDBASE_HH defined
//=============================================================================
