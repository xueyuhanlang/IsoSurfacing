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
//  CLASS Implicit
//
//=============================================================================


#ifndef ISOEX_IMPLICIT_HH
#define ISOEX_IMPLICIT_HH


//== INCLUDES =================================================================

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <IsoEx/Config/IsoExDefines.hh>
#include "IsoSurfGen/TinyVector.h"
#include <vector>

//== NAMESPACES ===============================================================

namespace IsoEx {

//== CLASS DEFINITION =========================================================

/** \class Implicit Implicit.hh <IsoEx/Implicits/Implicit.hh>
    
    This is the abstract base class for all objects representing implicit
    function. The three abstract virtual functions have to be overridden
    by derived classes.

    \ingroup implicits
*/
template< class Vec3 >
class ISOEXDLLEXPORT Implicit
{
public:
   
  typedef typename Vec3::value_type real;

  /// \name Constructor & destructor
  //@{

  /// constructor
  Implicit() {}
  /// destructor
  virtual ~Implicit() {}

  //@}



  /// \name Abstract interface of implicit objects
  //@{

  /** Is the point \b _point inside or outside w.r.t. the implicit
      function? The result should be the same as scalar_distance(_point)
      < 0, but it can be implemented more efficiently.
  */
  virtual bool is_inside(const Vec3& _point, const real _isovalue = 0) const = 0;


  /** Returns the <em>scalar distance</em> value of \b _point. Points
      inside the object should have negative distances, points outside
      positive distance.
  */
  virtual real scalar_distance(const Vec3& _point) const = 0;


  /** This method returns the <em>directed distance</em> of \b _p0 in
      the direction <b>_p1-_p0</b>. The resulting intersection point
      (casting a ray from _p0 in direction (_p1-_p0)) is stored in \b
      _point, the corresponding normal vector at this point is stored
      in \b _normal, the distance value is stored in _distance (again
      negative distance means \b _p0 is inside). Since the ray
      intersection may fail (because of precision problems) it is
      returned whether an intersection was found or not.
   */
  virtual bool directed_distance(const Vec3&  _p0,
				 const Vec3&  _p1,
				 Vec3&        _point,
				 Vec3&        _normal,
				 real&       _distance,
                 const real _isovalue = 0
      ) const = 0;


  virtual void directed_distance(
      const std::vector<TinyVector<real, 3>>& _p0vec,
      const std::vector<TinyVector<real, 3>>& _p1vec,
      const std::vector<real>& _p0vec_values,
      const std::vector<real>& _p1vec_values,
      std::vector<TinyVector<real, 3>>& _point,
      std::vector<TinyVector<real, 3>>& _normal,
      const real _isovalue = 0) const = 0;

  virtual void scalar_value(
      const std::vector<TinyVector<real, 3>>& pvec,
      std::vector<real>& vec_values) const = 0;
  //@}
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_IMPLICITSPHERE_HH defined
//=============================================================================

