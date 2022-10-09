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
//  CLASS ImplicitSphere
//
//=============================================================================


#ifndef ISOEX_IMPLICITSPHERE_HH
#define ISOEX_IMPLICITSPHERE_HH


//== INCLUDES =================================================================

#include "Implicit.hh"
#include <IsoEx/Config/IsoExDefines.hh>

//== NAMESPACES ===============================================================

namespace IsoEx {

//== CLASS DEFINITION =========================================================

	      
/** \class ImplicitSphere ImplicitSphere.hh <IsoEx/Implicits/ImplicitSphere.hh>
    This class implements a very simple implicit object: a sphere given its
    center and its radius.
    \see IsoEx::Implicit
    \ingroup implicits
*/	      
template< class Vec3 >
class ISOEXDLLEXPORT ImplicitSphere : public Implicit<Vec3>
{
public:

  typedef typename Vec3::value_type real;
   
  /// \name Constructor & destructor
  //@{

  /// Constructor: given sphere center and radius
  ImplicitSphere(const Vec3& _center, real _radius)
    : center_(_center), 
      radius_(_radius), 
      sqr_radius_(_radius*_radius)
  {}

  /// Empty destructor
  ~ImplicitSphere() {}

  //@}



  /// \name Abstract interface of implicit objects, see also IsoEx::Implicit.
  //@{

  //_isovalue is not implemented in this class //YANG LIU

  bool is_inside(const Vec3& _point, const Real _isovalue = 0) const
  {
    return (center_ - _point).sqrnorm() <= sqr_radius_;
  }

  real scalar_distance(const Vec3& _point) const
  {
    return (center_ - _point).norm() - radius_;
  }

  bool directed_distance(const Vec3&  _p0,
			 const Vec3&  _p1,
			 Vec3&        _point,
			 Vec3&        _normal,
			 real&        _distance,
            const real _isovalue = 0) const
  {
    Vec3 orig(_p0), dir(_p1-_p0);

    double a = dir.sqrnorm();
    double b = 2.0*(dir | (orig - center_));
    double c = (orig - center_).sqrnorm() - radius_*radius_;
    double d = b*b - 4.0*a*c;

    if (d >= 0)
    {
      d = sqrt(d);

      double t1 = (-b-d) / (2.0*a);
      double t2 = (-b+d) / (2.0*a);
      double t  = 1.00001;
      if (t1 >= 0.0 && t1 < t) t = t1;
      if (t2 >= 0.0 && t2 < t) t = t2;

      if (t != 1.00001)
      {
	_point    = orig + dir*t;
	_normal   = (_point - center_) / radius_;
	_distance = ((dir | _normal) < 0.0) ? dir.norm()*t : -dir.norm()*t;
	return true;
      }
    }

    return false;
  }
  
  //@}


private:

  Vec3  center_;
  real            radius_;
  real            sqr_radius_;
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_IMPLICITSPHERE_HH defined
//=============================================================================

