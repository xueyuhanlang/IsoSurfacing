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
//  CLASS RegularGrid
//
//=============================================================================


#ifndef ISOEX_REGULARGRID_HH
#define ISOEX_REGULARGRID_HH


//== INCLUDES =================================================================

#include "Grid.hh"
#include <IsoEx/Math/Matrix4x4T.hh>
#include <IsoEx/Config/IsoExDefines.hh>

//== NAMESPACES ===============================================================

namespace IsoEx
{

//== CLASS DEFINITION =========================================================


/** \class RegularGrid RegularGrid.hh <IsoEx/Grids/RegularGrid.hh>
    This class implements a regular 3D grid.
    \ingroup grids
*/
template <class Vec3 >
class ISOEXDLLEXPORT RegularGrid : public Grid< Vec3 >
{
public:

    typedef ACG::Matrix4x4T< typename Vec3::value_type >      Matrix;

    typedef typename Vec3::value_type real;

    /// CubeIdx can be used to refer to cubes
    typedef typename Grid< Vec3 >::CubeIdx CubeIdx;

    /// The grid points can be referred to by PointIdx
    typedef typename Grid< Vec3 >::PointIdx PointIdx;



    // ray plane intersections
    struct Intersection
    {
        real lambda;
        bool enter;
    };

    /** Constructor: given the implicit to be sampled, the grids extend
        in 3-space (origin and 3 axes) as well as the resolution (number
        of steps) of the axes. The grid will contain
        _x_res*_y_res*_z_res points and
        (_x_res-1)*(_y_res-1)*(_z_res-1) cubes.

        \note The resolution of each axis has to be less than 1024. This is
        to make sure that a cube or point can be represented by one
        integer.
    */
    RegularGrid( const Vec3&  _origin = Vec3( 0, 0, 0 ),
                 const Vec3&  _x_axis = Vec3( 1, 0, 0 ),
                 const Vec3&  _y_axis = Vec3( 0, 1, 0 ),
                 const Vec3&  _z_axis = Vec3( 0, 0, 1 ),
                 unsigned int            _x_res  = 1,
                 unsigned int            _y_res  = 1,
                 unsigned int            _z_res  = 1 )
    { initialize( _origin, _x_axis, _y_axis, _z_axis, _x_res, _y_res, _z_res ); }

    /// function to initialize the grid
    void initialize( const Vec3&  _origin,
                     const Vec3&  _x_axis,
                     const Vec3&  _y_axis,
                     const Vec3&  _z_axis,
                     unsigned int            _x_res,
                     unsigned int            _y_res,
                     unsigned int            _z_res );



    //------------------------------------------------------- mandatory interface

    /// Return number of cubes
    unsigned int n_cubes() const { return n_cubes_; }

    /// Return number of points
    unsigned int n_points() const { return n_points_; }

    /// Return the PointIdx of the \b _corners'th corner of the cube \b _idx
    PointIdx point_idx( CubeIdx _idx, unsigned char _corner ) const;

    /// Return the 3D point refered to by \b _idx.
    Vec3  point( PointIdx _idx ) const;

    /// Return the 3D point refered to by x,y,z
    Vec3  point( int x, int y, int z ) const;
    
    /// Return the nearest grid point
    PointIdx nearest_point( const Vec3 &_p );

    const Vec3& origin() const { return origin_; }
    const Vec3& x_axis() const { return x_axis_; }
    const Vec3& y_axis() const { return y_axis_; }
    const Vec3& z_axis() const { return z_axis_; }
    const Vec3& dx() const { return dx_; }
    const Vec3& dy() const { return dy_; }
    const Vec3& dz() const { return dz_; }

    unsigned int x_resolution() const { return x_res_; }
    unsigned int y_resolution() const { return y_res_; }
    unsigned int z_resolution() const { return z_res_; }

    /// transforms a point to local cube coordinates
    Vec3 to_local( const Vec3 &_pw );

    /// transforms a point for local cube coordinates to world coordinates
    Vec3 to_world( const Vec3 &_pl );

    /// function to intersect a Ray with the cube ( local coordinates )
    bool ray_intersect_local( const Vec3 &_o, const Vec3 &_d,
                              Vec3  &_entry,  Vec3 &_exit );

    ///  returns the volume of the grid
    real volume() { return this->x_axis().norm()*this->y_axis().norm()*this->z_axis().norm(); }

    /// returns the outer surface of the grid
    virtual real outer_surface()
    {
        return ( this->x_axis().norm()*this->y_axis().norm()*2.0f +
                 this->y_axis().norm()*this->z_axis().norm()*2.0f +
                 this->x_axis().norm()*this->z_axis().norm()*2.0f );
    }

protected:

    // matrices which transform points to local and world coordinates
    Matrix to_local_;
    Matrix to_world_;

    Vec3   	origin_, x_axis_, y_axis_, z_axis_, dx_, dy_, dz_;
    unsigned int x_res_, y_res_, z_res_, n_cubes_, n_points_;
    CubeIdx      offsets_[8];
};


//=============================================================================
} // namespace IsoEx
//=============================================================================

#if !defined(ISOEX_REGULAR_GRID_C)
#define ISOEX_REGULAR_GRID_TEMPLATES
#include "RegularGridT.cc"
#endif

#endif // ISOEX_REGULARGRID_HH defined
//=============================================================================

