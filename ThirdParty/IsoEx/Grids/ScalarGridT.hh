
//=============================================================================
//
//  CLASS ScalarGridT
//
//=============================================================================


#ifndef ISOEX_SCALARGRIDT_HH
#define ISOEX_SCALARGRIDT_HH


//== INCLUDES =================================================================

#include <OpenMesh/Core/IO/BinaryHelper.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>

#include "RegularGridT.hh"
#include "../Implicits/Implicit.hh"

#include <vector>
#include <iostream>

//== NAMESPACES ===============================================================

namespace IsoEx
{

//== CLASS DEFINITION =========================================================


/** \class ScalarGridT ScalarGridT.hh <IsoEx/Grids/ScalarGridT.hh>

    This grid just holds a vector of scalar values.

    \ingroup grids
*/
template <typename Scalar,class Vec3 >
class ScalarGridT : public RegularGrid<Vec3>
{
public:

  typedef typename Vec3::value_type real;

  /// CubeIdx can be used to refer to cubes
  typedef typename Grid< Vec3 >::CubeIdx CubeIdx;

  /// The grid points can be referred to by PointIdx
  typedef typename Grid< Vec3 >::PointIdx PointIdx;

    //typedef _Scalar Scalar;
    typedef std::vector<Scalar>  		Values;

    /// Default constructor
    ScalarGridT( const Vec3&  _origin = Vec3( 0, 0, 0 ),
                 const Vec3&  _x_axis = Vec3( 1, 0, 0 ),
                 const Vec3&  _y_axis = Vec3( 0, 1, 0 ),
                 const Vec3&  _z_axis = Vec3( 0, 0, 1 ),
                 unsigned int            _x_res  = 1,
                 unsigned int            _y_res  = 1,
                 unsigned int            _z_res  = 1 )
    {
      initialize( _origin, _x_axis, _y_axis, _z_axis, _x_res, _y_res, _z_res );
    }

    /// Destructor
    virtual ~ScalarGridT() {}

    /// function to initialize the grid
    void initialize( const Vec3&  _origin,
                     const Vec3&  _x_axis,
                     const Vec3&  _y_axis,
                     const Vec3&  _z_axis,
                     unsigned int            _x_res,
                     unsigned int            _y_res,
                     unsigned int            _z_res )
    {
        RegularGrid<Vec3>::initialize( _origin, _x_axis, _y_axis, _z_axis,
                                       _x_res, _y_res, _z_res );
        values_.resize( this->n_points_, 0.0 );
    }



    virtual real scalar_distance( PointIdx _pidx ) const
    {
        return values_[_pidx];
    }

    virtual bool is_inside( PointIdx _pidx, const real _isovalue = 0 ) const
    {
        return values_[_pidx] < _isovalue;
    }

    virtual bool directed_distance( const Vec3&  /*_p0*/,
                                    const Vec3&  /*_p1*/,
                                    Vec3&        /*_point*/,
                                    Vec3&        /*_normal*/,
                                    real&                  /*_distance*/,
                                    const real _isovalue = 0
                                   ) const
    {
        return false;
    }

    void sample( const Implicit<Vec3>& _implicit );


    virtual bool read( const char* _filename );
    virtual bool write( const char* _filename );
    virtual bool read( FILE* _in );
    virtual bool write( FILE* _out );

    /// data access
    Scalar& operator()( unsigned int x, unsigned int y, unsigned int z )
    {
        return values_[x + y*this->x_resolution() + z*this->x_resolution()*this->y_resolution()];
    }

    Scalar operator()( unsigned int x, unsigned int y, unsigned int z ) const
    {
        return values_[x + y*this->x_resolution() + z*this->x_resolution()*this->y_resolution()];
    }

    Scalar& value( unsigned int x, unsigned int y, unsigned int z )
    {
        return ( *this )( x, y, z );
    }

    Scalar value( unsigned int x, unsigned int y, unsigned int z ) const
    {
        return ( *this )( x, y, z );
    }

    /// get scalar value, returns 0.0 if position is not in range
    Scalar value_range( int x, int y, int z ) const;

//     void resize() {values_ = Values( n_points(), 0 );}  // changed

    /// function to lineary interpolate a scalar value at a local point
    Scalar lerp_local( Scalar _x, Scalar _y, Scalar _z );

    /// funciton to lineary interploate a scalar value at a world point
    Scalar lerp_world( Scalar _x, Scalar _y, Scalar _z );

    /// get intersections with isosurface in local coordinates
    void get_isosurface_intersections_local( const Vec3 &_o,
            const Vec3 &_d,
            Scalar _iso,
            std::vector< Vec3 > &_intersections );

    /// get intersections with isosurface in world coordinates
    void get_isosurface_intersections_world( const Vec3 &_o,
            const Vec3 &_d,
            Scalar _iso,
            std::vector< Vec3 > &_intersections );

    Values& get_values() { return values_; }
    const Values& get_values() const { return values_; }
protected:

    Values  values_;
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ISOEX_SCALARGRIDT_C)
#define ISOEX_SCALARGRIDT_TEMPLATES
#include "ScalarGridT.cc"
#endif
//=============================================================================
#endif // ISOEX_SCALARGRIDT_HH defined
//=============================================================================

