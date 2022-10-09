//=============================================================================
//
//  CLASS VectorFieldT
//
// Author:  Dominik Sibbing <sibbing@cs.rwth-aachen.de>
//
// Version: $Revision: 1$
// Date:    $Date: 2008$
//
//=============================================================================


#ifndef ACG_VECTORFIELDT_HH
#define ACG_VECTORFIELDT_HH


//== INCLUDES =================================================================

#include <OpenMesh/Core/Geometry/VectorT.hh>

#include "RegularGridT.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace IsoEx
{

//== CLASS DEFINITION =========================================================




/** \class VectorFieldT VectorFieldT.hh <ACG/.../VectorFieldT.hh>

    Brief Description.

    A more elaborate description follows.
*/

template <typename Scalar, int Dim, class Vec3>
class VectorFieldT : public RegularGrid<Vec3>
{
public:

  typedef typename Vec3::value_type real;

  /// CubeIdx can be used to refer to cubes
  typedef typename Grid< Vec3 >::CubeIdx CubeIdx;

  /// The grid points can be referred to by PointIdx
  typedef typename Grid< Vec3 >::PointIdx PointIdx;


    typedef OpenMesh::VectorT<Scalar, Dim> Vector;

    /// Default constructor
    VectorFieldT( const Vec3&  _origin = Vec3( 0, 0, 0 ),
                  const Vec3&  _x_axis = Vec3( 1, 0, 0 ),
                  const Vec3&  _y_axis = Vec3( 0, 1, 0 ),
                  const Vec3&  _z_axis = Vec3( 0, 0, 1 ),
                  unsigned int            _x_res  = 10,
                  unsigned int            _y_res  = 10,
                  unsigned int            _z_res  = 10 )
    {
        initialize( _origin, _x_axis, _y_axis, _z_axis, _x_res, _y_res, _z_res );
    }

    /// Destructor
    ~VectorFieldT() {}


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
        Vector vec;
        for ( int i = 0; i < vec.dim(); ++i ) vec[i] = 0.0;
        values_.resize( _x_res*_y_res*_z_res, vec );
    }


    virtual real scalar_distance( PointIdx /*_pidx*/ ) const
    {
        return 0.0;
    }

    virtual bool is_inside( PointIdx /*_pidx*/ ) const
    {
        return false;
    }

    virtual bool directed_distance( const Vec3&  /*_p0*/,
                                    const Vec3&  /*_p1*/,
                                    Vec3&        /*_point*/,
                                    Vec3&        /*_normal*/,
                                    real&                  /*_distance*/ ) const
    {
        return false;
    }


    /// data access
    Vector& operator()( unsigned int x, unsigned int y, unsigned int z )
    {
        return values_[x + y*this->x_resolution() + z*this->x_resolution()*this->y_resolution()];
    }

    Vector operator()( unsigned int x, unsigned int y, unsigned int z ) const
    {
        return values_[x + y*this->x_resolution() + z*this->x_resolution()*this->y_resolution()];
    }

    Vector& value( unsigned int x, unsigned int y, unsigned int z )
    {
        return ( *this )( x, y, z );
    }

    Vector value( unsigned int x, unsigned int y, unsigned int z ) const
    {
        return ( *this )( x, y, z );
    }


    /// get scalar value, returns zero vector if position is not in range
    Vector value_range( int x, int y, int z ) const;

    /// function to lineary interpolate a vector at a local point
    Vector lerp_local( Scalar _x, Scalar _y, Scalar _z );

    /// Function to lineary interploate a vector at a world point
    Vector lerp_world( Scalar _x, Scalar _y, Scalar _z );

private:

    /// vector of vectors
    std::vector< Vector > values_;
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_VECTORFIELDT_C)
#define ACG_VECTORFIELDT_TEMPLATES
#include "VectorFieldT.cc"
#endif
//=============================================================================
#endif // ACG_VECTORFIELDT_HH defined
//=============================================================================

