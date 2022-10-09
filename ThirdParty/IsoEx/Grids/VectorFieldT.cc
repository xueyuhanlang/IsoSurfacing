//=============================================================================
//
//  CLASS VectorFieldT - IMPLEMENTATION
//
//=============================================================================

#define ACG_VECTORFIELDT_C

//== INCLUDES =================================================================

#include "VectorFieldT.hh"

//== NAMESPACES ===============================================================

namespace IsoEx
{

//== IMPLEMENTATION ==========================================================



template <typename Scalar, int Dim, class Vec3>
typename VectorFieldT<Scalar, Dim, Vec3>::Vector
VectorFieldT< Scalar, Dim, Vec3 > ::
value_range( int x, int y, int z ) const
{
    Vector vec;
    if ( x < 0 || x >= ( int )this->x_resolution() ||
            y < 0 || y >= ( int )this->y_resolution() ||
            z < 0 || z >= ( int )this->z_resolution() )
    {
        for ( int i = 0; i < vec.dim(); ++i )
            vec[i] = 0;
    }
    else
        vec = value(( unsigned int )x, ( unsigned int )y, ( unsigned int )z );

    return vec;
}

//-----------------------------------------------------------------------------

template< typename Scalar, int Dim , class Vec3>
typename VectorFieldT<Scalar, Dim, Vec3>::Vector
VectorFieldT<Scalar, Dim, Vec3>::
lerp_local( Scalar _x, Scalar _y, Scalar _z )
{
    // grid spacing
    Scalar dxf = this->dx().norm();
    Scalar dyf = this->dy().norm();
    Scalar dzf = this->dz().norm();

    // get gird cell
    int x = ( int )( _x / dxf );
    int y = ( int )( _y / dyf );
    int z = ( int )( _z / dzf );

    // calculate interpalation parameters
    Scalar u = std::max(( _x / dxf - Scalar( x ) ), Scalar(0) );
    Scalar v = std::max(( _y / dyf - Scalar( y ) ), Scalar(0) );
    Scalar w = std::max(( _z / dzf - Scalar( z ) ), Scalar(0) );

    // get values
    Vector c0 = value_range( x  , y  , z );
    Vector c1 = value_range( x + 1, y  , z );
    Vector c2 = value_range( x  , y + 1, z );
    Vector c3 = value_range( x + 1, y + 1, z );

    Vector c4 = value_range( x  , y  , z + 1 );
    Vector c5 = value_range( x + 1, y  , z + 1 );
    Vector c6 = value_range( x  , y + 1, z + 1 );
    Vector c7 = value_range( x + 1, y + 1, z + 1 );

    // interpolate
    return   c0 * ( 1.0 - u ) * ( 1.0 - v ) * ( 1.0 - w )
             + c1 *     u * ( 1.0 - v ) * ( 1.0 - w )
             + c2 * ( 1.0 - u ) * v * ( 1.0 - w )
             + c3 * u * v * ( 1.0 - w ) +
             c4 * ( 1.0 - u ) * ( 1.0 - v ) * ( w )
             + c5 * u * ( 1.0 - v ) * ( w )
             + c6 * ( 1.0 - u ) * v * ( w )
             + c7 * u * v * ( w );
}

//-----------------------------------------------------------------------------


template< typename Scalar, int Dim , class Vec3>
typename VectorFieldT<Scalar, Dim, Vec3>::Vector
VectorFieldT<Scalar, Dim, Vec3>::
lerp_world( Scalar _x, Scalar _y, Scalar _z )
{
    Vec3 pl = to_local( Vec3( _x, _y, _z ) );
    return lerp_local( pl[0], pl[1], pl[2] );
}


//=============================================================================
} // namespace ACG
//=============================================================================
