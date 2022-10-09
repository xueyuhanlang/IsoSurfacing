//=============================================================================
//
//  CLASS ScalarGridT - IMPLEMENTATION
//
//=============================================================================

#define ISOEX_SCALARGRIDT_C

//== INCLUDES =================================================================

#include "ScalarGridT.hh"


//== NAMESPACES ===============================================================

namespace IsoEx
{

//== IMPLEMENTATION ==========================================================


template <class Scalar,class Vec3>
void
ScalarGridT<Scalar,Vec3>::
sample( const Implicit<Vec3>& _implicit )
{
    for ( unsigned int i = 0; i < this->n_points(); ++i )
        values_[i] = _implicit.scalar_distance( this->point( i ) );
}


//-----------------------------------------------------------------------------


template <class Scalar,class Vec3>
bool
ScalarGridT<Scalar,Vec3>::
read( const char* _filename )
{
    bool ok = false;
    FILE* in = fopen( _filename, "rb" );
    if ( in )
    {
        ok = read( in );
        fclose( in );
    }
    return ok;
}


template <class Scalar,class Vec3>
bool
ScalarGridT<Scalar,Vec3>::
read( FILE* _in )
{
    // header
    unsigned int     x_res, y_res, z_res;

    float a,b,c,d,e,f,g,h,i,j,k,l;
    fscanf( _in,
            "%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%d %d %d\n",
            &a, &b, &c,
            &d, &e, &f,
            &g, &h, &i,
            &j, &k, &l,
            &x_res, &y_res, &z_res );

    Vec3  origin( a,b,c );
    Vec3  x_axis( d,e,f );
    Vec3  y_axis( g,h,i );
    Vec3  z_axis( j,k,l );
    initialize( origin, x_axis, y_axis, z_axis, x_res, y_res, z_res );


    // values
    values_ = Values( this->n_points(), 0 );
    for ( unsigned int i = 0; i < this->n_points(); ++i )
        values_[i] = OpenMesh::IO::read_float( _in );


    return true;
}


//-----------------------------------------------------------------------------


template <class Scalar,class Vec3>
bool
ScalarGridT<Scalar,Vec3>::
write( const char* _filename )
{
    bool ok = false;
    FILE* out = fopen( _filename, "wb" );
    if ( out )
    {
        ok = write( out );
        fclose( out );
    }
    return ok;
}


template <class Scalar,class Vec3>
bool
ScalarGridT<Scalar,Vec3>::
write( FILE* _out )
{
    // header
    fprintf( _out,
             "%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%d %d %d\n",
             this->origin()[0], this->origin()[1], this->origin()[2],
             this->x_axis()[0], this->x_axis()[1], this->x_axis()[2],
             this->y_axis()[0], this->y_axis()[1], this->y_axis()[2],
             this->z_axis()[0], this->z_axis()[1], this->z_axis()[2],
             this->x_resolution(), this->y_resolution(), this->z_resolution() );


    // values
    for ( unsigned int i = 0; i < this->n_points(); ++i )
        OpenMesh::IO::write_float( values_[i], _out );


    return true;
}


//-----------------------------------------------------------------------------


template <class Scalar,class Vec3>
Scalar
ScalarGridT<Scalar,Vec3>::
value_range( int x, int y, int z ) const
{
    if ( x < 0 || x >= ( int )this->x_resolution() ||
            y < 0 || y >= ( int )this->y_resolution() ||
            z < 0 || z >= ( int )this->z_resolution() )
        return 0.0;
    else
        return value(( unsigned int )x, ( unsigned int )y, ( unsigned int )z );
}


//-----------------------------------------------------------------------------


template <class Scalar,class Vec3>
Scalar
ScalarGridT<Scalar,Vec3>::
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
    Scalar c0 = value_range( x  , y  , z );
    Scalar c1 = value_range( x + 1, y  , z );
    Scalar c2 = value_range( x  , y + 1, z );
    Scalar c3 = value_range( x + 1, y + 1, z );

    Scalar c4 = value_range( x  , y  , z + 1 );
    Scalar c5 = value_range( x + 1, y  , z + 1 );
    Scalar c6 = value_range( x  , y + 1, z + 1 );
    Scalar c7 = value_range( x + 1, y + 1, z + 1 );

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


template <class Scalar,class Vec3>
Scalar
ScalarGridT<Scalar,Vec3>::
lerp_world( Scalar _x, Scalar _y, Scalar _z )
{
    Vec3 pl = to_local( Vec3( _x, _y, _z ) );
    return lerp_local( pl[0], pl[1], pl[2] );
}


//-----------------------------------------------------------------------------

template <class Scalar,class Vec3>
void
ScalarGridT<Scalar,Vec3>::
get_isosurface_intersections_world( const Vec3 &_o, const Vec3 &_d, Scalar _iso,
                                    std::vector< Vec3 > &_intersections )
{
    Vec3 o = this->to_local_.transform_point( _o );
    Vec3 d = this->to_local_.transform_vector( _d );

    get_isosurface_intersections_local( o, d, _iso, _intersections );

    // transform to world coordinates
    for ( unsigned int i = 0; i < _intersections.size(); ++i )
        _intersections[i] = this->to_world_.transform_point( _intersections[i] );
}


//-----------------------------------------------------------------------------


template <class Scalar,class Vec3>
void
ScalarGridT<Scalar,Vec3>::
get_isosurface_intersections_local( const Vec3 &_o, const Vec3 &_d, Scalar _iso,
                                    std::vector< Vec3 > &_intersections )
{
    _intersections.clear();

    // find cube entry and exit points
    Vec3 entry, exit;
    if ( !ray_intersect_local( _o, _d, entry, exit ) )
    {
        std::cout << "Cube not hit. " << std::endl;
        return;
    }

    // grid spacing
    Scalar dxf = this->dx().norm();
    Scalar dyf = this->dy().norm();
    Scalar dzf = this->dz().norm();

    // get direction offset
    Vec3 ddir = _d;
    ddir.normalize();
    Scalar norm_ddir = std::min( dxf, std::min( dyf, dzf ) );
    ddir *= norm_ddir;

    // calculate number of steps needed to traverse the volume
    unsigned int n_steps(( exit-entry ).norm() / norm_ddir );

    Scalar prev_iso(0);
    Scalar cur_iso;

    Vec3 cur_pos = entry;

    // traverse volume
    for ( unsigned int i = 0; i < n_steps; ++i )
    {
        // get current scalar value
        cur_iso = lerp_local( Scalar( cur_pos[0] ), Scalar( cur_pos[1] ), Scalar( cur_pos[2] ) );
        if ( i == 0 ) prev_iso = cur_iso;

        // intersect isosurface?
        if ((( prev_iso < _iso ) && ( _iso < cur_iso ) ) ||
                (( prev_iso > _iso ) && ( _iso > cur_iso ) ) )
        {
            // refinement of position
            Scalar t = ( _iso - prev_iso )/( cur_iso - prev_iso );
            _intersections.push_back( cur_pos + ddir*(t-1.0) );
        }

        prev_iso = cur_iso;
        cur_pos += ddir;
    }
}


//=============================================================================
} // namespace IsoEx
//=============================================================================
