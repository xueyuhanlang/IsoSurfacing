//=============================================================================
//
//  CLASS ImplicitGrid
//
//=============================================================================


#ifndef ISOEX_IMPLICITGRID_HH
#define ISOEX_IMPLICITGRID_HH


//== INCLUDES =================================================================

#include <IsoEx/Grids/RegularGridT.hh>
#include <IsoEx/Implicits/Implicit.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <vector>
#include "IsoSurfGen/TinyVector.h"

//== NAMESPACES ===============================================================

namespace IsoEx {

//== CLASS DEFINITION =========================================================

	      
/** \class ImplicitGrid ImplicitGrid.hh <IsoEx/Grids/ImplicitGrid.hh>
    
    This is the base class for all grids representing implicit
    objects, i.e. they store a reference to an implicit.  All
    inside/outside tests and distance queries will be passed on to the
    implicit.

    In addition, this grid also provides caching of inside/outside tests and
    scalar distance queries.

    \ingroup grids
*/	      
template <class Vec3 >
class ImplicitGrid : public RegularGrid< Vec3 >
{
public:
   
  typedef typename Vec3::value_type real;

  /// CubeIdx can be used to refer to cubes
  typedef typename Grid< Vec3 >::CubeIdx CubeIdx;

  /// The grid points can be referred to by PointIdx
  typedef typename Grid< Vec3 >::PointIdx PointIdx;

  //------------------------------------------------------------ public methods


  /// Default constructor
  ImplicitGrid(const Implicit<Vec3>& _implicit,
	       const Vec3&  _origin,
	       const Vec3&  _x_axis,
	       const Vec3&  _y_axis,
	       const Vec3&  _z_axis,
	       unsigned int            _x_res,
	       unsigned int            _y_res,
	       unsigned int            _z_res) 
    : RegularGrid< Vec3 >(_origin, _x_axis, _y_axis, _z_axis, _x_res, _y_res, _z_res),
      implicit_(_implicit) 
  {}

  /// Destructor
  virtual ~ImplicitGrid() {}



  /// \name This function calls will be passed to the implicit object.
  //@{

  /// See IsoEx::Implicit::is_inside()
  virtual bool is_inside(PointIdx _pidx, const real _isovalue = 0) const {
    return (is_inside_cache_.empty() ? implicit_.is_inside(this->point(_pidx), _isovalue) : is_inside_cache_[_pidx]);
  }

  /// See IsoEx::Implicit::scalar_distance()
  virtual real scalar_distance(PointIdx _pidx) const {
    return (scalar_distance_cache_.empty() ? implicit_.scalar_distance(this->point(_pidx)) : scalar_distance_cache_[_pidx]);
  }

  /// See IsoEx::Implicit::directed_distance()
  virtual bool directed_distance(const Vec3&  _p0,
                                 const Vec3&  _p1,
                                 Vec3&        _point,
                                 Vec3&        _normal,
                                 real&       _distance,
                                 const real _isovalue = 0
      ) const {
    return implicit_.directed_distance(_p0, _p1, _point, _normal, _distance, _isovalue);
  }

  virtual void directed_distance(
      const std::vector<TinyVector<real, 3>>& _p0vec,
      const std::vector<TinyVector<real, 3>>& _p1vec,
      const std::vector<real>& _p0vec_values,
      const std::vector<real>& _p1vec_values,
      std::vector<TinyVector<real, 3>>& _point,
      std::vector<TinyVector<real, 3>>& _normal,
      const real _isovalue = 0) const
  {
      return implicit_.directed_distance(_p0vec, _p1vec, _p0vec_values, _p1vec_values, _point, _normal, _isovalue);
  }

  
  //@}



  /// \name Enable caching of inside/outside/distance computations
  //@{

  /// Cache results of is_inside()
  void build_is_inside_cache(const real _isovalue = 0) const
  {
    int i, np( this->n_points() );
    
    is_inside_cache_.clear();
    is_inside_cache_.resize( this->n_points());

    for (i=0; i<np; ++i)
      is_inside_cache_[i] = implicit_.is_inside(this->point(i), _isovalue);
  }

  /// Cache results of scalar_distance()
  void build_scalar_distance_cache() const
  {
    int i, np(this->n_points());

    scalar_distance_cache_.clear();
    scalar_distance_cache_.resize(np);

    for (i=0; i<np; ++i)
      scalar_distance_cache_[i] = implicit_.scalar_distance(this->point(i));
  }


  std::vector<bool>& get_is_inside_cache() { return  is_inside_cache_; }
  const std::vector<bool>& get_is_inside_cache() const { return  is_inside_cache_; }
  std::vector<real>& get_scalar_distance_cache() { return  scalar_distance_cache_; }
  const std::vector<real>& get_scalar_distance_cache() const { return  scalar_distance_cache_; }

  //@}
  

 
protected:

  const Implicit<Vec3>&     implicit_;

  mutable std::vector<bool>   is_inside_cache_;
  mutable std::vector<real>  scalar_distance_cache_;
};


//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif // ISOEX_IMPLICITGRID_HH defined
//=============================================================================
