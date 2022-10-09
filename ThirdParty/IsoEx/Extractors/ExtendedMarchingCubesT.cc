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
//  CLASS ExtendedMarchingCubesT - IMPLEMENTATION
//
//=============================================================================

#define ISOEX_EXTMARCHINGCUBEST_C

//== INCLUDES =================================================================

#include <IsoEx/Extractors/ExtendedMarchingCubesT.hh>
#include <IsoEx/Extractors/MCTablesIsoEx.hh>
#include <IsoEx/Math/MatrixT.hh>
#include <IsoEx/Math/svd.hh>
#include <vector>
#include <float.h>
#include <iostream>

//== NAMESPACES ===============================================================

namespace IsoEx {

//== IMPLEMENTATION ==========================================================


template <class Mesh, class Grid>
ExtendedMarchingCubesT<Mesh,Grid>::
ExtendedMarchingCubesT(const Grid&  _grid,
		       Mesh&        _mesh,
               double       _iso_value,
		       double       _feature_angle)
  : grid_(_grid),
    mesh_(_mesh),
    feature_angle_(_feature_angle / 180.0 * M_PI),
    n_edges_(0),
    n_corners_(0)
{
  CubeIterator cube_it(grid_.begin()), cube_end(grid_.end());
  for (; cube_it!=cube_end; ++cube_it)
    process_cube(*cube_it, _iso_value);

  flip_edges();

  std::cerr << "Found "
	    << n_edges_ << " edge features, " 
	    << n_corners_ << " corner features\n";
}


//-----------------------------------------------------------------------------


template <class Mesh, class Grid>
void
ExtendedMarchingCubesT<Mesh,Grid>::
process_cube(CubeIdx _idx, double _iso_value)
{
  PointIdx           corner[8];
  VertexHandle       samples[12];
  unsigned char      cubetype(0);
  unsigned int       i, j;
  unsigned int       n_components, n_vertices;
  int                *indices;
  VertexHandle       vh;
  std::vector<VertexHandle> vhandles;



  // get corner vertices
  for (i=0; i<8; ++i)
    corner[i] = grid_.point_idx(_idx, i);


  // determine cube type
  for (i=0; i<8; ++i)
    if (!grid_.is_inside(corner[i], _iso_value))
      cubetype |= (1<<i);


  // trivial reject ?
  if (cubetype == 0 || cubetype == 255)
    return;


  // compute samples on cube's edges
  if (edgeTable[cubetype]&1)    samples[0]  = add_vertex(corner[0], corner[1], _iso_value);
  if (edgeTable[cubetype]&2)    samples[1]  = add_vertex(corner[1], corner[2], _iso_value);
  if (edgeTable[cubetype]&4)    samples[2]  = add_vertex(corner[3], corner[2], _iso_value);
  if (edgeTable[cubetype]&8)    samples[3]  = add_vertex(corner[0], corner[3], _iso_value);
  if (edgeTable[cubetype]&16)   samples[4]  = add_vertex(corner[4], corner[5], _iso_value);
  if (edgeTable[cubetype]&32)   samples[5]  = add_vertex(corner[5], corner[6], _iso_value);
  if (edgeTable[cubetype]&64)   samples[6]  = add_vertex(corner[7], corner[6], _iso_value);
  if (edgeTable[cubetype]&128)  samples[7]  = add_vertex(corner[4], corner[7], _iso_value);
  if (edgeTable[cubetype]&256)  samples[8]  = add_vertex(corner[0], corner[4], _iso_value);
  if (edgeTable[cubetype]&512)  samples[9]  = add_vertex(corner[1], corner[5], _iso_value);
  if (edgeTable[cubetype]&1024) samples[10] = add_vertex(corner[2], corner[6], _iso_value);
  if (edgeTable[cubetype]&2048) samples[11] = add_vertex(corner[3], corner[7], _iso_value);



  // connect samples by triangles
  n_components =  triTable[cubetype][1][0];
  indices      = &triTable[cubetype][1][n_components+1];

  for (i=1; i<=n_components; ++i)  // sheets in this voxel
  {
    // current sheet contains n_vertices vertices
    n_vertices = triTable[cubetype][1][i];


    // collect vertices of n-gon
    vhandles.clear();
    for (j=0; j<n_vertices; ++j)
      vhandles.push_back(samples[indices[j]]);

    
    // look for a feature
    vh = find_feature(vhandles);


    // feature -> create triangle fan around feature vertex
    if (vh.is_valid())
    {
      vhandles.push_back(vhandles[0]);
      for (j=0; j<n_vertices; ++j)
	mesh_.add_face(vhandles[j], vhandles[j+1], vh);
    }


    // no feature -> old marching cubes triangle table
    else
    {
      for (j=0; polyTable[n_vertices][j] != -1; j+=3)
	mesh_.add_face( samples[indices[polyTable[n_vertices][  j]]],
			samples[indices[polyTable[n_vertices][j+1]]],
			samples[indices[polyTable[n_vertices][j+2]]] );
    }
    
    indices += n_vertices;
  }
}


//-----------------------------------------------------------------------------


template <class Mesh, class Grid>
typename ExtendedMarchingCubesT<Mesh,Grid>::VertexHandle
ExtendedMarchingCubesT<Mesh,Grid>::
add_vertex(PointIdx _p0, PointIdx _p1, double _iso_value)
{
  // find vertex if it has been computed already
  VertexHandle   vh = edge2vertex_.find(_p0, _p1);
  if (vh.is_valid())  return vh;



  // generate new vertex
  const typename Mesh::Point  p0(grid_.point(_p0));
  const typename Mesh::Point  p1(grid_.point(_p1));
  typename Mesh::Point        point, normal(0,0,0);
  typename Mesh::Point::value_type   distance;

  bool ok = grid_.directed_distance(p0, p1, point, normal, distance, _iso_value);
  if (!ok)
  {
    // should not happen, just in case of precision errors...
    float s0 = fabs(grid_.scalar_distance(_p0));
    float s1 = fabs(grid_.scalar_distance(_p1));
    float t  = s0 / (s0+s1);
    point = (1.0f-t)*p0 + t*p1;
  }

  
  // add vertex
  vh = mesh_.add_vertex(point);
  mesh_.set_normal(vh, normal);
  edge2vertex_.insert(_p0, _p1, vh);


  return vh;
}


//-----------------------------------------------------------------------------


template <class Mesh, class Grid>
typename ExtendedMarchingCubesT<Mesh,Grid>::VertexHandle
ExtendedMarchingCubesT<Mesh,Grid>::
find_feature(const VertexHandleVector& _vhandles)
{
  unsigned int i, j, nV = _vhandles.size(), rank;



  // collect point & normals;
  std::vector<typename Mesh::Point>  p, n;
  p.resize(nV); 
  n.resize(nV);
  for (i=0; i<nV; ++i)
  {
    p[i] = mesh_.point(_vhandles[i]);
    n[i] = mesh_.normal(_vhandles[i]);
  }
    


  // move barycenter of points into
  typename Mesh::Point cog(0,0,0);
  for (i=0; i<nV; ++i)  cog += p[i];
  cog /= (float)nV;
  for (i=0; i<nV; ++i)  p[i] -= cog;



  // normal angle criterion
  double  c, min_c, max_c;
  typename Mesh::Point  axis;
  for (min_c=1.0, i=0; i<nV; ++i)
    for (j=0; j<nV; ++j)
      if ((c = (n[i] | n[j])) < min_c)
      {
	min_c = c;
	axis  = n[i] % n[j];
      }



  // angle to small, no feature -> return invalid vertex handle
  if (min_c > cos(feature_angle_)) 
    return Mesh::InvalidVertexHandle; 



  // ok, we have a feature
  // is it edge or corner, i.e. rank 2 or 3 ?
  axis.normalize();
  for (min_c=1.0, max_c=-1.0, i=0; i<nV; ++i)
  {
    c = (axis | n[i]);
    if (c < min_c)  min_c = c;
    if (c > max_c)  max_c = c;
  }
  c = std::max(fabs(min_c),fabs(max_c));
  c = sqrt(1.0-c*c);
  rank = (c > cos(feature_angle_) ? 2 : 3);

  if (rank == 2)  ++n_edges_;
  else            ++n_corners_;




  // setup linear system (find intersection of tangent planes)
  Math::MatrixT<double>  A(nV, 3);
  Math::VectorT<double>  b(nV);

  for (i = 0; i < nV; ++i) {
    A(i, 0) = n[i][0];
    A(i, 1) = n[i][1];
    A(i, 2) = n[i][2];
    b(i)    = (p[i] | n[i]);
  }



  // SVD of matrix A
  Math::MatrixT<double>  V(3, 3);
  Math::VectorT<double>  S(nV);
  Math::svd_decomp(A, S, V);


  // rank == 2 -> suppress smallest singular value
  if (rank == 2) {
    double smin = FLT_MAX;
    unsigned int sminid = 0;
    unsigned int srank = std::min(nV, 3u);

    for (i = 0; i < srank; ++i)
      if (S(i) < smin) {
        smin = S(i);
        sminid = i;
      }

    S(sminid) = 0.0;
  }



  // SVD backsubstitution -> least squares, least norm solution x
  Math::VectorT<double>  x(3);
  Math::svd_backsub(A, S, V, b, x);



  // transform x to world coords
  typename Mesh::Point point(x(0), x(1), x(2));
  point += cog;



  // insert the feature-point 
  VertexHandle vh = mesh_.add_vertex(point);
  mesh_.status(vh).set_feature(true);


  return vh;
}


//-----------------------------------------------------------------------------


template<class Mesh, class Grid>
void ExtendedMarchingCubesT<Mesh,Grid>::flip_edges()
{
  VertexHandle v0, v1, v2, v3;
  typename Mesh::HalfedgeHandle he;

  typename Mesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end());

  for (; e_it != e_end; ++e_it) {
    if (mesh_.is_flip_ok(*e_it)) {
      he = mesh_.halfedge_handle(*e_it, 0);
      v0 = mesh_.to_vertex_handle(he);
      he = mesh_.next_halfedge_handle(he);
      v1 = mesh_.to_vertex_handle(he);
      he = mesh_.halfedge_handle(*e_it, 1);
      v2 = mesh_.to_vertex_handle(he);
      he = mesh_.next_halfedge_handle(he);
      v3 = mesh_.to_vertex_handle(he);

      // flip edge if it would connect two features (v1, v3)
      // and not disconnect two others (v0, v2) afterwards
      // (maybe we should check for flipping triangle normals)
      if (mesh_.status(v1).feature() && mesh_.status(v3).feature() && !mesh_.status(v0).feature() && !mesh_.status(v2).feature())
        mesh_.flip(*e_it);
    }
  }
}


//=============================================================================
} // namespace IsoEx
//=============================================================================
