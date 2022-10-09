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
//  CLASS ExtendedMarchingCubesT
//
//=============================================================================

#ifndef ISOEX_AEXTMARCHINGCUBEST_HH
#define ISOEX_AEXTMARCHINGCUBEST_HH

//== INCLUDES =================================================================

#include <IsoEx/Extractors/Edge2VertexMapT.hh>
#include <IsoEx/Grids/Grid.hh>
#include <vector>
#include "IsoSurfGen/ImplicitFunc.h"

//== NAMESPACES ===============================================================

namespace IsoEx
{

  //== CLASS DEFINITION =========================================================

  /** \class AExtendedMarchingCubesT ExtendedMarchingCubesT.hh <IsoEx/Extractors/ExtendedMarchingCubesT.hh>

      This class implements the Extended Marching Cubes of Kobbelt et al,
      Siggraph 2001.

      modified by Yang Liu for parallelizing the algorithm

      The 0-level iso-surface is extracted in the constructor. Use it through
      the convenience function
      <b>IsoEx::advanced_extended_marching_cubes()</b>.

      \ingroup extractors
  */
  template <class Mesh, class Grid, class Real>
  class AExtendedMarchingCubesT
  {
  public:
    AExtendedMarchingCubesT(
        const Grid &_grid,
        Mesh &_mesh,
        double _iso_value,
        double _feature_angle,
        const int buffer,
        bool verbose);

  private:
    typedef typename Grid::PointIdx PointIdx;
    typedef typename Grid::CubeIdx CubeIdx;
    typedef typename Grid::CubeIterator CubeIterator;
    typedef typename Mesh::VertexHandle VertexHandle;
    typedef std::vector<VertexHandle> VertexHandleVector;

    VertexHandle add_vertex(PointIdx _p0, PointIdx _p1, const TinyVector<Real, 3> &_point, const TinyVector<Real, 3> &_normal);
    VertexHandle find_feature(const VertexHandleVector &_vhandles);

    bool is_unprocessed_edge(PointIdx _p0, PointIdx _p1);
    int query_store_id(PointIdx _p0, PointIdx _p1);

    void flip_edges();

    const Grid &grid_;
    Mesh &mesh_;

    float feature_angle_;
    unsigned int n_edges_, n_corners_;

    // maps an edge to the sample vertex generated on it
    Edge2VertexMapT<PointIdx, VertexHandle> edge2vertex_;
    std::map<std::pair<PointIdx, PointIdx>, unsigned int> edgemap;
  };

  //-----------------------------------------------------------------------------

  /** Convenience wrapper for the Extended Marching Cubes algorithm.
      \see IsoEx::ExtendedMarchingCubesT
      \ingroup extractors
  */
  template <class Mesh, class Grid, class Real>
  void advanced_extended_marching_cubes(const Grid &_grid,
                                        Mesh &_mesh,
                                        double _iso_value,
                                        double _feature_angle,
                                        const int buffer,
                                        bool verbose)
  {
    AExtendedMarchingCubesT<Mesh, Grid, Real> emc(_grid, _mesh, _iso_value, _feature_angle, buffer, verbose);
  }

  //=============================================================================
} // namespace IsoEx
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ISOEX_AEXTMARCHINGCUBEST_C)
#define ISOEX_AEXTMARCHINGCUBEST_TEMPLATES
#include "AExtendedMarchingCubesT.cc"
#endif
//=============================================================================
#endif // ISOEX_AEXTMARCHINGCUBEST_HH defined
//=============================================================================
