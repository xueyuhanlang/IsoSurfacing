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

#define ISOEX_AOEXTMARCHINGCUBEST_C

//== INCLUDES =================================================================

#include <IsoEx/Extractors/AOExtendedMarchingCubesT.hh>
#include <IsoEx/Extractors/MCTablesIsoEx.hh>
#include <IsoEx/Math/MatrixT.hh>
#include <IsoEx/Math/svd.hh>
#include <vector>
#include <float.h>
#include <iostream>
#include <omp.h>
#include "IsoSurfGen/ProgressBar.h"

//== NAMESPACES ===============================================================

namespace IsoEx
{
	//== IMPLEMENTATION ==========================================================

	template <class Mesh, class Grid, class Real>
	AOExtendedMarchingCubesT<Mesh, Grid, Real>::
		AOExtendedMarchingCubesT(const Grid &_grid,
								 Mesh &_mesh,
								 double iso_value,
								 double _feature_angle,
								 const std::vector<int> &cache_cube,
								 const std::unordered_map<int, Real> &point_value_map,
								 const int buffer, bool verbose)
		: grid_(_grid),
		  mesh_(_mesh),
		  feature_angle_(_feature_angle / 180.0 * M_PI),
		  n_edges_(0),
		  n_corners_(0)
	{
		int xsize = grid_.x_resolution() - 1;
		int ysize = grid_.y_resolution() - 1;
		int zsize = grid_.z_resolution() - 1;

		/////////////////////////////////////////////////////////////////////
		std::vector<TinyVector<Real, 3>> edge_start, edge_end, intersection, normals;
		std::vector<Real> start_fvalue, end_fvalue;
		edge_start.reserve(buffer + 12);
		edge_end.reserve(buffer + 12);
		intersection.reserve(buffer + 12);
		normals.reserve(buffer + 12);
		start_fvalue.reserve(buffer + 12);
		end_fvalue.reserve(buffer + 12);

		print_progress_bar(0.0f, verbose);

		int counter = 0;
		int prev_start = 0;

		TinyVector<Real, 3> ipoint[8];
		int label[12] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
		int edgepair[12][2] = {{0, 1}, {1, 2}, {3, 2}, {0, 3}, {4, 5}, {5, 6}, {7, 6}, {4, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
		PointIdx corner[8];
		Real values[8];

		for (int s = 0; s < cache_cube.size(); s += 2)
		{
			int c = cache_cube[s];
			for (int i = 0; i < 8; ++i)
			{
				corner[i] = grid_.point_idx(c, i);
				values[i] = point_value_map.find(corner[i])->second;
			}

			unsigned char cubetype = cache_cube[s + 1];

			for (int i = 0; i < 8; i++)
			{
				const typename Mesh::Point &p = grid_.point(corner[i]);
				ipoint[i][0] = p[0];
				ipoint[i][1] = p[1];
				ipoint[i][2] = p[2];
			}

			for (int i = 0; i < 12; i++)
			{
				if ((edgeTable[cubetype] & label[i]) && is_unprocessed_edge(corner[edgepair[i][0]], corner[edgepair[i][1]]))
				{
					edge_start.push_back(ipoint[edgepair[i][0]]);
					edge_end.push_back(ipoint[edgepair[i][1]]);
					start_fvalue.push_back(values[edgepair[i][0]]);
					end_fvalue.push_back(values[edgepair[i][1]]);
					counter++;
				}
			}

			if (counter != 0 && (counter >= buffer || s + 2 == cache_cube.size()))
			{
				edge_start.resize(counter);
				edge_end.resize(counter);
				start_fvalue.resize(counter);
				end_fvalue.resize(counter);
				intersection.resize(counter);
				normals.resize(counter);

				grid_.directed_distance(edge_start, edge_end, start_fvalue, end_fvalue, intersection, normals, iso_value);

				int curcounter = 0;
				for (int n = prev_start; n <= s; n += 2)
				{
					CubeIdx _idx = cache_cube[n];
					PointIdx corner[8];
					VertexHandle samples[12];
					unsigned char cubetype(0);
					unsigned int i, j;
					unsigned int n_components, n_vertices;
					int *indices;
					VertexHandle vh;
					std::vector<VertexHandle> vhandles;

					// get corner vertices
					for (i = 0; i < 8; ++i)
						corner[i] = grid_.point_idx(_idx, i);

					cubetype = cache_cube[n + 1];

					for (int i = 0; i < 12; i++)
					{
						if ((edgeTable[cubetype] & label[i]))
						{
							int mid = query_store_id(corner[edgepair[i][0]], corner[edgepair[i][1]]);
							samples[i] = add_vertex(corner[edgepair[i][0]], corner[edgepair[i][1]], intersection[mid], normals[mid]);
						}
					}

					// connect samples by triangles
					n_components = triTable[cubetype][1][0];
					indices = &triTable[cubetype][1][n_components + 1];

					for (i = 1; i <= n_components; ++i) // sheets in this voxel
					{
						// current sheet contains n_vertices vertices
						n_vertices = triTable[cubetype][1][i];

						// collect vertices of n-gon
						vhandles.clear();
						for (j = 0; j < n_vertices; ++j)
							vhandles.push_back(samples[indices[j]]);

						// look for a feature
						vh = find_feature(vhandles);

						// feature -> create triangle fan around feature vertex
						if (vh.is_valid())
						{
							vhandles.push_back(vhandles[0]);
							for (j = 0; j < n_vertices; ++j)
								mesh_.add_face(vhandles[j], vhandles[j + 1], vh);
						}
						// no feature -> old marching cubes triangle table
						else
						{
							for (j = 0; polyTable[n_vertices][j] != -1; j += 3)
								mesh_.add_face(samples[indices[polyTable[n_vertices][j]]],
											   samples[indices[polyTable[n_vertices][j + 1]]],
											   samples[indices[polyTable[n_vertices][j + 2]]]);
						}

						indices += n_vertices;
					}
				}
				prev_start = s + 2;
				counter = 0;

				edgemap.clear();
				edge_start.resize(0);
				edge_end.resize(0);
				start_fvalue.resize(0);
				end_fvalue.resize(0);
				intersection.resize(0);
				normals.resize(0);

				print_progress_bar((s + 2.0f) / cache_cube.size(), verbose);
				if (s + 2 == cache_cube.size())
					printf("\n");
			}
		}
		/////////////////////////////////////////////////////////////////////

		flip_edges();

		std::cerr << "\nFound "
				  << n_edges_ << " edge features, "
				  << n_corners_ << " corner features\n";
	}
	//-----------------------------------------------------------------------------

	template <class Mesh, class Grid, class Real>
	typename AOExtendedMarchingCubesT<Mesh, Grid, Real>::VertexHandle
	AOExtendedMarchingCubesT<Mesh, Grid, Real>::
		add_vertex(PointIdx _p0, PointIdx _p1, const TinyVector<Real, 3> &_point, const TinyVector<Real, 3> &_normal)
	{
		// find vertex if it has been computed already
		VertexHandle vh = edge2vertex_.find(_p0, _p1);
		if (vh.is_valid())
			return vh;

		typename Mesh::Point point(_point[0], _point[1], _point[2]), normal(_normal[0], _normal[1], _normal[2]);

		// add vertex
		vh = mesh_.add_vertex(point);
		mesh_.set_normal(vh, normal);
		edge2vertex_.insert(_p0, _p1, vh);

		return vh;
	}

	//-----------------------------------------------------------------------------
	template <class Mesh, class Grid, class Real>
	bool
	AOExtendedMarchingCubesT<Mesh, Grid, Real>::
		is_unprocessed_edge(PointIdx _p0, PointIdx _p1)
	{
		std::pair<PointIdx, PointIdx> edge(_p0, _p1);
		auto iter = edgemap.find(edge);
		if (iter == edgemap.end())
		{
			unsigned int size = (unsigned int)edgemap.size();
			edgemap[edge] = size;
			return true;
		}

		return false;
	}

	//-----------------------------------------------------------------------------
	template <class Mesh, class Grid, class Real>
	int
	AOExtendedMarchingCubesT<Mesh, Grid, Real>::
		query_store_id(PointIdx _p0, PointIdx _p1)
	{
		std::pair<PointIdx, PointIdx> edge(_p0, _p1);
		auto iter = edgemap.find(edge);
		if (iter == edgemap.end())
		{
			std::cout << "This query is illegal!" << std::endl;
			return -1;
		}

		return iter->second;
	}

	//-----------------------------------------------------------------------------

	template <class Mesh, class Grid, class Real>
	typename AOExtendedMarchingCubesT<Mesh, Grid, Real>::VertexHandle
	AOExtendedMarchingCubesT<Mesh, Grid, Real>::
		find_feature(const VertexHandleVector &_vhandles)
	{
		unsigned int i, j, nV = _vhandles.size(), rank;

		// collect point & normals;
		std::vector<typename Mesh::Point> p, n;
		p.resize(nV);
		n.resize(nV);
		for (i = 0; i < nV; ++i)
		{
			p[i] = mesh_.point(_vhandles[i]);
			n[i] = mesh_.normal(_vhandles[i]);
		}

		// move barycenter of points into
		typename Mesh::Point cog(0, 0, 0);
		for (i = 0; i < nV; ++i)
			cog += p[i];
		cog /= (float)nV;
		for (i = 0; i < nV; ++i)
			p[i] -= cog;

		// normal angle criterion
		double c, min_c, max_c;
		typename Mesh::Point axis;
		for (min_c = 1.0, i = 0; i < nV; ++i)
			for (j = 0; j < nV; ++j)
				if ((c = (n[i] | n[j])) < min_c)
				{
					min_c = c;
					axis = n[i] % n[j];
				}

		// angle to small, no feature -> return invalid vertex handle
		if (min_c > cos(feature_angle_))
			return Mesh::InvalidVertexHandle;

		// ok, we have a feature
		// is it edge or corner, i.e. rank 2 or 3 ?
		axis.normalize();
		for (min_c = 1.0, max_c = -1.0, i = 0; i < nV; ++i)
		{
			c = (axis | n[i]);
			if (c < min_c)
				min_c = c;
			if (c > max_c)
				max_c = c;
		}
		c = std::max(fabs(min_c), fabs(max_c));
		c = sqrt(1.0 - c * c);
		rank = (c > cos(feature_angle_) ? 2 : 3);

		if (rank == 2)
			++n_edges_;
		else
			++n_corners_;

		// setup linear system (find intersection of tangent planes)
		Math::MatrixT<double> A(nV, 3);
		Math::VectorT<double> b(nV);

		for (i = 0; i < nV; ++i)
		{
			A(i, 0) = n[i][0];
			A(i, 1) = n[i][1];
			A(i, 2) = n[i][2];
			b(i) = (p[i] | n[i]);
		}

		// SVD of matrix A
		Math::MatrixT<double> V(3, 3);
		Math::VectorT<double> S(nV);
		Math::svd_decomp(A, S, V);

		// rank == 2 -> suppress smallest singular value
		if (rank == 2)
		{
			double smin = FLT_MAX;
			unsigned int sminid = 0;
			unsigned int srank = std::min(nV, 3u);

			for (i = 0; i < srank; ++i)
				if (S(i) < smin)
				{
					smin = S(i);
					sminid = i;
				}

			S(sminid) = 0.0;
		}

		// SVD backsubstitution -> least squares, least norm solution x
		Math::VectorT<double> x(3);
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

	template <class Mesh, class Grid, class Real>
	void AOExtendedMarchingCubesT<Mesh, Grid, Real>::flip_edges()
	{
		VertexHandle v0, v1, v2, v3;
		typename Mesh::HalfedgeHandle he;

		typename Mesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end());

		for (; e_it != e_end; ++e_it)
		{
			if (mesh_.is_flip_ok(*e_it))
			{
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