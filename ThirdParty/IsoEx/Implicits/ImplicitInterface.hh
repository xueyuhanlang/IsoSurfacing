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
#define ISOEX_IMPLICITSINTERFACE_HH

//== INCLUDES =================================================================

#include "Implicit.hh"
#include <IsoEx/Config/IsoExDefines.hh>
#include <IsoSurfGen/ImplicitFunc.h>

//== NAMESPACES ===============================================================

namespace IsoEx {
	//== CLASS DEFINITION =========================================================

	template< class Vec3 >
	class ISOEXDLLEXPORT ImplicitInterface : public Implicit<Vec3>
	{
	public:

		typedef typename Vec3::value_type real;

		/// \name Constructor & destructor
		//@{
		ImplicitInterface(ImplicitFunc<real>* func = 0)
			: m_func(func)
		{}

		/// Empty destructor
		~ImplicitInterface() {}

		//@}

		/// \name Abstract interface of implicit objects, see also IsoEx::Implicit.
		//@{
		bool is_inside(const Vec3& _point, const real _isovalue) const
		{
			return m_func->is_inside(TinyVector<real, 3>(_point[0], _point[1], _point[2]), _isovalue);
		}

		real scalar_distance(const Vec3& _point) const
		{
			return m_func->scalar_value(TinyVector<real, 3>(_point[0], _point[1], _point[2]));
		}

		bool directed_distance(const Vec3& _p0,
			const Vec3& _p1,
			Vec3& _point,
			Vec3& _normal,
			real& _distance,
			const real _isovalue = 0) const
		{
			TinyVector<real, 3> p, n;
			if (m_func->directed_distance(
				TinyVector<real, 3>(_p0[0], _p0[1], _p0[2]),
				TinyVector<real, 3>(_p1[0], _p1[1], _p1[2]),
				p, n, _distance, _isovalue)
				)
			{
				_point[0] = p[0], _point[1] = p[1], _point[2] = p[2];
				_normal[0] = n[0], _normal[1] = n[1], _normal[2] = n[2];
				return true;
			}
			return false;
		}

		void directed_distance(
			const std::vector<TinyVector<real, 3>>& _p0vec,
			const std::vector<TinyVector<real, 3>>& _p1vec,
			const std::vector<real>& _p0vec_values,
			const std::vector<real>& _p1vec_values,
			std::vector<TinyVector<real, 3>>& _point,
			std::vector<TinyVector<real, 3>>& _normal,
			const real _isovalue = 0) const
		{
			m_func->directed_distance(_p0vec, _p1vec, _p0vec_values, _p1vec_values, _point, _normal, _isovalue);
		}

		void scalar_value(
			const std::vector<TinyVector<real, 3>>& pvec,
			std::vector<real>& vec_values) const
		{
			m_func->scalar_value(pvec, vec_values);
		}
		//@}

	private:

		ImplicitFunc<real>* m_func;
	};

	//=============================================================================
} // namespace IsoEx
//=============================================================================
#endif
//=============================================================================