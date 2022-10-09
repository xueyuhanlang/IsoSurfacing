#pragma once

#include "TinyVector.h"
#include <vector>
#include <cmath>

template <typename Real>
class ImplicitFunc
{
public:
	ImplicitFunc()
	{}
	virtual ~ImplicitFunc() {}
	virtual bool is_inside(const TinyVector<Real, 3>& p, const Real _isovalue = 0) const
	{
		return false;
	}
	virtual Real scalar_value(const TinyVector<Real, 3>& p) const
	{
		return 0;
	}
	virtual void gradient(const TinyVector<Real, 3>& p, TinyVector<Real, 3>& gradient) const
	{
	}
	bool directed_distance(const TinyVector<Real, 3>& _p0,
		const TinyVector<Real, 3>& _p1,
		TinyVector<Real, 3>& _point,
		TinyVector<Real, 3>& _normal,
		Real& _distance,
		const Real _isovalue = 0) const
	{
		TinyVector<Real, 3> world_start = _p0, world_end = _p1;
		Real fstart = scalar_value(world_start) - _isovalue, fend = scalar_value(world_end) - _isovalue;
		Real store_f = fstart;

		//recursive bisection
		int iter = 0;
		int max_iter = 15; Real threshold = (Real)1.0e-8;
		Real value = 0;
		while (iter < max_iter)
		{
			Real s0 = fabs(fstart);
			Real s1 = fabs(fend);
			Real t = s0 / (s0 + s1);
			_point = (1.0f - t) * world_start + t * world_end;
			value = scalar_value(_point) - _isovalue;
			if (fabs(value) < threshold)
			{
				break;
			}
			else
			{
				if (value * fstart < 0)
				{
					world_end = _point;
					fend = value;
				}
				else if (value * fend < 0)
				{
					world_start = _point;
					fstart = value;
				}
				else
				{
					//the case should not happen
					break;
				}
			}
			iter++;
		}
		
		///////////////////
		gradient(_point, _normal);
		_normal.Normalize();
		_distance = store_f < 0 ? -(_point - world_start).Length() : (_point - world_start).Length();
		return true;
	}

	////////////////////////////////////
	virtual void directed_distance(
		const std::vector<TinyVector<Real, 3>>& _p0vec,
		const std::vector<TinyVector<Real, 3>>& _p1vec,
		const std::vector<Real>& _p0vec_values,
		const std::vector<Real>& _p1vec_values,
		std::vector<TinyVector<Real, 3>>& _point,
		std::vector<TinyVector<Real, 3>>& _normal,
		const Real _isovalue = 0) const
	{
	}
	////////////////////////////////////
	virtual void scalar_value(
		const std::vector<TinyVector<Real, 3>>& pvec,
		std::vector<Real>& vec_values) const
	{
	}
};