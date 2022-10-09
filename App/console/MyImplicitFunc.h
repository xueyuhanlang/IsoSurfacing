#pragma once

#include "ThirdParty/IsoSurfGen/ImplicitFunc.h"


template <typename Real>
class MyImplicitFunc : public ImplicitFunc<Real>
{
public:
	MyImplicitFunc()
		:ImplicitFunc<Real>()
	{}
	~MyImplicitFunc() {}
	bool is_inside(const TinyVector<Real, 3>& p, const Real _isovalue = 0) const
	{
		return scalar_value(p) <= _isovalue;
	}
	Real scalar_value(const TinyVector<Real, 3>& p) const
	{
		return std::max(p.SquaredLength() - (Real)1 + p[0] * p[1], (p - TinyVector<Real, 3>((Real)1, (Real)0, (Real)0)).SquaredLength() - p[1] * p[2] - (Real)1);
		//return std::max(p.SquaredLength() - 1, (p - TinyVector<Real, 3>(1, 0, 0)).SquaredLength() - 1);

		//return p.SquaredLength() - 1;
		//return p.Length() - 1;
	}
	void gradient(const TinyVector<Real, 3>& p, TinyVector<Real, 3>& gradient) const
	{
		if (p.SquaredLength() - (Real)1 + p[0] * p[1] > (p - TinyVector<Real, 3>((Real)1, (Real)0, (Real)0)).SquaredLength() - p[1] * p[2] - (Real)1)
			gradient = p + TinyVector<Real, 3>((Real)0.5 * p[1], (Real)0.5 * p[0], (Real)0);
		else
			gradient = p - TinyVector<Real, 3>((Real)1, (Real)0, (Real)0) - TinyVector<Real, 3>((Real)0, (Real)0.5 * p[2], (Real)0.5 * p[1]);
	}
	void directed_distance(
		const std::vector<TinyVector<Real, 3>>& _p0vec,
		const std::vector<TinyVector<Real, 3>>& _p1vec,
		const std::vector<Real>& _p0vec_values,
		const std::vector<Real>& _p1vec_values,
		std::vector<TinyVector<Real, 3>>& _point,
		std::vector<TinyVector<Real, 3>>& _normal,
		const Real _isovalue = 0) const
	{
		int size_edges = (int)_p0vec.size();
		_point.resize(size_edges);
		_normal.resize(size_edges);

		std::vector<TinyVector<Real, 3>> world_start, world_end;
		std::vector<Real> fstart, fend, value;
		int max_iter = 5;
		world_start.assign(_p0vec.begin(), _p0vec.end());
		world_end.assign(_p1vec.begin(), _p1vec.end());
		fstart.assign(_p0vec_values.begin(), _p0vec_values.end());
		fend.assign(_p1vec_values.begin(), _p1vec_values.end());
		value.resize(size_edges);

		if (_isovalue != (Real)0.0)
		{
#pragma omp parallel for
			for (int i = 0; i < (int)fstart.size(); i++)
			{
				fstart[i] -= _isovalue;
				fend[i] -= _isovalue;
			}
		}

		for (int iter = 0; iter < max_iter; iter++)
		{
#pragma omp parallel for
			for (int i = 0; i < size_edges; i++)
			{
				Real s0 = fabs(fstart[i]);
				Real s1 = fabs(fend[i]);
				Real t = s0 / (s0 + s1);
				_point[i] = ((Real)1.0 - t) * world_start[i] + t * world_end[i];
				value[i] = scalar_value(_point[i]) - _isovalue;
				if (value[i] * fstart[i] < 0)
				{
					world_end[i] = _point[i];
					fend[i] = value[i];
				}
				else if (value[i] * fend[i] < 0)
				{
					world_start[i] = _point[i];
					fstart[i] = value[i];
				}
			}
		}

#pragma omp parallel for
		for (int i = 0; i < size_edges; i++)
		{
			gradient(_point[i], _normal[i]);
			_normal[i].Normalize();
		}
	}
};
