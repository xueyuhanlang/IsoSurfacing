#pragma once

#include <torch/torch.h>
#include <torch/script.h>
#include <omp.h>
#include <cmath>
#include <algorithm>
#include "ThirdParty/IsoSurfGen/ImplicitFunc.h"

template<typename Real>
class MyUtil
{
public:
	MyUtil()
	{
		max_iter = 50; threshold = 1.0e-7f;
		buffer = 131072;// 262144;
		use_GPU = true;
		dim_type = 2;
		verbose = true;
	}
public:
	int max_iter;
	int buffer;
	float threshold;
	bool use_GPU;
	int dim_type;
	at::Tensor direct_input, direct_output, direct_grad_input, batch_point_input, batch_output;
	std::vector<TinyVector<Real, 3>> world_start, world_end;
	std::vector<Real> fstart, fend, value;
	std::vector<Real> values;
	std::vector<int> point_id_total_cache, point_id_cache;
	bool verbose;
};

template <typename Real>
class MyImplicitFunc : public ImplicitFunc<Real>
{
public:
	torch::jit::script::Module* net;
	MyUtil<Real>* util;
public:
	MyImplicitFunc(torch::jit::script::Module* _net, MyUtil<Real>* _util)
		:ImplicitFunc<Real>(), net(_net), util(_util)
	{
	}
	/////////////////////////////////////////
	~MyImplicitFunc() {}
	/// //////////////////////////////////////
	bool is_inside(const TinyVector<Real, 3>& p, const Real _isovalue = 0) const
	{
		return scalar_value(p) < _isovalue;
	}
	/////////////////////////////////////////
	Real scalar_value(const TinyVector<Real, 3>& p) const
	{
		std::vector<torch::jit::IValue> inputs;
		at::Tensor m;
		if (util->dim_type == 3)
			m = torch::zeros({ 1, 1, 3 });
		else
			m = torch::zeros({ 1, 3 });

		auto mpstorage = static_cast<float*>(m.storage().data());
		mpstorage[0] = p[0], mpstorage[1] = p[1], mpstorage[2] = p[2];
		if (util->use_GPU)
			inputs.push_back(m.to(at::kCUDA));
		else
			inputs.push_back(m);

		at::Tensor output = net->forward(inputs).toTensor().cpu();
		auto outputstorage = static_cast<float*>(output.storage().data());
		return outputstorage[0];
	}
	/////////////////////////////////////////
	void gradient(const TinyVector<Real, 3>& p, TinyVector<Real, 3>& gradient) const
	{
		std::vector<torch::jit::IValue> inputs;
		at::Tensor m;
		if (util->dim_type == 3)
			m = torch::zeros({ 1, 1, 3 }, torch::requires_grad(true));
		else
			m = torch::zeros({ 1, 3 }, torch::requires_grad(true));

		auto mpstorage = static_cast<float*>(m.storage().data());
		mpstorage[0] = p[0], mpstorage[1] = p[1], mpstorage[2] = p[2];
		if (util->use_GPU)
			inputs.push_back(m.to(at::kCUDA));
		else
			inputs.push_back(m);

		at::Tensor output = net->forward(inputs).toTensor();
		output.sum().backward();
		auto gradstorage = static_cast<float*>(m.grad().cpu().storage().data());
		gradient[0] = gradstorage[0], gradient[1] = gradstorage[1], gradient[2] = gradstorage[2];
	}
	/// ////////////////////////////////////////
	void scalar_value(const std::vector<TinyVector<Real, 3>>& pvec, std::vector<Real>& vec_values) const
	{
		if (util->dim_type == 3)
			util->batch_point_input = torch::zeros({ 1, (int)pvec.size(), 3 });
		else
			util->batch_point_input = torch::zeros({ (int)pvec.size(), 3 });

		auto mpstorage = static_cast<float*>(util->batch_point_input.storage().data());
#pragma omp parallel for
		for (int i = 0; i < (int)pvec.size(); i++)
		{
			mpstorage[3 * i] = pvec[i][0], mpstorage[3 * i + 1] = pvec[i][1], mpstorage[3 * i + 2] = pvec[i][2];
		}
		std::vector<torch::jit::IValue> inputs;
		if (util->use_GPU)
			inputs.push_back(util->batch_point_input.to(at::kCUDA));
		else
			inputs.push_back(util->batch_point_input);

		util->batch_output = net->forward(inputs).toTensor().cpu();
		auto outputstorage = static_cast<float*>(util->batch_output.storage().data());
		vec_values.resize(pvec.size());
#pragma omp parallel for
		for (int i = 0; i < (int)pvec.size(); i++)
		{
			vec_values[i] = outputstorage[i];
		}
	}
	//////////////////////////////////////////////////
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

		if (util->dim_type == 3)
		{
			util->direct_input = torch::zeros({ 1, size_edges, 3 });
			util->direct_output = torch::zeros({ 1, size_edges });
			util->direct_grad_input = torch::zeros({ 1, size_edges, 3 }, torch::requires_grad(true));
		}
		else
		{
			util->direct_input = torch::zeros({ size_edges, 3 });
			util->direct_output = torch::zeros({ size_edges, 1 });
			util->direct_grad_input = torch::zeros({ size_edges, 3 }, torch::requires_grad(true));
		}

		util->world_start.assign(_p0vec.begin(), _p0vec.end());
		util->world_end.assign(_p1vec.begin(), _p1vec.end());
		util->fstart.assign(_p0vec_values.begin(), _p0vec_values.end());
		util->fend.assign(_p1vec_values.begin(), _p1vec_values.end());
		util->value.resize(size_edges);

		if (_isovalue != (Real)0.0)
		{
#pragma omp parallel for
			for (int i = 0; i < (int)util->fstart.size(); i++)
			{
				util->fstart[i] -= _isovalue;
				util->fend[i] -= _isovalue;
			}
		}

		auto mpstorage = static_cast<float*>(util->direct_input.storage().data());
		Real tiny_value = 1.0e-12;
		for (int iter = 0; iter < util->max_iter; iter++)
		{
#pragma omp parallel for
			for (int i = 0; i < size_edges; i++)
			{
				Real s0 = fabs(util->fstart[i]);
				Real s1 = fabs(util->fend[i]);
				Real t = s0 / std::max(tiny_value, s0 + s1);
				_point[i] = ((Real)1.0 - t) * util->world_start[i] + t * util->world_end[i];
				mpstorage[3 * i] = _point[i][0], mpstorage[3 * i + 1] = _point[i][1], mpstorage[3 * i + 2] = _point[i][2];
			}
			if (iter < util->max_iter - 1)
			{
				std::vector<torch::jit::IValue> inputs;
				if (util->use_GPU)
					inputs.push_back(util->direct_input.to(at::kCUDA));
				else
					inputs.push_back(util->direct_input);

				util->direct_output = net->forward(inputs).toTensor().cpu();
				auto outputstorage = static_cast<float*>(util->direct_output.storage().data());
#pragma omp parallel for
				for (int i = 0; i < size_edges; i++)
				{
					util->value[i] = outputstorage[i] - _isovalue;
					if (util->value[i] * util->fstart[i] < 0)
					{
						util->world_end[i] = _point[i];
						util->fend[i] = util->value[i];
					}
					else if (util->value[i] * util->fend[i] < 0)
					{
						util->world_start[i] = _point[i];
						util->fstart[i] = util->value[i];
					}
				}
				float vmax = 0;
				for (int i = 0; i < size_edges; i++)
					if ((_p0vec_values[i] - _isovalue) * (_p1vec_values[i] - _isovalue) < 0 && fabs(util->value[i]) > vmax)
						vmax = fabs(util->value[i]);
				if (vmax < util->threshold)
					break;
			}
		}

		auto gradstorage = static_cast<float*>(util->direct_grad_input.storage().data());
#pragma omp parallel for
		for (int i = 0; i < size_edges; i++)
		{
			gradstorage[3 * i] = _point[i][0]; gradstorage[3 * i + 1] = _point[i][1]; gradstorage[3 * i + 2] = _point[i][2];
		}
		std::vector<torch::jit::IValue> inputs;
		if (util->use_GPU)
			inputs.push_back(util->direct_grad_input.to(at::kCUDA));
		else
			inputs.push_back(util->direct_grad_input);

		at::Tensor output = net->forward(inputs).toTensor();
		output.sum().backward();
		auto gradoutputstorage = static_cast<float*>(util->direct_grad_input.grad().cpu().storage().data());

#pragma omp parallel for
		for (int i = 0; i < size_edges; i++)
		{
			_normal[i][0] = gradoutputstorage[3 * i];
			_normal[i][1] = gradoutputstorage[3 * i + 1];
			_normal[i][2] = gradoutputstorage[3 * i + 2];
			_normal[i].Normalize();
		}
	}
	//////////////////////////////////////////////////
};
