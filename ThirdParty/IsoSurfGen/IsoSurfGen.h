#pragma once

#include "ImplicitFunc.h"

#include <vector>
#include <unordered_map>

enum class ImpGenType { MC, MC33_CPLUS, MC33_C, EXTENDED_MC, DUAL_MC, EXTENDED_MC_GPU};

bool IsoSurfGen(
	float start_x, float end_x,
	float start_y, float end_y,
	float start_z, float end_z,
	unsigned int nGridX, unsigned int nGridY, unsigned int nGridZ,
	std::vector<float>& ptScalarField,
	ImpGenType method_type = ImpGenType::MC33_C,
	float iso_value = 0, float feature_angle = 30,
	ImplicitFunc<float>* impfunc = 0,
	const char* outobjfile = 0,
	std::vector<TinyVector<float, 3>>* vertices = 0,
	std::vector<unsigned int>* trifacets = 0, 
	const int buffer = 131072,
	bool verbose = true
);

bool IsoSurfGen(
	float start_x, float end_x,
	float start_y, float end_y,
	float start_z, float end_z,
	unsigned int nGridX, unsigned int nGridY, unsigned int nGridZ,
	ImpGenType method_type = ImpGenType::MC33_C,
	float iso_value = 0, float feature_angle = 30,
	ImplicitFunc<float>* impfunc = 0,
	const char* outobjfile = 0,
	std::vector<TinyVector<float, 3>>* vertices = 0,
	std::vector<unsigned int>* trifacets = 0,
	const int buffer = 131072,
	bool verbose = true
);

bool IsoSurfGen(
	const char* inputfile,
	ImpGenType method_type = ImpGenType::MC33_C,
	float iso_value = 0, float feature_angle = 30,
	const char* outputfile = 0,
	std::vector<TinyVector<float, 3>>* vertices = 0,
	std::vector<unsigned int>* trifacets = 0,
	const int buffer = 131072,
	bool verbose = true
);


bool IsoSurfGen_FeatureAware(
	float start_x, float end_x,
	float start_y, float end_y,
	float start_z, float end_z,
	unsigned int gridsize,
	ImplicitFunc<float>* impfunc,
	const std::vector<int>& cubes, const std::vector<int>& cubetypes,
	const std::unordered_map<int, float>& point_value_map,
	const char* outobjfile,
	float iso_value = 0,
	float feature_angle = 30,
	const int buffer = 131072,
	bool verbose = true
);