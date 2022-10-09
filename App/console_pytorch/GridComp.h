#pragma once

#include "MyImplicitFunc.h"
#include "ThirdParty/IsoSurfGen/ProgressBar.h"
#include "ThirdParty/termcolor/termcolor.hpp"
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <cassert>
#include "ThirdParty/DualContour/octree.h"

int edgeTables[256] =
	{
		0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
		0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
		0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
		0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
		0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
		0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
		0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
		0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
		0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c,
		0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
		0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc,
		0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
		0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c,
		0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
		0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
		0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
		0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
		0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
		0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
		0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
		0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
		0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
		0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
		0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460,
		0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
		0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0,
		0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
		0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230,
		0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
		0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99, 0x190,
		0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
		0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0};

const int cubestatusmap[256] = {
	0, 1, 16, 17, 64, 65, 80, 81, 4, 5, 20, 21, 68, 69, 84, 85, 2, 3,
	18, 19, 66, 67, 82, 83, 6, 7, 22, 23, 70, 71, 86, 87, 32, 33, 48,
	49, 96, 97, 112, 113, 36, 37, 52, 53, 100, 101, 116, 117, 34, 35,
	50, 51, 98, 99, 114, 115, 38, 39, 54, 55, 102, 103, 118, 119, 128,
	129, 144, 145, 192, 193, 208, 209, 132, 133, 148, 149, 196, 197,
	212, 213, 130, 131, 146, 147, 194, 195, 210, 211, 134, 135, 150,
	151, 198, 199, 214, 215, 160, 161, 176, 177, 224, 225, 240, 241,
	164, 165, 180, 181, 228, 229, 244, 245, 162, 163, 178, 179, 226,
	227, 242, 243, 166, 167, 182, 183, 230, 231, 246, 247, 8, 9, 24,
	25, 72, 73, 88, 89, 12, 13, 28, 29, 76, 77, 92, 93, 10, 11, 26,
	27, 74, 75, 90, 91, 14, 15, 30, 31, 78, 79, 94, 95, 40, 41, 56,
	57, 104, 105, 120, 121, 44, 45, 60, 61, 108, 109, 124, 125, 42,
	43, 58, 59, 106, 107, 122, 123, 46, 47, 62, 63, 110, 111, 126,
	127, 136, 137, 152, 153, 200, 201, 216, 217, 140, 141, 156, 157,
	204, 205, 220, 221, 138, 139, 154, 155, 202, 203, 218, 219, 142,
	143, 158, 159, 206, 207, 222, 223, 168, 169, 184, 185, 232, 233,
	248, 249, 172, 173, 188, 189, 236, 237, 252, 253, 170, 171, 186,
	187, 234, 235, 250, 251, 174, 175, 190, 191, 238, 239, 254, 255};

class IntersectionInfo
{
public:
	IntersectionInfo(int _start, int _end)
	{
		vert[0] = _start;
		vert[1] = _end;
	}
	inline bool operator<(const IntersectionInfo &p) const
	{
		for (int i = 0; i < 2; i++)
		{
			if (vert[i] < p.vert[i])
			{
				return true;
			}
			else if (vert[i] > p.vert[i])
			{
				return false;
			}
		}
		return false;
	}
	inline bool operator==(const IntersectionInfo &p) const
	{
		for (int i = 0; i < 2; i++)
		{
			if (vert[i] != p.vert[i])
			{
				return false;
			}
		}
		return true;
	}
	inline bool operator!=(const IntersectionInfo &p) const
	{
		for (int i = 0; i < 2; i++)
		{
			if (vert[i] != p.vert[i])
			{
				return true;
			}
		}
		return false;
	}

public:
	int vert[2];
	float shift, nx, ny, nz;
};
//////////////////////////////////////
template <typename Real>
class MyOctreeLayer
{
private:
	MyImplicitFunc<Real> *myfunc;
	int depth, gridsize, offsets_[8];
	Real start_x, end_x, start_y, end_y, start_z, end_z;
	Real fCellLengthX, fCellLengthY, fCellLengthZ;
	std::vector<int> valid_cell_id;
	std::vector<int> valid_cell_types;
	std::unordered_map<int, Real> point_value_map;
	Real isovalue;
	bool verbose;

public:
	MyOctreeLayer<Real>(
		int _depth,
		MyImplicitFunc<Real> *_myfunc,
		Real _start_x, Real _end_x,
		Real _start_y, Real _end_y,
		Real _start_z, Real _end_z,
		Real _isovalue = (Real)0,
		bool _verbose = true)
		: depth(_depth), myfunc(_myfunc),
		  start_x(_start_x), end_x(_end_x),
		  start_y(_start_y), end_y(_end_y),
		  start_z(_start_z), end_z(_end_z), isovalue(_isovalue), verbose(_verbose)
	{
		assert(depth > 0);
		gridsize = pow(2, depth) + 1;

		fCellLengthX = (end_x - start_x) / (gridsize - 1);
		fCellLengthY = (end_y - start_y) / (gridsize - 1);
		fCellLengthZ = (end_z - start_z) / (gridsize - 1);

		offsets_[0] = 0;
		offsets_[1] = 1;
		offsets_[2] = 1 + gridsize;
		offsets_[3] = gridsize;
		offsets_[4] = gridsize * gridsize;
		offsets_[5] = 1 + gridsize * gridsize;
		offsets_[6] = 1 + gridsize + gridsize * gridsize;
		offsets_[7] = gridsize + gridsize * gridsize;
	}
	/////////////////////////////////////////////////
	const std::vector<int> &get_valid_cells() const { return valid_cell_id; }
	/////////////////////////////////////////////////
	const std::vector<int> &get_valid_cell_types() const { return valid_cell_types; }
	/////////////////////////////////////////////////
	const int get_grid_size() const { return gridsize; }
	/////////////////////////////////////////////////
	const int get_depth() const { return depth; }
	/////////////////////////////////////////////////
	const int get_point_id(int cell_id, int vertex_id) const
	{
		unsigned int X(gridsize - 1), Y(gridsize - 1);
		unsigned int x = cell_id % X;
		cell_id /= X;
		unsigned int y = cell_id % Y;
		cell_id /= Y;
		unsigned int z = cell_id;
		unsigned int p = x + y * gridsize + z * gridsize * gridsize;
		return p + offsets_[vertex_id];
	}
	/////////////////////////////////////////////////
	const void get_point(int globa_vertex_id, TinyVector<Real, 3> &p)
	{
		int k = globa_vertex_id / (gridsize * gridsize);
		int j = (globa_vertex_id - k * gridsize * gridsize) / gridsize;
		int i = globa_vertex_id % gridsize;
		p[2] = start_z + k * fCellLengthZ;
		p[1] = start_y + j * fCellLengthY;
		p[0] = start_x + i * fCellLengthX;
	}
	/////////////////////////////////////////////////
	const std::vector<int> &get_valid_cell_id() const { return valid_cell_id; }
	/////////////////////////////////////////////////
	const int get_gridsize() const { return gridsize; }
	/////////////////////////////////////////////////
	const std::unordered_map<int, Real> &get_point_value_map() const { return point_value_map; }
	/////////////////////////////////////////////////
	void scan_octants(MyOctreeLayer<Real> *parent = 0, bool gather_neighbors = true)
	{
		valid_cell_id.reserve(10 * gridsize * gridsize);
		valid_cell_id.resize(0);
		valid_cell_types.reserve(10 * gridsize * gridsize);
		valid_cell_types.resize(0);

		MyUtil<Real> *util = myfunc->util;

		std::vector<Real> &values = myfunc->util->values;
		std::vector<int> &point_id_total_cache = myfunc->util->point_id_total_cache;
		std::vector<int> &point_id_cache = myfunc->util->point_id_cache;

		values.resize(0);
		point_value_map.clear();
		point_id_total_cache.resize(0);
		point_id_cache.resize(0);

		////////////////////

		std::cout << termcolor::yellow << "---------"
				  << "process " << termcolor::blue << termcolor::on_white << depth << termcolor::reset << termcolor::yellow << "-th octree layer!" << std::endl;
		std::cout << termcolor::cyan;
		print_progress_bar(0.0f, verbose);

		if (parent == 0)
		{
			compute_grid_point_values(gridsize, values);

			std::vector<int> total_cell_types((gridsize - 1) * (gridsize - 1) * (gridsize - 1));

#pragma omp parallel for
			for (int c = 0; c < (gridsize - 1) * (gridsize - 1) * (gridsize - 1); c++)
			{
				unsigned char cubetype(0);
				for (int s = 0; s < 8; ++s)
					if (values[get_point_id(c, s)] > isovalue)
						cubetype |= (1 << s);
				total_cell_types[c] = cubetype;
			}
			for (int n = 0; n < (gridsize - 1) * (gridsize - 1) * (gridsize - 1); n++)
			{
				int k = n / ((gridsize - 1) * (gridsize - 1));
				int j = (n - k * (gridsize - 1) * (gridsize - 1)) / (gridsize - 1);
				int i = n % (gridsize - 1);
				int cubetype = total_cell_types[n];
				bool find = true;
				if ((cubetype == 0 || cubetype == 255))
				{
					find = false;
					if (gather_neighbors)
					{
						for (int kshift = -1; kshift <= 1; kshift++)
						{
							int kk = k + kshift;
							if (kk < 0 || kk >= gridsize - 1)
								continue;

							int cid_k = kk * (gridsize - 1) * (gridsize - 1);
							for (int jshift = -1; jshift <= 1; jshift++)
							{
								int jj = j + jshift;
								if (jj < 0 || jj >= gridsize - 1)
									continue;

								int cid_j = cid_k + jj * (gridsize - 1);

								for (int ishift = -1; ishift <= 1; ishift++)
								{
									int ii = i + ishift;
									if (ii < 0 || ii >= gridsize - 1)
										continue;
									int c_id = cid_j + ii;

									if (total_cell_types[c_id] != cubetype)
									{
										find = true;
										break;
									}
								}
								if (find)
									break;
							}
							if (find)
								break;
						}
					}
				}
				if (find)
				{
					valid_cell_types.push_back(cubetype);
					valid_cell_id.push_back(n);
				}
			}

			for (size_t i = 0; i < values.size(); i++)
				point_value_map[i] = values[i];
		}
		else
		{
			assert(parent->get_depth() == depth - 1);
			for (size_t c = 0; c < parent->get_valid_cell_id().size(); c++)
			{
				int n = parent->get_valid_cell_id()[c];
				int k = n / ((parent->get_gridsize() - 1) * (parent->get_gridsize() - 1));
				int j = (n - k * (parent->get_gridsize() - 1) * (parent->get_gridsize() - 1)) / (parent->get_gridsize() - 1);
				int i = n % (parent->get_gridsize() - 1);
				k *= 2, j *= 2, i *= 2;
				int mid[8];
				mid[0] = k * (gridsize - 1) * (gridsize - 1) + j * (gridsize - 1) + i;
				mid[1] = k * (gridsize - 1) * (gridsize - 1) + j * (gridsize - 1) + i + 1;
				mid[2] = k * (gridsize - 1) * (gridsize - 1) + (j + 1) * (gridsize - 1) + i;
				mid[3] = k * (gridsize - 1) * (gridsize - 1) + (j + 1) * (gridsize - 1) + i + 1;
				mid[4] = (k + 1) * (gridsize - 1) * (gridsize - 1) + j * (gridsize - 1) + i;
				mid[5] = (k + 1) * (gridsize - 1) * (gridsize - 1) + j * (gridsize - 1) + i + 1;
				mid[6] = (k + 1) * (gridsize - 1) * (gridsize - 1) + (j + 1) * (gridsize - 1) + i;
				mid[7] = (k + 1) * (gridsize - 1) * (gridsize - 1) + (j + 1) * (gridsize - 1) + i + 1;

				for (int s = 0; s < 8; s++)
				{
					for (int v = 0; v < 8; v++)
						point_value_map[get_point_id(mid[s], v)] = 1;
				}
			}

			typename std::unordered_map<int, Real>::iterator iter = point_value_map.begin();
			point_id_total_cache.resize(point_value_map.size());
			for (int i = 0; iter != point_value_map.end(); iter++, i++)
				point_id_total_cache[i] = iter->first;
			at::Tensor m;
			if (util->dim_type == 3)
				m = torch::zeros({1, std::min((int)point_id_total_cache.size(), util->buffer), 3});
			else
				m = torch::zeros({std::min((int)point_id_total_cache.size(), util->buffer), 3});

			auto mpstorage = static_cast<float *>(m.storage().data());
			at::Tensor output;

			int counter = 0;
			for (size_t pi = 0; pi < point_id_total_cache.size(); pi++)
			{
				int n = point_id_total_cache[pi];
				int k = n / (gridsize * gridsize);
				int j = (n - k * gridsize * gridsize) / gridsize;
				int i = n % gridsize;
				mpstorage[counter + 2] = start_z + k * fCellLengthZ;
				mpstorage[counter + 1] = start_y + j * fCellLengthY;
				mpstorage[counter] = start_x + i * fCellLengthX;

				point_id_cache.push_back(n);
				counter += 3;

				if (counter == 3 * util->buffer || pi + 1 == point_id_total_cache.size())
				{
					std::vector<torch::jit::IValue> inputs;
					if (util->use_GPU)
					{
						inputs.push_back(m.to(at::kCUDA));
					}
					else
						inputs.push_back(m);

					output = myfunc->net->forward(inputs).toTensor().cpu();
					auto outputstorage = static_cast<float *>(output.storage().data());

					for (size_t s = 0; s < point_id_cache.size(); s++)
					{
						point_value_map[point_id_cache[s]] = outputstorage[s];
					}

					print_progress_bar((pi + 1.0f) / point_id_total_cache.size(), verbose);
					if (pi + 1 == point_id_total_cache.size())
						printf("\n");

					point_id_cache.resize(0);
					counter = 0;
				}
			}

			valid_cell_types.reserve(8 * parent->get_valid_cell_id().size());
			valid_cell_types.resize(0);

			// create its own validcells
			std::vector<int> total_cell_types(8 * parent->get_valid_cell_id().size());
			std::vector<int> total_cells(8 * parent->get_valid_cell_id().size());
#pragma omp parallel for
			for (int c = 0; c < (int)parent->get_valid_cell_id().size(); c++)
			{
				int n = parent->get_valid_cell_id()[c];
				int k = n / ((parent->get_gridsize() - 1) * (parent->get_gridsize() - 1));
				int j = (n - k * (parent->get_gridsize() - 1) * (parent->get_gridsize() - 1)) / (parent->get_gridsize() - 1);
				int i = n % (parent->get_gridsize() - 1);
				k *= 2, j *= 2, i *= 2;
				int *mid = &total_cells[8 * c];
				mid[0] = k * (gridsize - 1) * (gridsize - 1) + j * (gridsize - 1) + i;
				mid[1] = k * (gridsize - 1) * (gridsize - 1) + j * (gridsize - 1) + i + 1;
				mid[2] = k * (gridsize - 1) * (gridsize - 1) + (j + 1) * (gridsize - 1) + i;
				mid[3] = k * (gridsize - 1) * (gridsize - 1) + (j + 1) * (gridsize - 1) + i + 1;
				mid[4] = (k + 1) * (gridsize - 1) * (gridsize - 1) + j * (gridsize - 1) + i;
				mid[5] = (k + 1) * (gridsize - 1) * (gridsize - 1) + j * (gridsize - 1) + i + 1;
				mid[6] = (k + 1) * (gridsize - 1) * (gridsize - 1) + (j + 1) * (gridsize - 1) + i;
				mid[7] = (k + 1) * (gridsize - 1) * (gridsize - 1) + (j + 1) * (gridsize - 1) + i + 1;

				for (int ci = 0; ci < 8; ci++)
				{
					unsigned char cubetype(0);
					for (int s = 0; s < 8; ++s)
						if (point_value_map.find(get_point_id(mid[ci], s))->second > isovalue)
							cubetype |= (1 << s);

					total_cell_types[8 * c + ci] = cubetype;
				}
			}

			std::unordered_map<int, int> cell_map;
			for (size_t c = 0; c < total_cells.size(); c++)
			{
				cell_map[total_cells[c]] = c;
			}

			for (size_t c = 0; c < total_cells.size(); c++)
			{
				int n = total_cells[c];
				int k = n / ((gridsize - 1) * (gridsize - 1));
				int j = (n - k * (gridsize - 1) * (gridsize - 1)) / (gridsize - 1);
				int i = n % (gridsize - 1);
				int cubetype = total_cell_types[c];
				bool find = true;
				if ((cubetype == 0 || cubetype == 255))
				{
					find = false;
					if (gather_neighbors)
					{
						for (int kshift = -1; kshift <= 1; kshift++)
						{
							int kk = k + kshift;
							if (kk < 0 || kk >= gridsize - 1)
								continue;

							int cid_k = kk * (gridsize - 1) * (gridsize - 1);
							for (int jshift = -1; jshift <= 1; jshift++)
							{
								int jj = j + jshift;
								if (jj < 0 || jj >= gridsize - 1)
									continue;

								int cid_j = cid_k + jj * (gridsize - 1);

								for (int ishift = -1; ishift <= 1; ishift++)
								{
									int ii = i + ishift;
									if (ii < 0 || ii >= gridsize - 1)
										continue;
									int c_id = cid_j + ii;

									typename std::unordered_map<int, int>::iterator miter = cell_map.find(c_id);

									if (miter != cell_map.end())
									{
										if (total_cell_types[miter->second] != cubetype)
										{
											find = true;
											break;
										}
									}
								}
								if (find)
									break;
							}
							if (find)
								break;
						}
					}
				}
				if (find)
				{
					valid_cell_types.push_back(cubetype);
					valid_cell_id.push_back(n);
				}
			}
		}
		// std::cout << termcolor::reset << "valid_cells: " << valid_cell_id.size() << std::endl;
	}
	/////////////////////////////////////////////////
	void compute_grid_point_values(int specified_size, std::vector<Real> &values)
	{
		MyUtil<Real> *util = myfunc->util;

		values.resize(specified_size * specified_size * specified_size);

		Real _fCellLengthX = (end_x - start_x) / (specified_size - 1);
		Real _fCellLengthY = (end_y - start_y) / (specified_size - 1);
		Real _fCellLengthZ = (end_z - start_z) / (specified_size - 1);
		at::Tensor m;
		if (util->dim_type == 3)
			m = torch::zeros({1, std::min((int)values.size(), util->buffer), 3});
		else // USE_2_DIM_NETWORK
			m = torch::zeros({std::min((int)values.size(), util->buffer), 3});

		at::Tensor output;
		auto mpstorage = static_cast<float *>(m.storage().data());
		int counter = 0;
		int start_n = 0;
		for (int n = 0; n < specified_size * specified_size * specified_size; n++)
		{
			int k = n / (specified_size * specified_size);
			int j = (n - k * specified_size * specified_size) / specified_size;
			int i = n % specified_size;
			mpstorage[3 * counter + 2] = start_z + k * _fCellLengthZ;
			mpstorage[3 * counter + 1] = start_y + j * _fCellLengthY;
			mpstorage[3 * counter] = start_x + i * _fCellLengthX;
			counter++;
			if (counter == util->buffer || n + 1 == specified_size * specified_size * specified_size)
			{
				std::vector<torch::jit::IValue> inputs;
				if (util->use_GPU)
					inputs.push_back(m.to(at::kCUDA));
				else
					inputs.push_back(m);

				output = myfunc->net->forward(inputs).toTensor().cpu();
				auto outputstorage = static_cast<float *>(output.storage().data());

#pragma omp parallel for
				for (int l = start_n; l <= n; l++)
					values[l] = outputstorage[l - start_n];

				start_n = n + 1;
				counter = 0;

				print_progress_bar((n + 1.0f) / (specified_size * specified_size * specified_size), verbose);
				if (n + 1 == specified_size * specified_size * specified_size)
					printf("\n");
			}
		}
	}
	/////////////////////////////////////////////////
	void save_cubes(const char filename[]) // *.obj format
	{
		std::ofstream mfile(filename);

		int cubfacets[][4] = {{3, 2, 1, 0}, {6, 5, 1, 2}, {4, 0, 1, 5}, {4, 7, 3, 0}, {4, 5, 6, 7}, {2, 3, 7, 6}};
		int counter = 1;
		for (size_t c = 0; c < valid_cell_id.size(); c++)
		{
			int n = valid_cell_id[c];
			int k = n / ((gridsize - 1) * (gridsize - 1));
			int j = (n - k * (gridsize - 1) * (gridsize - 1)) / (gridsize - 1);
			int i = n % (gridsize - 1);

			Real x = start_x + i * fCellLengthX;
			Real y = start_y + j * fCellLengthY;
			Real z = start_z + k * fCellLengthZ;

			mfile << "v " << x << ' ' << y << ' ' << z << std::endl;
			mfile << "v " << x + fCellLengthX << ' ' << y << ' ' << z << std::endl;
			mfile << "v " << x + fCellLengthX << ' ' << y + fCellLengthY << ' ' << z << std::endl;
			mfile << "v " << x << ' ' << y + fCellLengthY << ' ' << z << std::endl;
			mfile << "v " << x << ' ' << y << ' ' << z + fCellLengthZ << std::endl;
			mfile << "v " << x + fCellLengthX << ' ' << y << ' ' << z + fCellLengthZ << std::endl;
			mfile << "v " << x + fCellLengthX << ' ' << y + fCellLengthY << ' ' << z + fCellLengthZ << std::endl;
			mfile << "v " << x << ' ' << y + fCellLengthY << ' ' << z + fCellLengthZ << std::endl;
			for (int i = 0; i < 6; i++)
				mfile << "f "
					  << counter + cubfacets[i][0] << ' '
					  << counter + cubfacets[i][1] << ' '
					  << counter + cubfacets[i][2] << ' '
					  << counter + cubfacets[i][3] << std::endl;
			counter += 8;
		}
		mfile.close();
	}
	/////////////////////////////////////////////////
	void compute_intersection(std::set<IntersectionInfo> &edgeintersections)
	{
		edgeintersections.clear();
		////////////////////////////////////////////////////////////////////
		std::vector<int> cache_cube(valid_cell_id.size() * 2);
		for (size_t i = 0; i < valid_cell_id.size(); i++)
		{
			cache_cube[2 * i] = valid_cell_id[i];
			cache_cube[2 * i + 1] = valid_cell_types[i];
		}
		MyUtil<Real> *util = myfunc->util;
		/////////////////////////////////////////////////////////////////////
		int buffer = myfunc->util->buffer;
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

		TinyVector<Real, 3> ipoint[8];
		int label[12] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
		int edgepair[12][2] = {{0, 1}, {1, 2}, {3, 2}, {0, 3}, {4, 5}, {5, 6}, {7, 6}, {4, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
		int corner[8];
		Real values[8];
		std::set<std::pair<int, int>> edgeset;
		std::vector<int> edgestore;
		edgestore.reserve(2 * (buffer + 12));
		Real unit_length = (end_z - start_z) / (gridsize - 1);

		for (int s = 0; s < cache_cube.size(); s += 2)
		{
			int c = cache_cube[s];
			for (int i = 0; i < 8; ++i)
			{
				corner[i] = get_point_id(c, i);
				values[i] = point_value_map.find(corner[i])->second;
			}

			unsigned char cubetype = cache_cube[s + 1];

			for (int i = 0; i < 8; i++)
			{
				get_point(corner[i], ipoint[i]);
			}

			for (int i = 0; i < 12; i++)
			{
				if (edgeTables[cubetype] & label[i])
				{
					std::pair<int, int> mpair(corner[edgepair[i][0]], corner[edgepair[i][1]]);
					if (edgeset.find(mpair) != edgeset.end())
					{
						;
					}
					else
					{
						edgeset.insert(mpair);
						edgestore.push_back(corner[edgepair[i][0]]);
						edgestore.push_back(corner[edgepair[i][1]]);

						edge_start.push_back(ipoint[edgepair[i][0]]);
						edge_end.push_back(ipoint[edgepair[i][1]]);
						start_fvalue.push_back(values[edgepair[i][0]]);
						end_fvalue.push_back(values[edgepair[i][1]]);
						counter++;
					}
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

				myfunc->directed_distance(edge_start, edge_end, start_fvalue, end_fvalue, intersection, normals, isovalue);

				for (size_t h = 0; h < intersection.size(); h++)
				{
					IntersectionInfo info(edgestore[2 * h], edgestore[2 * h + 1]);
					info.nx = (float)normals[h][0], info.ny = (float)normals[h][1], info.nz = (float)normals[h][2];
					info.shift = (float)((intersection[h] - edge_start[h]).Length() / unit_length);
					edgeintersections.insert(info);
				}

				counter = 0;

				edgestore.resize(0);
				edgeset.clear();
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
	}
	/////////////////////////////////////////////////
	void save_dcf(const char filename[], bool binary) // modified dcf format
	{
		if (valid_cell_id.empty())
		{
			std::cerr << "scan_octants() should be called first!" << std::endl;
			return;
		}
		//--- compute intersection---
		std::set<IntersectionInfo> edgeintersections;
		compute_intersection(edgeintersections);
		//---------------------------

		class Node
		{
		public:
			Node(int _id, int _depth, bool _is_leaf = false, int _cell_type = 0)
				: id_(_id), depth_(_depth), is_leaf_(_is_leaf)
			{
				parent = 0;
				cell_type_ = cubestatusmap[_cell_type];
				for (int i = 0; i < 8; i++)
					child[i] = 0;
				for (int i = 0; i < 12; i++)
					interinfo[i] = 0;
			}
			~Node()
			{
				for (int i = 0; i < 12; i++)
					if (interinfo[i])
						delete interinfo[i];
			}
			void travel_print(std::ofstream &out, bool binary)
			{
				bool is_empty = true;
				for (int i = 0; i < 8; i++)
				{
					if (child[i])
					{
						is_empty = false;
						break;
					}
				}

				int type = is_leaf_ ? 2 : (is_empty ? 1 : 0);
				if (binary)
					out.write((char *)&type, sizeof(int));
				else
				{
					for (int k = 0; k < depth_; k++)
						out << "*";
					out << "type: " << type << std::endl;
				}

				if (type == 0)
				{
					for (int i = 0; i < 8; i++)
					{
						if (child[i])
							child[i]->travel_print(out, binary);
						else
						{
							if (binary)
							{
								int ctype = 1;
								out.write((char *)&ctype, sizeof(int));
							}
							else
							{
								out << "type: 1" << std::endl;
							}
						}
					}
				}
				else if (type == 1)
				{
					// short sg = 1; // it is not a true sign value, faked for DCF.
					// out.write((char*)&sg, sizeof(short));
				}
				else
				{
					if (binary)
						out.write((char *)&cell_type_, sizeof(int)); // different with orginal dcf format
					else
					{
						for (int k = 0; k < depth_; k++)
							out << "*";
						out << "c" << cell_type_ << std::endl;
					}
					const int zero = 0, one = 1;
					for (int l = 0; l < 12; l++)
					{
						if (binary)
						{
							if (interinfo[l])
							{
								out.write((char *)&one, sizeof(int));
								if (fabs((*interinfo[l])[0] - 0.524244) < 0.00001)
								{
									int k = 0;
									k = 1;
								}
								float h = (*interinfo[l])[0], x = (*interinfo[l])[1], y = (*interinfo[l])[2], z = (*interinfo[l])[3];
								// out.write((char*)&(*interinfo[l]), sizeof(float)*4);
								out.write((char *)&h, sizeof(float));
								out.write((char *)&x, sizeof(float));
								out.write((char *)&y, sizeof(float));
								out.write((char *)&z, sizeof(float));
							}
							else
								out.write((char *)&zero, sizeof(int));
						}
						else
						{
							for (int k = 0; k < depth_; k++)
								out << "*";
							out << "i";
							if (interinfo[l])
							{
								out << 1 << ' ' << *interinfo[l] << std::endl;
							}
							else
								out << 0 << std::endl;
						}
					}
				}
			}

		public:
			int id_, depth_, cell_type_;
			bool is_leaf_;
			Node *parent;
			Node *child[8];
			TinyVector<float, 4> *interinfo[12];
		};
		//--------------------

		std::vector<std::unordered_map<int, Node *>> layernodes(depth + 1);
		// intialize
		const int vertmap2dc[8] = {0, 4, 3, 7, 1, 5, 2, 6};
		const int edgevmap[12][2] =
			{{0, 4}, {1, 5}, {2, 6}, {3, 7}, {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 1}, {2, 3}, {4, 5}, {6, 7}};
		int corner[8];
		for (size_t c = 0; c < valid_cell_id.size(); c++)
		{
			int n = valid_cell_id[c];
			Node *node = new Node(n, depth, true, valid_cell_types[c]);
			layernodes[depth][n] = node;
			// todo add 12 edges here
			for (int i = 0; i < 8; ++i)
			{
				corner[i] = get_point_id(n, i);
			}
			for (int i = 0; i < 12; i++)
			{
				auto iter = edgeintersections.find(
					IntersectionInfo(corner[vertmap2dc[edgevmap[i][0]]], corner[vertmap2dc[edgevmap[i][1]]]));
				if (iter != edgeintersections.end())
				{
					node->interinfo[i] = new TinyVector<float, 4>(iter->shift, iter->nx, iter->ny, iter->nz);
				}
			}
		}

		int curdepth = depth;
		int res = gridsize - 1;

		while (curdepth > 0)
		{
			int half_res = res / 2;
			auto iter = layernodes[curdepth].begin();
			for (; iter != layernodes[curdepth].end(); iter++)
			{
				int n = iter->first;
				int k = n / (res * res);
				int j = (n - k * res * res) / res;
				int i = n % res;
				int inck = k % 2, incj = j % 2, inci = i % 2;
				k = k / 2, j = j / 2, i = i / 2;

				int up_n = k * half_res * half_res + j * half_res + i;
				Node *node = 0;
				auto siter = layernodes[curdepth - 1].find(up_n);
				if (siter == layernodes[curdepth - 1].end())
				{
					node = new Node(up_n, curdepth - 1);
					layernodes[curdepth - 1][up_n] = node;
				}
				else
				{
					node = siter->second;
				}
				assert(node->child[inci * 4 + incj * 2 + inck] == 0);
				node->child[inci * 4 + incj * 2 + inck] = iter->second;
			}
			curdepth--;
			res /= 2;
		}

		///--------------------
		std::ofstream dcffile;
		if (binary)
			dcffile.open(filename, std::ofstream::binary);
		else
			dcffile.open(filename);
		const char version[10] = "multisign";
		int dim[3];
		dim[0] = dim[1] = dim[2] = gridsize - 1;
		if (binary)
		{
			dcffile.write(version, sizeof(char) * 10);
			dcffile.write((char *)dim, sizeof(int) * 3);
		}
		else
		{
			dcffile << version << std::endl;
			dcffile << dim[0] << ' ' << dim[1] << ' ' << dim[2] << std::endl;
		}

		layernodes[0][0]->travel_print(dcffile, binary);
		dcffile.close();

		//--------------------

		for (int i = 0; i < depth + 1; i++)
		{
			for (auto iter = layernodes[i].begin(); iter != layernodes[i].end(); iter++)
				delete iter->second;
		}
	}
	/////////////////////////////////////////////////
	void dcf_meshing(const char plyname[])
	{
		if (valid_cell_id.empty())
		{
			std::cerr << "scan_octants() should be called first!" << std::endl;
			return;
		}
		//--- compute intersection---
		std::set<IntersectionInfo> edgeintersections;
		compute_intersection(edgeintersections);
		//---------------------------
		class Node
		{
		public:
			Node(int _id, int _depth, bool _is_leaf = false, int _cell_type = 0)
				: id_(_id), depth_(_depth), is_leaf_(_is_leaf)
			{
				parent = 0;
				cell_type_ = cubestatusmap[_cell_type];
				for (int i = 0; i < 8; i++)
					child[i] = 0;
				for (int i = 0; i < 12; i++)
					interinfo[i] = 0;
			}
			~Node()
			{
				for (int i = 0; i < 12; i++)
					if (interinfo[i])
						delete interinfo[i];
			}
			OctreeNode *travel(int *st, int len, int ht)
			{
				OctreeNode *rvalue = 0;
				bool is_empty = true;
				for (int i = 0; i < 8; i++)
				{
					if (child[i])
					{
						is_empty = false;
						break;
					}
				}

				int type = is_leaf_ ? 2 : (is_empty ? 1 : 0);

				if (type == 0)
				{
					rvalue = new InternalNode();
					int nlen = len / 2;
					int nst[3];

					for (int i = 0; i < 8; i++)
					{
						if (child[i])
						{
							nst[0] = st[0] + vertMap[i][0] * nlen;
							nst[1] = st[1] + vertMap[i][1] * nlen;
							nst[2] = st[2] + vertMap[i][2] * nlen;
							((InternalNode *)rvalue)->child[i] = child[i]->travel(nst, nlen, ht - 1);
						}
					}
				}
				else if (type == 1)
				{
					return 0;
				}
				else if (type == 2)
				{
					float inters[12][3], norms[12][3];
					int numinters = 0;

					for (int l = 0; l < 12; l++)
					{
						if (interinfo[l])
						{
							int dir = l / 4;
							int base = edgevmap[l][0];
							inters[numinters][0] = st[0] + vertMap[base][0] * len;
							inters[numinters][1] = st[1] + vertMap[base][1] * len;
							inters[numinters][2] = st[2] + vertMap[base][2] * len;
							inters[numinters][dir] += (*interinfo[l])[0];

							norms[numinters][0] = (*interinfo[l])[1];
							norms[numinters][1] = (*interinfo[l])[2];
							norms[numinters][2] = (*interinfo[l])[3];

							numinters++;
						}
					}

					if (numinters > 0)
						rvalue = new LeafNode(ht, (unsigned char)cell_type_, st, len, numinters, inters, norms);
				}

				return rvalue;
			}

		public:
			int id_, depth_, cell_type_;
			bool is_leaf_;
			Node *parent;
			Node *child[8];
			TinyVector<float, 4> *interinfo[12];
		};
		//--------------------

		std::vector<std::unordered_map<int, Node *>> layernodes(depth + 1);
		// intialize
		const int vertmap2dc[8] = {0, 4, 3, 7, 1, 5, 2, 6};
		const int edgevmap[12][2] =
			{{0, 4}, {1, 5}, {2, 6}, {3, 7}, {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 1}, {2, 3}, {4, 5}, {6, 7}};
		int corner[8];
		for (size_t c = 0; c < valid_cell_id.size(); c++)
		{
			int n = valid_cell_id[c];
			Node *node = new Node(n, depth, true, valid_cell_types[c]);
			layernodes[depth][n] = node;
			// todo add 12 edges here
			for (int i = 0; i < 8; ++i)
			{
				corner[i] = get_point_id(n, i);
			}
			for (int i = 0; i < 12; i++)
			{
				auto iter = edgeintersections.find(
					IntersectionInfo(corner[vertmap2dc[edgevmap[i][0]]], corner[vertmap2dc[edgevmap[i][1]]]));
				if (iter != edgeintersections.end())
				{
					node->interinfo[i] = new TinyVector<float, 4>(iter->shift, iter->nx, iter->ny, iter->nz);
				}
			}
		}

		int curdepth = depth;
		int res = gridsize - 1;

		while (curdepth > 0)
		{
			int half_res = res / 2;
			auto iter = layernodes[curdepth].begin();
			for (; iter != layernodes[curdepth].end(); iter++)
			{
				int n = iter->first;
				int k = n / (res * res);
				int j = (n - k * res * res) / res;
				int i = n % res;
				int inck = k % 2, incj = j % 2, inci = i % 2;
				k = k / 2, j = j / 2, i = i / 2;

				int up_n = k * half_res * half_res + j * half_res + i;
				Node *node = 0;
				auto siter = layernodes[curdepth - 1].find(up_n);
				if (siter == layernodes[curdepth - 1].end())
				{
					node = new Node(up_n, curdepth - 1);
					layernodes[curdepth - 1][up_n] = node;
				}
				else
				{
					node = siter->second;
				}
				assert(node->child[inci * 4 + incj * 2 + inck] == 0);
				node->child[inci * 4 + incj * 2 + inck] = iter->second;
			}
			curdepth--;
			res /= 2;
		}

		///--------------------
		Octree myOctree;
		myOctree.dimen = gridsize - 1;
		myOctree.hasQEF = 0;
		myOctree.maxDepth = depth;
		myOctree.scale = std::max(std::max(end_x - start_x, end_y - start_y), end_z - start_z) / pow(2, depth);
		myOctree.shift[0] = end_x;
		myOctree.shift[1] = end_y;
		myOctree.shift[2] = end_z;
		int st[3] = {0, 0, 0};
		int len = myOctree.dimen;
		int ht = myOctree.maxDepth;
		myOctree.root = layernodes[0][0]->travel(st, len, ht);
		// myOctree.genContour(plyname);
		myOctree.genContourNoInter2(plyname);
		//--------------------
		Node *root = layernodes[0].begin()->second;
		delete root;
		// for (int i = 0; i < depth + 1; i++)
		//{
		//	for (auto iter = layernodes[i].begin(); iter != layernodes[i].end(); iter++)
		//		delete iter->second;
		// }
	}
	/////////////////////////////////////////////////
	void save_vtk(const char filename[])
	{
		std::ofstream mfile(filename);

		mfile << "# vtk DataFile Version 3.0\nImplicit Function\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		mfile << "POINTS " << point_value_map.size() << " float\n";

		std::map<int, int> idmap;
		TinyVector<Real, 3> pos;
		int counter = 0;
		for (auto iter = point_value_map.begin(); iter != point_value_map.end(); iter++)
		{
			int id = iter->first;
			get_point(iter->first, pos);
			mfile << pos << std::endl;
			idmap[id] = counter++;
		}
		mfile << "CELLS " << valid_cell_id.size() << ' ' << valid_cell_id.size() * 9 << std::endl;
		for (size_t c = 0; c < valid_cell_id.size(); c++)
		{
			mfile << 8;

			int n = valid_cell_id[c];
			int k = n / ((gridsize - 1) * (gridsize - 1));
			int j = (n - k * (gridsize - 1) * (gridsize - 1)) / (gridsize - 1);
			int i = n % (gridsize - 1);
			n = k * gridsize * gridsize + j * gridsize + i;
			mfile << ' ' << idmap[n];
			mfile << ' ' << idmap[n + offsets_[1]];
			mfile << ' ' << idmap[n + offsets_[2]];
			mfile << ' ' << idmap[n + offsets_[3]];
			mfile << ' ' << idmap[n + offsets_[4]];
			mfile << ' ' << idmap[n + offsets_[5]];
			mfile << ' ' << idmap[n + offsets_[6]];
			mfile << ' ' << idmap[n + offsets_[7]] << std::endl;
		}

		mfile << "CELL_TYPES " << valid_cell_id.size() << std::endl;
		for (size_t c = 0; c < valid_cell_id.size(); c++)
			mfile << 12 << std::endl;

		mfile << "POINT_DATA " << point_value_map.size() << std::endl;
		mfile << "SCALARS funcvalue float 1\nLOOKUP_TABLE default\n";
		for (auto iter = point_value_map.begin(); iter != point_value_map.end(); iter++)
		{
			mfile << iter->second << std::endl;
		}

		mfile.close();
	}
};
//////////////////////////////////////
