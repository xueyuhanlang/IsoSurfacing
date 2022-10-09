#include "ThirdParty/cxxopts/cxxopts.hpp"
#include "ThirdParty/IsoSurfGen/IsoSurfGen.h"
#include <chrono>
#include "GridComp.h"
#include "ThirdParty/Happly/happly.h"
#include "ThirdParty/tmd/TriangleMeshDistance.h"
#include "mesh_vis.h"

//////////////////////////////////////
int main(int argc, char** argv)
{
	float isovalue = 0, feature_angle = 30;
	int depth = 7;
	int gridsize = pow(2, depth) + 1;
	float start_x = -1, end_x = 1, start_y = -1, end_y = 1, start_z = -1, end_z = 1;
	bool dump_all_layer_mesh = false, output_vtk = false;
	bool use_octree = false;
	MyUtil<float> util;
	try
	{
		cxxopts::Options options("ISG_pytorch", "IsoSurfaceGenator (author: Yang Liu, Email: yangliu@microsoft.com)");
		options
			.positional_help("[optional args]")
			.show_positional_help()
			.allow_unrecognised_options()
			.add_options()
			("i,input", "input model(*.pt)", cxxopts::value<std::string>())
			("d,depth", "octree depth (2-10)(default: 7)", cxxopts::value<int>())
			("b,box", "bounding box size (default: 2)", cxxopts::value<float>())
			("a,angle", "feature angle threshold for EMC  (default: 30 degree)", cxxopts::value<float>())
			("m,method", "Method: EMC, DC, MC (default: EMC)", cxxopts::value<std::string>())
			//("o,output", "output mesh (obj/ply format) or vtk format (for volumetric visualization)", cxxopts::value<std::string>())
			("o,output", "output mesh (ply format)", cxxopts::value<std::string>())
			("t,threshold", "threshold value for computing intersection. (default: 1e-7)", cxxopts::value<float>())
			("n,maxiter", "max iteration number for computing intersection. (default: 50)", cxxopts::value<int>())
			("g,gpu", "use GPU model (default: true)", cxxopts::value<bool>())
			("s,setbuf", "set buffer size (default: 131072)", cxxopts::value<int>())
			("x,type", "set model type (2dim or 3dim) (default: 2)", cxxopts::value<int>())
			("l,alldepth", "generate meshes for each depth layer (default, false)", cxxopts::value<bool>())
			("v,isovalue", "isovalue (default: 0)", cxxopts::value<float>())
			("k,vtkoput", "output vtk format", cxxopts::value<bool>())
			("c,useoctree", "use octree speedup (default, false)", cxxopts::value<bool>())
			("y,verbose", "print progress (default, true)", cxxopts::value<bool>())
			("h,help", "Print help")
			("compare", "GT mesh for SDF comparison (*.ply)", cxxopts::value<std::string>())
			("compres", "resolution for for SDF comparison (default: 32)", cxxopts::value<int>());
		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}
		std::string inputfilename, outputfilename, configname;
		if (result.count("y"))
		{
			util.verbose = result["y"].as<bool>();
		}
		if (result.count("g"))
		{
			util.use_GPU = result["g"].as<bool>();
		}
		if (result.count("s"))
		{
			util.buffer = std::max(1, result["s"].as<int>());
		}
		if (result.count("x"))
		{
			util.dim_type = result["x"].as<int>();
			if (util.dim_type != 2 && util.dim_type != 3)
			{
				std::cerr << "type should be 2 or 3!" << std::endl;
				return -1;
			}
		}
		if (result.count("v"))
		{
			isovalue = result["v"].as<float>();
		}

		torch::jit::script::Module module;
		if (result.count("i"))
		{
			inputfilename = result["i"].as<std::string>();
			try
			{
				if (util.use_GPU)
					module = torch::jit::load(inputfilename, torch::kCUDA);
				else
					module = torch::jit::load(inputfilename);
			}
			catch (const c10::Error& e)
			{
				std::cerr << termcolor::red << "error: loading the model.\n" << termcolor::reset;
				return -1;
			}
		}
		else
			throw 1;

		if (result.count("n"))
		{
			util.max_iter = std::max(1, result["n"].as<int>());
		}

		if (result.count("t"))
		{
			util.threshold = std::max(0.0f, result["t"].as<float>());
		}

		if (result.count("r"))
		{
			int r = std::max(2, result["r"].as<int>());
			gridsize = 1 + pow(2, (int)ceil(log2f((float)r)));
		}

		if (result.count("d"))
		{
			depth = std::max(1, result["d"].as<int>());
			gridsize = 1 + pow(2, depth);
		}

		if (result.count("b"))
		{
			float boxsize = std::max(result["b"].as<float>(), 0.001f);
			start_x = -boxsize / 2, end_x = boxsize / 2, start_y = -boxsize / 2, end_y = boxsize / 2, start_z = -boxsize / 2, end_z = boxsize / 2;
		}

		if (result.count("a"))
			feature_angle = std::max(0.0f, result["a"].as<float>());

		if (result.count("l"))
			dump_all_layer_mesh = result["l"].as<bool>();
		if (result.count("k"))
			output_vtk = result["k"].as<bool>();
		if (result.count("c"))
		{
			use_octree = result["c"].as<bool>();
		}
		bool use_DC = false, use_MC = false;
		if (result.count("m"))
		{
			auto method = result["m"].as<std::string>();
			if (method == "DC")
			{
				use_DC = true;
			}
			else if (method == "MC")
				use_MC = true;

			if (use_MC && use_octree)
			{
				std::cout << "Please set --useoctree=false for using -m use_MC !\n";
				return -1;
			}
		}

		if (result.count("compare") && result.count("i"))
		{
			MyImplicitFunc<float> impfunc(&module, &util);
			std::string gtmesh = result["compare"].as<std::string>();
			int res = 64;
			if (result.count("compres"))
			{
				res = result["compres"].as<int>();
			}
			util.point_id_total_cache.reserve(10 * res * res);
			util.point_id_cache.reserve(util.buffer);
			std::vector<float> ptScalarField(res * res * res), gtField(res * res * res);
			MyOctreeLayer<float> mygrid(1, &impfunc, start_x, end_x, start_y, end_y, start_z, end_z, 0, false);
			mygrid.compute_grid_point_values(res, ptScalarField);

			happly::PLYData plyIn(gtmesh);
			std::vector<std::array<double, 3>> vertices = plyIn.getVertexPositions();
			std::vector<std::vector<int>> facets = plyIn.getFaceIndices<int>();
			std::vector<std::array<int, 3>> triangles(facets.size());
			for (size_t i = 0; i < facets.size(); i++)
			{
				triangles[i][0] = facets[i][0], triangles[i][1] = facets[i][1], triangles[i][2] = facets[i][2];
			}
			tmd::TriangleMeshDistance mesh_distance(vertices, triangles);
			double _fCellLengthX = (end_x - start_x) / (res - 1);
			double _fCellLengthY = (end_y - start_y) / (res - 1);
			double _fCellLengthZ = (end_z - start_z) / (res - 1);
			double derror = 0, derror2 = 0;
			std::vector<double> errorvec(res * res * res), errorvec2(res * res * res);
#pragma omp parallel for reduction(+:derror, derror2)
			for (int n = 0; n < res * res * res; n++)
			{
				int k = n / (res * res);
				int j = (n - k * res * res) / res;
				int i = n % res;
				double z = start_z + k * _fCellLengthZ;
				double y = start_y + j * _fCellLengthY;
				double x = start_x + i * _fCellLengthX;
				// std::cout << ptScalarField[n] << ' ' << mesh_distance.signed_distance({ x, y, z }).distance << std::endl;
				double gtdis = mesh_distance.signed_distance({ x, y, z }).distance;
				derror += fabs(ptScalarField[n] - gtdis) / (fabs(gtdis) + 1.0e-9);
				derror2 += fabs(ptScalarField[n] - gtdis);
				errorvec2[n] = fabs(ptScalarField[n] - gtdis);
				errorvec[n] = fabs(ptScalarField[n] - gtdis) / (fabs(gtdis) + 1.0e-9);
			}
			std::cout << "derror: " << derror / (res * res * res) << " derror2: " << derror2 / (res * res * res)<< std::endl;
			std::cout << *std::max_element(errorvec.begin(), errorvec.end()) << ' ' << *std::max_element(errorvec2.begin(), errorvec2.end()) << std::endl;
// 			int nbsamples = 100000;
// 			std::vector<TinyVector<float, 3>> samples(nbsamples);
// 			for (size_t i = 0; i < nbsamples; i++)
// 			{
// 				size_t vid = rand() % vertices.size();
// 				samples[i][0] = (float)(vertices[vid][0] + (rand() * 1.0 / RAND_MAX - 0.5) * 0.2);
// 				samples[i][1] = (float)(vertices[vid][1] + (rand() * 1.0 / RAND_MAX - 0.5) * 0.2);
// 				samples[i][2] = (float)(vertices[vid][2] + (rand() * 1.0 / RAND_MAX - 0.5) * 0.2);
// 			}
		
// 			std::vector<float> values;
// 			impfunc.scalar_value(samples, values);
// 			derror = 0, derror2 = 0;
// 			errorvec.resize(nbsamples), errorvec2.resize(nbsamples);
// #pragma omp parallel for reduction(+:derror, derror2)
// 			for (int n = 0; n < nbsamples; n++)
// 			{
// 				double gtdis = mesh_distance.signed_distance({ samples[n][0], samples[n][1], samples[n][2]}).distance;
// 				derror += fabs(values[n] - gtdis) / (fabs(gtdis) + 1.0e-9);
// 				derror2 += fabs(values[n] - gtdis);
// 				errorvec[n] = fabs(values[n] - gtdis);
// 				errorvec2[n] = fabs(values[n] - gtdis) / (fabs(gtdis) + 1.0e-9);
// 			}
// 			std::cout << "sampling derror: " << derror / nbsamples << " derror2: " << derror2 / nbsamples << std::endl;
// 			std::cout << *std::max_element(errorvec.begin(), errorvec.end()) << ' ' << *std::max_element(errorvec2.begin(), errorvec2.end()) << std::endl;

			return 0;
		}

		if (result.count("o"))
		{
			outputfilename = result["o"].as<std::string>();

			util.point_id_total_cache.reserve(10 * gridsize * gridsize);
			util.point_id_cache.reserve(util.buffer);

			MyImplicitFunc<float> impfunc(&module, &util);

			if (use_octree)
			{
				using namespace std::chrono;
				auto start = high_resolution_clock::now();

				int init_depth = depth, final_depth = depth;
				if (util.use_GPU)
					init_depth = std::min(7, depth);
				else
					init_depth = std::min(7, depth);

				std::vector<MyOctreeLayer<float>*> layer;
				layer.push_back(new MyOctreeLayer<float>(init_depth, &impfunc, start_x, end_x, start_y, end_y, start_z, end_z, isovalue, util.verbose));
				layer.back()->scan_octants(0, init_depth == final_depth ? false : true);
				//layer.back()->scan_octants(0, true);
				for (int step = init_depth + 1; step <= final_depth; step++)
				{
					layer.push_back(new MyOctreeLayer<float>(layer.back()->get_depth() + 1, &impfunc, start_x, end_x, start_y, end_y, start_z, end_z, isovalue, util.verbose));
					layer.back()->scan_octants(layer[layer.size() - 2], step == final_depth ? false : true);
					//layer.back()->scan_octants(layer[layer.size() - 2], true);
				}
				MyOctreeLayer<float>* mygrid = layer.back();

				std::cout << termcolor::yellow << "---------process grid edges!" << std::endl;
				std::cout << termcolor::cyan;

				std::string outputnamewithoutext = outputfilename;
				if (outputfilename.find_last_of(".") != std::string::npos)
					outputnamewithoutext = outputfilename.substr(0, outputfilename.find_last_of("."));

				char name[256];
				if (dump_all_layer_mesh)
				{
					for (int step = init_depth; step <= final_depth; step++)
					{
						mygrid = layer[step - init_depth];
						int d = pow(2, step);
						sprintf(name, "%s-%d.ply", outputnamewithoutext.c_str(), d);
						if (use_DC)
							mygrid->dcf_meshing(name);
						else
						{
							IsoSurfGen_FeatureAware(
								start_x, end_x, start_y, end_y, start_z, end_z,
								mygrid->get_grid_size(),
								&impfunc,
								mygrid->get_valid_cells(), mygrid->get_valid_cell_types(),
								mygrid->get_point_value_map(),
								name, isovalue, feature_angle, util.buffer, util.verbose
							);
						}
					}
				}
				else
				{
					if (use_DC)
						mygrid->dcf_meshing(outputfilename.c_str());
					else
					{
						IsoSurfGen_FeatureAware(
							start_x, end_x, start_y, end_y, start_z, end_z,
							mygrid->get_grid_size(),
							&impfunc,
							mygrid->get_valid_cells(), mygrid->get_valid_cell_types(),
							mygrid->get_point_value_map(),
							outputfilename.c_str(), isovalue, feature_angle, util.buffer, util.verbose
						);
					}
				}

				if (output_vtk)
				{
					for (int step = init_depth; step <= final_depth; step++)
					{
						mygrid = layer[step - init_depth];
						int d = pow(2, step);
						sprintf(name, "%s-%d.vtk", outputnamewithoutext.c_str(), d);
						mygrid->save_vtk(name);
					}
				}

				auto stop = high_resolution_clock::now();
				auto duration = duration_cast<microseconds>(stop - start);
				std::cout << termcolor::reset << "computation & file exporting time: " << termcolor::blue << termcolor::on_white << duration.count() / 1000000.0 << termcolor::reset << " secs" << std::endl;
				//char name[256];
				for (size_t i = 0; i < layer.size(); i++)
				{
					//sprintf_s(name, "c:\\Work\\Code\\DL\\libtorch\\lib\\grid-%d.obj", layer[i]->get_depth());
					//layer[i]->save_cubes(name);
					delete layer[i];
				}
			}
			else
			{
				//prepare grid data
				std::vector<float> ptScalarField(gridsize * gridsize * gridsize);
				MyOctreeLayer<float> mygrid(1, &impfunc, start_x, end_x, start_y, end_y, start_z, end_z, isovalue, util.verbose);
				std::cout << termcolor::red << "Step 1: " << termcolor::yellow << "process grid corners!" << std::endl;
				std::cout << termcolor::cyan;
				print_progress_bar(0.0f, util.verbose);

				using namespace std::chrono;
				auto start = high_resolution_clock::now();

				mygrid.compute_grid_point_values(gridsize, ptScalarField);

				//-----------------------------------------------------
				std::string extname = "";
				if (outputfilename.find_last_of(".") != std::string::npos)
					extname = outputfilename.substr(outputfilename.find_last_of(".") + 1);
				if (extname == "vtk")
				{
					std::ofstream mfile(outputfilename);
					mfile << "# vtk DataFile Version 3.0" << std::endl
						<< "Implicit function" << std::endl
						<< "ASCII" << std::endl
						<< "DATASET RECTILINEAR_GRID" << std::endl;
					mfile << "DIMENSIONS " << gridsize << ' ' << gridsize << ' ' << gridsize << std::endl;
					mfile << "X_COORDINATES " << gridsize << " float" << std::endl;
					for (int i = 0; i < gridsize - 1; i++)
					{
						mfile << start_x + i * (end_x - start_x) / (gridsize - 1) << ' ';
					}
					mfile << end_x << std::endl;
					mfile << "Y_COORDINATES " << gridsize << " float" << std::endl;
					for (int i = 0; i < gridsize - 1; i++)
					{
						mfile << start_y + i * (end_y - start_y) / (gridsize - 1) << ' ';
					}
					mfile << end_y << std::endl;
					mfile << "Z_COORDINATES " << gridsize << " float" << std::endl;
					for (int i = 0; i < gridsize - 1; i++)
					{
						mfile << start_z + i * (end_z - start_z) / (gridsize - 1) << ' ';
					}
					mfile << end_z << std::endl;

					mfile << "POINT_DATA " << ptScalarField.size() << std::endl;
					mfile << "SCALARS funcvalue float 1" << std::endl;
					mfile << "LOOKUP_TABLE default" << std::endl;
					for (size_t i = 0; i < ptScalarField.size() - 1; i++)
						mfile << ptScalarField[i] << ' ';
					mfile << ptScalarField.back() << std::endl;
					mfile.close();
				}
				else // write obj
				{
					if (util.use_GPU)
					{
						std::cout << termcolor::red << "Step 2: " << termcolor::yellow << "process grid edges!" << std::endl;
						if (util.verbose) std::cout << termcolor::cyan;
						IsoSurfGen(start_x, end_x, start_y, end_y, start_z, end_z, gridsize, gridsize, gridsize,
							ptScalarField,
							use_MC ? ImpGenType::MC33_C : ImpGenType::EXTENDED_MC_GPU,
							isovalue, feature_angle, &impfunc,
							outputfilename.c_str(), 0, 0, util.buffer, util.verbose);
					}
					else

						IsoSurfGen(start_x, end_x, start_y, end_y, start_z, end_z, gridsize, gridsize, gridsize,
							ptScalarField,
							use_MC ? ImpGenType::MC33_C : ImpGenType::EXTENDED_MC,
							isovalue, feature_angle, &impfunc,
							outputfilename.c_str(), 0, 0, util.buffer, util.verbose);
				}
				auto stop = high_resolution_clock::now();
				auto duration = duration_cast<microseconds>(stop - start);
				std::cout << termcolor::reset << "computation & file exporting time: " << termcolor::blue << termcolor::on_white << duration.count() / 1000000.0 << termcolor::reset << " secs" << std::endl;
			}

			Mesh_Vis mv(outputfilename, feature_angle);
		}
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << termcolor::reset << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
	return 0;
}
//////////////////////////////////////