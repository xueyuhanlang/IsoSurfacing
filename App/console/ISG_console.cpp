
#include "ThirdParty/cxxopts/cxxopts.hpp"
#include <algorithm>
#include "MyImplicitFunc.h"
#include "ThirdParty/IsoSurfGen/IsoSurfGen.h"

#define TEST_IMP

#ifdef TEST_IMP

int main()
{
	MyImplicitFunc<float> impfunc;
	IsoSurfGen(-2.0f, 2.0f, -2.0f, 2.0f, -2.0f, 2.0f, 64, 64, 64,
		ImpGenType::EXTENDED_MC,
		0.0f, 30.0f, &impfunc,
		"test_imp1.obj");

	//MyImplicitFunc<float> impfunc;
	//IsoSurfGen(-2.0f, 2.0f, -2.0f, 2.0f, -2.0f, 2.0f, 64, 64, 64,
	//	ImpGenType::EXTENDED_MC_GPU,
	//	0.0f, 30.0f, &impfunc,
	//	"test_imp.obj");
	return 0;
}

#else

int main(int argc, char** argv)
{
	float feature_angle = 30;
	ImpGenType method_type = ImpGenType::MC33_C;

	try
	{
		cxxopts::Options options("IsoSurfaceGen", "IsoSurfaceGenator (author: Yang Liu, Email: yangliu@microsoft.com)");
		options
			.positional_help("[optional args]")
			.show_positional_help()
			.allow_unrecognised_options()
			.add_options()
			("i,input", "input volume data (*.dat)", cxxopts::value<std::string>())
			("m,method", "extraction method (MC: original MC, MC1: improved MC; MC2: fastest improved MC (default); EMC: feature-aware MC; DMC: dual MC", cxxopts::value<std::string>())
			("v,isovalue", "isovalue, (default : 0.0), valid for MC, MC2, DMC only", cxxopts::value<float>())
			("a,angle", "feature angle threshold for EMC  (default: 30 degree)", cxxopts::value<float>())
			("o,output", "output mesh (obj format) or vtk format (for volumetric visualization)", cxxopts::value<std::string>())
			("h,help", "Print help");

		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}
		float iso_value = 0.0;
		std::string inputfilename, outputfilename;

		if (result.count("i"))
		{
			inputfilename = result["i"].as<std::string>();
		}

		if (result.count("m"))
		{
			auto& method = result["m"].as<std::string>();
			if (method == "MC")
				method_type = ImpGenType::MC;
			else if (method == "MC1")
				method_type = ImpGenType::MC33_CPLUS;
			else if (method == "MC2")
				method_type = ImpGenType::MC33_C;
			else if (method == "EMC")
				method_type = ImpGenType::EXTENDED_MC;
			else if (method == "DMC")
				method_type = ImpGenType::DUAL_MC;
		}

		if (result.count("a"))
			feature_angle = std::max(0.0f, result["a"].as<float>());

		if (result.count("v"))
			iso_value = std::max(0.0f, result["v"].as<float>());

		if (result.count("o"))
		{
			outputfilename = result["o"].as<std::string>();
			IsoSurfGen(inputfilename.c_str(), method_type, iso_value, feature_angle, outputfilename.c_str(), 0, 0);
		}
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}

	return 0;
}
#endif