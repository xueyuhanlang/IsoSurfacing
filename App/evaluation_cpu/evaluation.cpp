#include <omp.h>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>

#include "ThirdParty/cxxopts/cxxopts.hpp"
#include "ThirdParty/Happly/happly.h"
#include "ThirdParty/tmd/TriangleMeshDistance.h"

inline std::string GetFileExtension(const std::string &FileName)
{
    if (FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".") + 1);
    return "";
}

inline std::string GetFileNameWithoutExtension(const std::string &FileName)
{
    return FileName.substr(0, FileName.find_last_of("."));
}

bool load_mesh(const char meshfilename[],
              std::vector<std::array<double, 3>>& vertices,
              std::vector<std::array<int, 3>>& triangles)
{
    auto ext = GetFileExtension(meshfilename);
    if (ext == "ply")
    {
        happly::PLYData plyIn(meshfilename);
        vertices = plyIn.getVertexPositions();
        auto facets = plyIn.getFaceIndices<int>();

        triangles.resize(facets.size());
#pragma omp parallel for
        for (ptrdiff_t i = 0; i < (ptrdiff_t)facets.size(); i++)
        {
            triangles[i][0] = facets[i][0], triangles[i][1] = facets[i][1], triangles[i][2] = facets[i][2];
        }

        return true;
    }
    else if (ext == "obj")
    {
        vertices.resize(0), triangles.resize(0);
        FILE *m_pFile = fopen(meshfilename, "r");
        if (!m_pFile)
            return false;
        char temp[128];
        char *tok;

        fseek(m_pFile, 0, SEEK_SET);
        char pLine[1024];
        std::array<double, 3> v;
        while (fgets(pLine, 1024, m_pFile))
        {
            if (pLine[0] == 'v' && pLine[1] == ' ')
            {

                tok = strtok(pLine, " ");
                for (int i = 0; i < 3; i++)
                {
                    tok = strtok(NULL, " ");
                    strcpy(temp, tok);
                    temp[strcspn(temp, " ")] = 0;
                    v[i] = (double)atof(temp);
                }
                vertices.push_back(v);
            }
        }
        fseek(m_pFile, 0, SEEK_SET);
        std::vector<int> s_faceid;
        while (fgets(pLine, 1024, m_pFile))
        {
            if (pLine[0] == 'f')
            {
                s_faceid.resize(0);
                tok = strtok(pLine, " ");
                while ((tok = strtok(NULL, " ")) != NULL)
                {
                    strcpy(temp, tok);
                    temp[strcspn(temp, "/")] = 0;
                    int id = (int)strtol(temp, NULL, 10) - 1;
                    s_faceid.push_back(id);
                }
                std::array<int, 3> f;
                f[0] = s_faceid[0], f[1] = s_faceid[1], f[2] = s_faceid[2];
                triangles.push_back(f);
            }
        }
        fclose(m_pFile);
        return true;
    }
    else
        return false;
}

int main(int argc, char **argv)
{
    try
    {
        unsigned int npoints = 131072;
        cxxopts::Options options("Evaluation", "Evaluation (author: Yang Liu, Email: yangliu@microsoft.com)");
        options
            .positional_help("[optional args]")
            .show_positional_help()
            .allow_unrecognised_options()
            .add_options()
            ("g,gtmesh", "GT mesh (*.obj,*.ply)", cxxopts::value<std::string>())
            ("r,rmesh", "reconstruction mesh (*.obj,*.ply)", cxxopts::value<std::string>())
            ("n,npoints", "number of points for IoU computation (default: 131072)", cxxopts::value<unsigned int>());
        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cerr << options.help({"", "Group"}) << std::endl;
            exit(0);
        }

		if (result.count("n"))
		{
			npoints = result["n"].as<unsigned int>();
		}

        if (result.count("r") && result.count("g"))
        {
            auto gt_filename = result["g"].as<std::string>();
            auto recon_filename = result["r"].as<std::string>();

            std::vector<std::array<double, 3>> gt_vertices, recon_vertices;
            std::vector<std::array<int, 3>> gt_triangles, recon_triangles;
            if (!load_mesh(gt_filename.c_str(), gt_vertices, gt_triangles))
            {
                std::cout << "the mesh file format is not supported!" << std::endl;
                exit(-1);
            }
            if (!load_mesh(recon_filename.c_str(), recon_vertices, recon_triangles))
            {
                std::cout << "the mesh file format is not supported!" << std::endl;
                exit(-1);
            }

            tmd::TriangleMeshDistance gt_mesh_distance(gt_vertices, gt_triangles);
            tmd::TriangleMeshDistance recon_mesh_distance(recon_vertices, recon_triangles);

            std::vector<std::array<double, 3>> pvec(npoints);
            std::mt19937 eng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
            std::uniform_real_distribution<> unif(-1, 1);
            for (size_t i = 0; i < npoints; i++)
            {
                pvec[i][0] = unif(eng), pvec[i][1] = unif(eng), pvec[i][2] = unif(eng);
            }
            double derror = 0;
            ptrdiff_t n_intersect = 0, n_union = 0;
#pragma omp parallel for reduction(+ \
                                   : derror, n_intersect, n_union)
            for (int n = 0; n < npoints; n++)
            {
                double gtdis = gt_mesh_distance.signed_distance({pvec[n][0], pvec[n][1], pvec[n][2]}).distance;
                double rtdis = recon_mesh_distance.signed_distance({pvec[n][0], pvec[n][1], pvec[n][2]}).distance;

                derror += fabs(rtdis - gtdis) / (fabs(gtdis) + 1.0e-9);
                if (rtdis <= 0 && gtdis <= 0)
                {
                    n_intersect++;
                }
                if (rtdis <= 0 || gtdis <= 0)
                    n_union++;
            }

            std::ofstream outputfile(GetFileNameWithoutExtension(recon_filename) + "_eval.txt");
            outputfile << derror / npoints << std::endl;
            outputfile << (double)n_intersect / (n_union == 0 ? 1 : n_union) << std::endl;
            outputfile.close();
            std::cout << derror / npoints << ' ' << (double)n_intersect / (n_union == 0 ? 1 : n_union) << std::endl;
        }
        else
        {
            std::cerr << "The inputs are not valid!" << std::endl;
            exit(0);
        }
    }
    catch (const cxxopts::OptionException &e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}