#pragma once
#include "ThirdParty/Happly/happly.h"
#include <cmath>
#include <omp.h>
#include "ThirdParty/IsoSurfGen/TinyVector.h"
#include <map>

class Mesh_Vis
{
public:
    Mesh_Vis(const std::string &plyfilename, const float sharp_angle_in_degree = 10)
    {
        std::string filenamewithoutext = plyfilename.substr(0, plyfilename.find_last_of("."));

        happly::PLYData plyIn(plyfilename);
        std::vector<TinyVector<double, 3>> vertices = plyIn.getVertexPositions<TinyVector<double, 3>>();
        std::vector<std::vector<int>> facets = plyIn.getFaceIndices<int>();
        const double cos_threshold = cos(sharp_angle_in_degree * M_PI / 180);
        std::vector<std::pair<int, int>> sharp_edges;
        std::vector<TinyVector<double, 3>> face_normals(facets.size());
#pragma omp parallel for
        for (ptrdiff_t i = 0; i < (ptrdiff_t)facets.size(); i++)
        {
            for (size_t k = 0; k < facets[i].size(); k++)
            {
                face_normals[i] += vertices[facets[i][k]].Cross(vertices[facets[i][(k + 1) % facets[i].size()]]);
            }
            face_normals[i].Normalize();
        }
        std::map<std::pair<int, int>, std::vector<size_t>> edge2facets;
        for (ptrdiff_t i = 0; i < (ptrdiff_t)facets.size(); i++)
        {
            for (size_t k = 0; k < facets[i].size(); k++)
            {
                auto v0 = facets[i][k], v1 = facets[i][(k + 1) % facets[i].size()];
                if (v0 > v1)
                    std::swap(v0, v1);
                edge2facets[std::make_pair(v0, v1)].push_back(i);
            }
        }

        std::vector<bool> vertex_tag(vertices.size(), false);
        for (const auto &e : edge2facets)
        {
            if (e.second.size() == 2 && fabs(face_normals[e.second[0]].Dot(face_normals[e.second[1]])) <= cos_threshold)
            {
                sharp_edges.push_back(e.first);
                vertex_tag[e.first.first] = true, vertex_tag[e.first.second] = true;
            }
        }

        std::string sharpedgefilename = filenamewithoutext + "_sharpedge" + ".obj";
        std::ofstream edgestream(sharpedgefilename);
        std::vector<ptrdiff_t> vertex_indice(vertices.size(), -1);
        ptrdiff_t vcounter = 0;
        for (size_t i = 0; i < vertices.size(); i++)
        {
            if (vertex_tag[i])
            {
                edgestream << "v " << vertices[i] << std::endl;
                vertex_indice[i] = vcounter;
                vcounter++;
            }
        }
        for (const auto &e : sharp_edges)
        {
            edgestream << "l " << vertex_indice[e.first] + 1 << ' ' << vertex_indice[e.second] + 1 << std::endl;
        }
        edgestream.close();

        std::string meshlabprojname = filenamewithoutext + ".mlp";
        std::ofstream mlpfile(meshlabprojname);
        mlpfile << "<!DOCTYPE MeshLabDocument>\r"
                << "<MeshLabProject>\r"
                << "<MeshGroup>\r";

        mlpfile << "<MLMesh label=\"" << plyfilename << "\" visible=\"" << 1 << "\" filename=\"" << plyfilename<< "\">\r";
        mlpfile << "<MLMatrix44>\r"
                << "1 0 0 0 \r"
                << "0 1 0 0 \r"
                << "0 0 1 0 \r"
                << "0 0 0 1 \r"
                << "</MLMatrix44>\r";
        mlpfile << "<RenderingOption solidColor=\"0 155 233 255\">"
                << "100001000000000000000000000001010100000010100000000001111011110000001001"
                << "</RenderingOption>\r";
        mlpfile << "</MLMesh>\r";

        mlpfile << "<MLMesh label=\"" << sharpedgefilename << "\" visible=\"" << 1 << "\" filename=\"" << sharpedgefilename<< "\">\r";
        mlpfile << "<MLMatrix44>\r"
                << "1 0 0 0 \r"
                << "0 1 0 0 \r"
                << "0 0 1 0 \r"
                << "0 0 0 1 \r"
                << "</MLMatrix44>\r";
        mlpfile << "<RenderingOption wireWidth=\"5\" wireColor=\"0 0 0 255\">"
                << "001001000000010000000100000001000100000010100000110100111011110000001001"
                << "</RenderingOption>\r";
        mlpfile << "</MLMesh>\r";


        mlpfile << "</MeshGroup>\r"
                << "<RasterGroup/>\r"
                << "</MeshLabProject>\r";
        mlpfile.close();
    }
};