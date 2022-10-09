#include "IsoSurfGen.h"
#include <fstream>
#include <iostream>
#include <omp.h>
#include "CIsoSurface/CIsoSurface.h"
#include "marchingcube33/MarchingCubes.h"
#include "marchingcube33/marching_cubes_33.h"
#include "dualmc/dualmc.h"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <IsoEx/Implicits/ImplicitInterface.hh>
#include <IsoEx/Grids/ScalarGridT.hh>
#include <IsoEx/Grids/ImplicitGridT.hh>
#include <IsoEx/Extractors/MarchingCubesT.hh>
#include <IsoEx/Extractors/ExtendedMarchingCubesT.hh>
#include <IsoEx/Extractors/AExtendedMarchingCubesT.hh>
#include <IsoEx/Extractors/AOExtendedMarchingCubesT.hh>

#define VectorType OpenMesh::Vec3f

//-----------------------------------------------------------------------------

using namespace IsoEx;
using namespace OpenMesh;
using namespace OpenMesh::IO;

//-----------------------------------------------------------------------------

// Define the mesh to be used: need vertex normal and status for EMC
struct MyTraits : public DefaultTraits
{
	VertexAttributes(Attributes::Normal | Attributes::Status);
	HalfedgeAttributes(Attributes::PrevHalfedge);

	/// Use float precision points
	typedef VectorType Point;
	/// Use float precision Normals
	typedef VectorType Normal;
};
typedef TriMesh_ArrayKernelT<MyTraits> MyMesh;

bool IsoSurfGen(
	float start_x, float end_x,
	float start_y, float end_y,
	float start_z, float end_z,
	unsigned int nGridX, unsigned int nGridY, unsigned int nGridZ,
	std::vector<float> &ptScalarField,
	ImpGenType method_type,
	float iso_value, float feature_angle,
	ImplicitFunc<float> *impfunc,
	const char *outobjfile,
	std::vector<TinyVector<float, 3>> *vertices,
	std::vector<unsigned int> *trifacets,
	const int buffer,
	bool verbose)
{
	if (start_x > end_x || start_y > end_y || start_z > end_z ||
		nGridX < 2 || nGridY < 2 || nGridZ < 2)
		return false;
	if (outobjfile == 0 && vertices == 0 && trifacets == 0)
		return false;

	if (ptScalarField.size() != nGridX * nGridY * nGridZ)
		return false;

	unsigned int nCellsX = nGridZ - 1, nCellsY = nGridY - 1, nCellsZ = nGridZ - 1;
	float fCellLengthX = (end_x - start_x) / nCellsX;
	float fCellLengthY = (end_y - start_y) / nCellsY;
	float fCellLengthZ = (end_z - start_z) / nCellsZ;

	std::string filename(outobjfile);
	std::string extname = "";
	if (filename.find_last_of(".") != std::string::npos)
		extname = filename.substr(filename.find_last_of(".") + 1);

	std::ofstream meshfile;
	if (outobjfile)
		meshfile.open(outobjfile);

	if (method_type == ImpGenType::MC)
	{
		CIsoSurface<float> ciso;

		ciso.GenerateSurface(&ptScalarField[0], (float)iso_value, nCellsX, nCellsY, nCellsZ, fCellLengthX, fCellLengthY, fCellLengthZ);

		// construct the mesh
		if (ciso.IsSurfaceValid() && ciso.m_nTriangles > 0)
		{
			if (meshfile.is_open())
			{
				if (extname == "obj")
				{
					for (unsigned int i = 0; i < ciso.m_nVertices; i++)
					{
						meshfile << "v "
								 << ciso.m_ppt3dVertices[i][0] + start_x << ' '
								 << ciso.m_ppt3dVertices[i][1] + start_y << ' '
								 << ciso.m_ppt3dVertices[i][2] + start_z << std::endl;
					}

					for (unsigned int i = 0; i < ciso.m_nTriangles; i++)
					{
						meshfile << "f "
								 << ciso.m_piTriangleIndices[3 * i] + 1 << ' '
								 << ciso.m_piTriangleIndices[3 * i + 1] + 1 << ' '
								 << ciso.m_piTriangleIndices[3 * i + 2] + 1 << std::endl;
					}
				}
				else
				{
					meshfile << "ply\nformat ascii 1.0\nelement vertex " << ciso.m_nVertices << "\nproperty float x\nproperty float y\nproperty float z\n"
							 << "element face " << ciso.m_nTriangles << "\nproperty list uchar int vertex_index\nend_header\n";
					for (unsigned int i = 0; i < ciso.m_nVertices; i++)
					{
						meshfile << ciso.m_ppt3dVertices[i][0] + start_x << ' '
								 << ciso.m_ppt3dVertices[i][1] + start_y << ' '
								 << ciso.m_ppt3dVertices[i][2] + start_z << std::endl;
					}

					for (unsigned int i = 0; i < ciso.m_nTriangles; i++)
					{
						meshfile << "3 "
								 << ciso.m_piTriangleIndices[3 * i] << ' '
								 << ciso.m_piTriangleIndices[3 * i + 1] << ' '
								 << ciso.m_piTriangleIndices[3 * i + 2] << std::endl;
					}
				}
			}
			if (vertices && trifacets)
			{
				vertices->resize(ciso.m_nVertices);
				// trifacets->resize(3*ciso.m_nTriangles);
				for (unsigned int i = 0; i < ciso.m_nVertices; i++)
				{
					(*vertices)[i][0] = ciso.m_ppt3dVertices[i][0] + start_x;
					(*vertices)[i][1] = ciso.m_ppt3dVertices[i][1] + start_y;
					(*vertices)[i][2] = ciso.m_ppt3dVertices[i][2] + start_z;
				}
				trifacets->assign(ciso.m_piTriangleIndices, ciso.m_piTriangleIndices + 3 * ciso.m_nTriangles);
			}
		}
	}
	else if (method_type == ImpGenType::MC33_CPLUS)
	{
		MarchingCubes mc;
		mc.set_resolution(nCellsX + 1, nCellsY + 1, nCellsZ + 1);
		mc.init_all();

#pragma omp parallel for
		for (ptrdiff_t i = 0; i < (ptrdiff_t)ptScalarField.size(); i++)
			ptScalarField[i] -= iso_value;

		mc.set_ext_data(&ptScalarField[0]);

		mc.run();
		mc.clean_temps();

		if (meshfile.is_open())
		{
			if (extname == "obj")
			{
				for (int i = 0; i < mc.nverts(); i++)
				{
					meshfile << "v "
							 << mc.vertices()[i].x * fCellLengthX + start_x << ' '
							 << mc.vertices()[i].y * fCellLengthY + start_y << ' '
							 << mc.vertices()[i].z * fCellLengthZ + start_z << std::endl;
				}

				for (int i = 0; i < mc.ntrigs(); i++)
				{
					meshfile << "f "
							 << mc.triangles()[i].v1 + 1 << ' '
							 << mc.triangles()[i].v2 + 1 << ' '
							 << mc.triangles()[i].v3 + 1 << std::endl;
				}
			}
			else
			{
				meshfile << "ply\nformat ascii 1.0\nelement vertex " << mc.nverts() << "\nproperty float x\nproperty float y\nproperty float z\n"
						 << "element face " << mc.ntrigs() << "\nproperty list uchar int vertex_index\nend_header\n";
				for (int i = 0; i < mc.nverts(); i++)
				{
					meshfile
						<< mc.vertices()[i].x * fCellLengthX + start_x << ' '
						<< mc.vertices()[i].y * fCellLengthY + start_y << ' '
						<< mc.vertices()[i].z * fCellLengthZ + start_z << std::endl;
				}

				for (int i = 0; i < mc.ntrigs(); i++)
				{
					meshfile << "3 "
							 << mc.triangles()[i].v1 << ' '
							 << mc.triangles()[i].v2 << ' '
							 << mc.triangles()[i].v3 << std::endl;
				}
			}
		}

		if (vertices && trifacets)
		{
			vertices->resize(mc.nverts());
			trifacets->resize(3 * mc.ntrigs());

			for (int i = 0; i < mc.nverts(); i++)
			{
				(*vertices)[i][0] = mc.vertices()[i].x * fCellLengthX + start_x;
				(*vertices)[i][1] = mc.vertices()[i].y * fCellLengthY + start_y;
				(*vertices)[i][2] = mc.vertices()[i].z * fCellLengthZ + start_z;
			}

			for (int i = 0; i < mc.ntrigs(); i++)
			{
				(*trifacets)[3 * i] = mc.triangles()[i].v1;
				(*trifacets)[3 * i + 1] = mc.triangles()[i].v2;
				(*trifacets)[3 * i + 2] = mc.triangles()[i].v3;
			}
		}
		mc.clean_all();
	}
	else if (method_type == ImpGenType::MC33_C)
	{
		_GRD *Z;

		unsigned short int nx, ny, nz;
		Z = (_GRD *)malloc(sizeof(_GRD));
		Z->periodic = 0;
		memset(Z->r0, 0, sizeof(Z->r0));
		for (int i = 0; i < 3; i++)
			Z->d[i] = 1;
		Z->d[0] = fCellLengthX, Z->d[1] = fCellLengthY, Z->d[2] = fCellLengthZ;
		nx = nCellsX + 1, ny = nCellsY + 1, nz = nCellsZ + 1;
		Z->L[0] = Z->N[0] = nCellsX;
		Z->L[1] = Z->N[1] = nCellsY;
		Z->L[2] = Z->N[2] = nCellsZ;
		Z->F = (GRD_data_type ***)malloc(nz * sizeof(void *));
		for (int k = 0; k < nz; k++)
		{
			Z->F[k] = (GRD_data_type **)malloc(ny * sizeof(void *));
			for (int j = 0; j < ny; j++)
				Z->F[k][j] = (GRD_data_type *)malloc(nx * sizeof(GRD_data_type));
		}

		int size = (nCellsX + 1) * (nCellsY + 1) * (nCellsZ + 1);

#pragma omp parallel for
		for (int n = 0; n < size; n++)
		{
			int k = n / ((nCellsY + 1) * (nCellsX + 1));
			int j = (n - k * (nCellsY + 1) * (nCellsX + 1)) / (nCellsX + 1);
			int i = n % (nCellsX + 1);
			Z->F[k][j][i] = ptScalarField[n];
		}

		surface *S;
		S = calculate_isosurface(Z, (float)iso_value);

		//----
		if (Z->F)
		{
			for (int k = 0; k <= Z->N[2]; k++)
			{
				if (Z->F[k])
					for (int j = 0; j <= Z->N[1]; j++)
						free(Z->F[k][j]);
				free(Z->F[k]);
			}
			free(Z->F);
		}
		free(Z);

		// construct the mesh
		if (S && S->nT > 0)
		{
			if (meshfile.is_open())
			{
				if (extname == "obj")
				{
					for (int i = 0; i <= S->nV; i++)
					{
						meshfile << "v "
								 << S->V[i >> _MC_N][i & _MC_A][0] + start_x << ' '
								 << S->V[i >> _MC_N][i & _MC_A][1] + start_y << ' '
								 << S->V[i >> _MC_N][i & _MC_A][2] + start_z << std::endl;
					}

					for (int i = 0; i <= S->nT; i++)
					{
						meshfile << "f "
								 << S->T[i >> _MC_N][i & _MC_A][1] + 1 << ' '
								 << S->T[i >> _MC_N][i & _MC_A][0] + 1 << ' '
								 << S->T[i >> _MC_N][i & _MC_A][2] + 1 << std::endl;
					}
				}
				else
				{
					meshfile << "ply\nformat ascii 1.0\nelement vertex " << S->nV + 1 << "\nproperty float x\nproperty float y\nproperty float z\n"
							 << "element face " << S->nT + 1 << "\nproperty list uchar int vertex_index\nend_header\n";
					for (int i = 0; i <= S->nV; i++)
					{
						meshfile
							<< S->V[i >> _MC_N][i & _MC_A][0] + start_x << ' '
							<< S->V[i >> _MC_N][i & _MC_A][1] + start_y << ' '
							<< S->V[i >> _MC_N][i & _MC_A][2] + start_z << std::endl;
					}

					for (int i = 0; i <= S->nT; i++)
					{
						meshfile << "3 "
								 << S->T[i >> _MC_N][i & _MC_A][1] << ' '
								 << S->T[i >> _MC_N][i & _MC_A][0] << ' '
								 << S->T[i >> _MC_N][i & _MC_A][2] << std::endl;
					}
				}
			}

			if (vertices && trifacets)
			{
				vertices->resize(S->nV + 1);
				trifacets->resize((S->nT + 1) * 3);

				for (int i = 0; i <= S->nV; i++)
				{
					(*vertices)[i][0] = S->V[i >> _MC_N][i & _MC_A][0] + start_x;
					(*vertices)[i][1] = S->V[i >> _MC_N][i & _MC_A][1] + start_y;
					(*vertices)[i][2] = S->V[i >> _MC_N][i & _MC_A][2] + start_z;
				}

				for (int i = 0; i <= S->nT; i++)
				{
					(*trifacets)[3 * i] = S->T[i >> _MC_N][i & _MC_A][0];
					(*trifacets)[3 * i + 1] = S->T[i >> _MC_N][i & _MC_A][1];
					(*trifacets)[3 * i + 2] = S->T[i >> _MC_N][i & _MC_A][2];
				}
			}
			free_surface_memory(&S);
		}
	}
	else if (method_type == ImpGenType::EXTENDED_MC || method_type == ImpGenType::EXTENDED_MC_GPU)
	{
		MyMesh mesh;

		if (impfunc == 0)
		{
			ScalarGridT<float, VectorType> grid(
				VectorType(start_x, start_y, start_z), // origin
				VectorType(end_x - start_x, 0, 0),	   // x-axis
				VectorType(0, end_y - start_y, 0),	   // y-axis
				VectorType(0, 0, end_z - start_z),	   // z-axis
				nGridX, nGridY, nGridZ);			   // resolution

#pragma omp parallel for
			for (ptrdiff_t i = 0; i < (ptrdiff_t)ptScalarField.size(); i++)
				ptScalarField[i] -= iso_value;

			grid.get_values().assign(ptScalarField.begin(), ptScalarField.end());

			extended_marching_cubes(grid, mesh, iso_value, feature_angle);
		}
		else
		{
			ImplicitInterface<VectorType> mfunc(impfunc);
			ImplicitGrid<VectorType> grid(mfunc,								 // implicit
										  VectorType(start_x, start_y, start_z), // origin
										  VectorType(end_x - start_x, 0, 0),	 // x-axis
										  VectorType(0, end_y - start_y, 0),	 // y-axis
										  VectorType(0, 0, end_z - start_z),	 // z-axis
										  nGridX, nGridY, nGridZ);				 // resolution

			grid.get_is_inside_cache().resize(ptScalarField.size());
			grid.get_scalar_distance_cache().resize(ptScalarField.size());
#pragma omp parallel for
			for (int i = 0; i < (int)ptScalarField.size(); i++)
			{
				grid.get_scalar_distance_cache()[i] = ptScalarField[i];
				grid.get_is_inside_cache()[i] = (ptScalarField[i] <= iso_value);
			}

			if (method_type == ImpGenType::EXTENDED_MC)
			{
				extended_marching_cubes(grid, mesh, iso_value, feature_angle);
			}
			else
			{
				advanced_extended_marching_cubes<MyMesh, ImplicitGrid<VectorType>, float>(grid, mesh, iso_value, feature_angle, buffer, verbose);
			}
		}

		unsigned int scale = std::max(std::max(nCellsX, nCellsY), nCellsZ);

		if (meshfile.is_open())
		{
			meshfile.close();
			std::string filename(outobjfile);
			std::string extname = "";
			if (filename.find_last_of(".") != std::string::npos)
				extname = filename.substr(filename.find_last_of(".") + 1);
			if (extname == "obj")
			{
				write_mesh(mesh, outobjfile);
			}
			else if (extname == "ply")
			{
				OpenMesh::IO::Options op;
				op += IO::Options::Binary;
				op += IO::Options::Swap;
				write_mesh(mesh, outobjfile, op);
			}

			/*
			for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
			{
				meshfile << "v "
					<< mesh.point(*v_it)[0] << ' '
					<< mesh.point(*v_it)[1] << ' '
					<< mesh.point(*v_it)[2] << std::endl;
			}
			for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
			{
				meshfile << "vn "
					<< mesh.normal(*v_it)[0] << ' '
					<< mesh.normal(*v_it)[1] << ' '
					<< mesh.normal(*v_it)[2] << std::endl;
			}
			//for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
			//{
			//	MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(f_it);
			//	meshfile << "f";
			//	for (; fh_it.is_valid(); ++fh_it) {
			//		meshfile << " " << mesh.from_vertex_handle(fh_it).idx() + 1 << "//" << mesh.from_vertex_handle(fh_it).idx() + 1;
			//	}
			//	meshfile << std::endl;
			//}
			for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
			{
				MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(f_it);
				meshfile << "f";
				for (; fh_it.is_valid(); ++fh_it) {
					meshfile << " " << mesh.from_vertex_handle(fh_it).idx() + 1;
				}
				meshfile << std::endl;
			}
			*/
		}
		if (vertices && trifacets)
		{
			vertices->resize(mesh.n_vertices());
			trifacets->resize(3 * mesh.n_faces());

			int i = 0;
			for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it, i++)
			{
				(*vertices)[i][0] = mesh.point(*v_it)[0];
				(*vertices)[i][1] = mesh.point(*v_it)[1];
				(*vertices)[i][2] = mesh.point(*v_it)[2];
			}
			i = 0;
			for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
			{
				MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(f_it);
				for (; fh_it.is_valid(); ++fh_it)
				{
					(*trifacets)[i] = mesh.from_vertex_handle(fh_it).idx();
					i++;
				}
			}
		}
	}
	else if (method_type == ImpGenType::DUAL_MC)
	{
		std::vector<dualmc::Vertex> dvertices;
		std::vector<dualmc::Quad> quads;
		dualmc::DualMC<float> builder;
		builder.build(&ptScalarField[0], nGridX, nGridY, nGridZ,
					  iso_value, true, false, dvertices, quads);
		unsigned int scale = std::max(std::max(nCellsX, nCellsY), nCellsZ);
		// write vertices
		if (meshfile.is_open())
		{
			if (extname == "obj")
			{
				for (auto const &v : dvertices)
				{
					meshfile << "v " << v.x * fCellLengthX + start_x << ' ' << v.y * fCellLengthY + start_y << ' ' << v.z * fCellLengthZ + start_z << std::endl;
				}

				// write quad indices
				for (auto const &q : quads)
				{
					meshfile << "f " << (q.i0 + 1) << ' ' << (q.i1 + 1) << ' ' << (q.i2 + 1) << ' ' << (q.i3 + 1) << '\n';
				}
			}
			else
			{
				meshfile << "ply\nformat ascii 1.0\nelement vertex " << dvertices.size() << "\nproperty float x\nproperty float y\nproperty float z\n"
						 << "element face " << quads.size() << "\nproperty list uchar int vertex_index\nend_header\n";
				for (auto const &v : dvertices)
				{
					meshfile << v.x * fCellLengthX + start_x << ' ' << v.y * fCellLengthY + start_y << ' ' << v.z * fCellLengthZ + start_z << std::endl;
				}

				// write quad indices
				for (auto const &q : quads)
				{
					meshfile << "4 " << q.i0 << ' ' << q.i1 << ' ' << q.i2 << ' ' << q.i3 << '\n';
				}
			}
		}
		if (vertices && trifacets)
		{
			vertices->resize(dvertices.size());
			trifacets->resize(quads.size() * 4);

			int i = 0;
			for (auto const &v : dvertices)
			{
				(*vertices)[i][0] = v.x * fCellLengthX + start_x;
				(*vertices)[i][1] = v.y * fCellLengthY + start_y;
				(*vertices)[i][2] = v.z * fCellLengthZ + start_z;
				i++;
			}

			i = 0;
			// write quad indices
			for (auto const &q : quads)
			{
				(*trifacets)[i] = q.i0;
				(*trifacets)[i + 1] = q.i1;
				(*trifacets)[i + 2] = q.i2;
				(*trifacets)[i + 3] = q.i3;
				i += 4;
			}
		}
	}

	if (meshfile.is_open())
		meshfile.close();

	return true;
}

bool IsoSurfGen(
	float start_x, float end_x,
	float start_y, float end_y,
	float start_z, float end_z,
	unsigned int nGridX, unsigned int nGridY, unsigned int nGridZ,
	ImpGenType method_type,
	float iso_value, float feature_angle,
	ImplicitFunc<float> *impfunc,
	const char *outobjfile,
	std::vector<TinyVector<float, 3>> *vertices,
	std::vector<unsigned int> *trifacets,
	const int buffer,
	bool verbose)
{
	if (impfunc == 0)
		return false;
	if (start_x > end_x || start_y > end_y || start_z > end_z ||
		nGridX < 2 || nGridY < 2 || nGridZ < 2)
		return false;
	if (outobjfile == 0 && vertices == 0 && trifacets == 0)
		return false;
	// assert(method_type != ImpGenType::EXTENDED_MC_GPU);

	unsigned int nCellsX = nGridZ - 1, nCellsY = nGridY - 1, nCellsZ = nGridZ - 1;
	float fCellLengthX = (end_x - start_x) / nCellsX;
	float fCellLengthY = (end_y - start_y) / nCellsY;
	float fCellLengthZ = (end_z - start_z) / nCellsZ;
	std::vector<float> ptScalarField;

	ptScalarField.resize((nCellsX + 1) * (nCellsY + 1) * (nCellsZ + 1));
#pragma omp parallel for
	for (int n = 0; n < (int)ptScalarField.size(); n++)
	{
		int k = n / ((nCellsY + 1) * (nCellsX + 1));
		int j = (n - k * (nCellsY + 1) * (nCellsX + 1)) / (nCellsX + 1);
		int i = n % (nCellsX + 1);
		float z = start_z + k * fCellLengthZ;
		float y = start_y + j * fCellLengthY;
		float x = start_x + i * fCellLengthX;
		ptScalarField[n] = impfunc->scalar_value(TinyVector<float, 3>(x, y, z));
	}

	IsoSurfGen(start_x, end_x, start_y, end_y, start_z, end_z,
			   nGridX, nGridY, nGridZ,
			   ptScalarField, method_type,
			   iso_value, feature_angle, impfunc,
			   outobjfile, vertices, trifacets, buffer, verbose);

	return true;
}

bool IsoSurfGen(
	const char *inputfile,
	ImpGenType method_type,
	float iso_value, float feature_angle,
	const char *outputfile,
	std::vector<TinyVector<float, 3>> *vertices,
	std::vector<unsigned int> *trifacets,
	const int buffer,
	bool verbose)
{
	if (inputfile == 0)
		return false;
	if (outputfile == 0 && vertices == 0 && trifacets == 0)
		return false;

	float start_x, end_x, start_y, end_y, start_z, end_z;
	unsigned int nGridX, nGridY, nGridZ;
	std::vector<float> ptScalarField;

	try
	{
		std::ifstream mfile(inputfile);
		mfile >> start_x >> end_x >> start_y >> end_y >> start_z >> end_z >> nGridX >> nGridY >> nGridZ;
		ptScalarField.resize(nGridX * nGridY * nGridZ);
		for (size_t i = 0; i < ptScalarField.size(); i++)
			mfile >> ptScalarField[i];
		mfile.close();

		std::string outfilename(outputfile);
		std::string extname = "";
		if (outfilename.find_last_of(".") != std::string::npos)
			extname = outfilename.substr(outfilename.find_last_of(".") + 1);
		if (extname == "vtk")
		{
			std::ofstream mfile(outputfile);
			mfile << "# vtk DataFile Version 3.0" << std::endl
				  << "Implicit function" << std::endl
				  << "ASCII" << std::endl
				  << "DATASET RECTILINEAR_GRID" << std::endl;
			mfile << "DIMENSIONS " << nGridX << ' ' << nGridY << ' ' << nGridZ << std::endl;
			mfile << "X_COORDINATES " << nGridX << " float" << std::endl;
			for (unsigned int i = 0; i < nGridX - 1; i++)
			{
				mfile << start_x + i * (end_x - start_x) / (nGridX - 1) << ' ';
			}
			mfile << end_x << std::endl;
			mfile << "Y_COORDINATES " << nGridY << " float" << std::endl;
			for (unsigned int i = 0; i < nGridY - 1; i++)
			{
				mfile << start_y + i * (end_y - start_y) / (nGridY - 1) << ' ';
			}
			mfile << end_y << std::endl;
			mfile << "Z_COORDINATES " << nGridZ << " float" << std::endl;
			for (unsigned int i = 0; i < nGridZ - 1; i++)
			{
				mfile << start_z + i * (end_z - start_z) / (nGridZ - 1) << ' ';
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
		else if (extname == "obj")
		{
			IsoSurfGen(
				start_x, end_x, start_y, end_y, start_z, end_z,
				nGridX, nGridY, nGridZ,
				ptScalarField, method_type,
				iso_value, feature_angle, 0,
				outputfile, vertices, trifacets, buffer, verbose);
		}
		else
		{
			std::cout << "please specify correct output filename!" << std::endl;
		}
	}
	catch (...)
	{
		std::cout << "Encountered unknown error!" << std::endl;
		return false;
	}

	return true;
}

bool IsoSurfGen_FeatureAware(
	float start_x, float end_x,
	float start_y, float end_y,
	float start_z, float end_z,
	unsigned int gridsize,
	ImplicitFunc<float> *impfunc,
	const std::vector<int> &cubes, const std::vector<int> &cubetypes,
	const std::unordered_map<int, float> &point_value_map,
	const char *outobjfile,
	float iso_value,
	float feature_angle,
	const int buffer,
	bool verbose)
{
	ImplicitInterface<VectorType> mfunc(impfunc);
	ImplicitGrid<VectorType> grid(mfunc,								 // implicit
								  VectorType(start_x, start_y, start_z), // origin
								  VectorType(end_x - start_x, 0, 0),	 // x-axis
								  VectorType(0, end_y - start_y, 0),	 // y-axis
								  VectorType(0, 0, end_z - start_z),	 // z-axis
								  gridsize, gridsize, gridsize);		 // resolution

	std::vector<int> cache_cube(cubes.size() * 2);
	for (size_t i = 0; i < cubes.size(); i++)
	{
		cache_cube[2 * i] = cubes[i];
		cache_cube[2 * i + 1] = cubetypes[i];
	}
	MyMesh mesh;
	advanced_octree_extended_marching_cubes<MyMesh, ImplicitGrid<VectorType>, float>(grid, mesh, iso_value, feature_angle, cache_cube, point_value_map, buffer, verbose);

	std::string filename(outobjfile);
	std::string extname = "";
	if (filename.find_last_of(".") != std::string::npos)
		extname = filename.substr(filename.find_last_of(".") + 1);
	if (extname == "obj")
	{
		write_mesh(mesh, outobjfile);
	}
	else if (extname == "ply")
	{
		OpenMesh::IO::Options op;
		op += IO::Options::Binary;
		op += IO::Options::Swap;
		write_mesh(mesh, outobjfile, op);
	}
	return true;
}