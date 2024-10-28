#ifndef _TOOLS_MESHIO_H
#define _TOOLS_MESHIO_H

#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>
#include "../include/mesh_data_struct.h"

#define DEBUG_INFO() printf("File:%s, Line:%d, Function:%s\n", __FILE__, __LINE__, __FUNCTION__)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
namespace TiGER
{
	// Futher: Attribute should contain these data structure: double,int, vector(double,3), tensor(double,9), and no quantity constrain should be in there
	struct Mesh
	{
		Eigen::MatrixXd Vertex;
		Eigen::MatrixXi Topo;
		Eigen::MatrixXi Masks;
		Eigen::MatrixXd f_normal;
		Eigen::MatrixXd v_normal;
	};

	namespace MESHIO
	{
        static int readVTK_WYF_HW(std::string filename,
                                  std::vector<std::array<double, 3>> &V,
                                  std::vector<std::array<int, 3>> &F,
                                  std::vector<std::array<int, 4>> &T,
                                  std::vector<int> &F_geo,
                                  std::vector<int> &T_geo);
		static int readVTK(std::string filename, Mesh &mesh, std::string mark_pattern = "");
		static int readEPS(std::string filename, int &cou, std::map<int, double> &mpd, std::map<int, std::vector<int>> &mpi);
        static int readVTK(std::string filename, VolumeMesh& mesh, std::string mark_pattern);
		static int readVTK_newer_version(std::string filename, Mesh &mesh, std::string mark_pattern);
		static int readMESH(std::string filename, Mesh &mesh);
		static int readPLS(std::string filename, Mesh &mesh);
		static int readPLY(std::string filename, Mesh &mesh); // TODO
		static int readOBJ(std::string filename, Mesh &mesh);
		static int readVTK_toVolumemesh(std::string filename, TiGER::VolumeMesh &mesh, std::string mark_pattern = "");
		static int readFacet(std::string filename, Mesh &mesh);
		static int readTetgen(std::string nodefilename, std::string elefilename, Mesh &mesh);
		static int readSTL(std::string filename, Mesh &mesh);

		static int writeVTK(std::string filename, const Mesh &mesh, std::string mark_pattern = "");
		static int writeVTK(std::string filename, const TiGER::SurfaceMesh& mesh,std::string mark_pattern = "");
		static int writeVTK(std::string filename, const TiGER::VolumeMesh &mesh, std::string mark_pattern = "");
		static int writeVTK(std::string filename, const TiGER::LineMesh &mesh, std::string mark_pattern = " ");
		static int writeEpsVTK(std::string filename, const Mesh &mesh, int &cou, std::map<int, double> &mpd, std::map<int, std::vector<int>> &mpi, std::string mark_pattern = "");
		static int writeMESH(std::string filename, const Mesh &mesh);
		static int writePLY(std::string filename, const Mesh &mesh);
		static int writePLS(std::string filename, const Mesh &mesh);
		static int writeFacet(std::string filename, const Mesh &mesh);
		static int writeOBJ(std::string filename, const Mesh &mesh);
		static int writeStlIn(std::string filename, const Mesh &mesh);
		static int writeLineVTK(std::string filename,
						 Eigen::MatrixXd &V,
						 std::vector<std::array<int, 2>> &Lines);
	}

}
#include"meshIO.cpp"
#endif
