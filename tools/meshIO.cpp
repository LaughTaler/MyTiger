#include "./meshIO.h"
#include "stl_reader.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <time.h>
#include <array>

//#include "IO/Legacy/vtkPolyDataReader.h"
//#include "Common/Core/vtkSmartPointer.h"
//#include "Common/DataModel/vtkPolyData.h"
//#include "IO/Legacy/vtkUnstructuredGridReader.h"
//#include "IO/Legacy/vtkPolyDataWriter.h"

#define BUFFER_LENGTH 256

using namespace std;

namespace TiGER{
/**
 * @brief Seperate string origin by given a set of patterns.
 * 
 * @param origin 
 * @return std::vector<std::string> 
 */
static std::vector<std::string> seperate_string(std::string origin) {
    std::vector<std::string> result;
    stringstream ss(origin);
    while(ss >> origin) result.push_back(origin);
    return result;
}

int MESHIO::writeVTK(std::string filename, const TiGER::LineMesh& mesh, 
                            std::string mark_pattern) {
    
    std::ofstream vtkfile(filename);
    vtkfile << "# vtk DataFile Version 2.0" << std::endl;
    vtkfile << "Line Mesh" << std::endl;
    vtkfile << "ASCII" << std::endl;
    vtkfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    vtkfile << "POINTS " << mesh.coord.size() << " double" << std::endl;
    for (const auto& point : mesh.coord)
    {
        vtkfile << point[0] << " " << point[1] << " " << point[2] << std::endl;
    }

    vtkfile << "CELLS " << mesh.segments.size() << " " << mesh.segments.size() * 3 << std::endl;
    for (const auto& segment : mesh.segments)
    {
        vtkfile << "2 " << segment[0] << " " << segment[1] << std::endl;
    }

    vtkfile << "CELL_TYPES " << mesh.segments.size() << std::endl;
    for (size_t i = 0; i < mesh.segments.size(); ++i)
    {
        vtkfile << "3 ";
    }
    vtkfile << std::endl;
    vtkfile.close();
    return 0;
}
int MESHIO::writeVTK(std::string filename, const TiGER::SurfaceMesh& mesh,
    std::string mark_pattern) {
    FILE* f = fopen(filename.c_str(), "w");
    if (!f) {
        std::cout << "Write VTK file failed. - " << filename << std::endl;
        return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "TetWild Mesh\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(f, "POINTS %d double\n", mesh.coords.size());
    for (int i = 0; i < mesh.coords.size(); i++) {
        fprintf(f, "%.17g %.17g %.17g\n", mesh.coords[i][0], mesh.coords[i][1],
            mesh.coords[i][2]);
    }
    fprintf(f, "CELLS %d %d\n", mesh.tris.size() + mesh.quads.size(),mesh.tris.size()*4+mesh.quads.size()*5);
    for (int i = 0; i < mesh.tris.size(); i++) {
        fprintf(f, "%d ",3);
        for (int j = 0; j < 3; j++) {
            fprintf(f, "%d ", mesh.tris[i][j]);
        }
        fprintf(f, "\n");
    }
    for (int i = 0; i < mesh.quads.size(); i++) {
        fprintf(f, "%d ", 4);
        for (int j = 0; j < 4; j++) {
            fprintf(f, "%d ", mesh.quads[i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "CELL_TYPES %d\n",
        mesh.tris.size() + mesh.quads.size());
    for (int i = 0; i < mesh.tris.size(); i++) {
        fprintf(f, "5\n");
    }
    for (int i = 0; i < mesh.quads.size(); i++) {
        fprintf(f, "9\n");
    }

    if (mesh.attribute_int.size() == mesh.tris.size() + mesh.quads.size()) {
        fprintf(f, "CELL_DATA %d\n", mesh.tris.size() + mesh.quads.size());
        fprintf(f, "SCALARS surface_id int %d\n", 1);
        fprintf(f, "LOOKUP_TABLE default\n");
        for (int i = 0; i < mesh.tris.size() + mesh.quads.size(); i++) {
            fprintf(f, "%d", mesh.attribute_int[i]);
            fprintf(f, "\n");
        }
    }


  
    fclose(f);

}
int MESHIO::writeVTK(std::string filename, const TiGER::VolumeMesh& mesh,
                            std::string mark_pattern) {
    
    FILE* f = fopen(filename.c_str(), "w");
    if (!f) {
      std::cout << "Write VTK file failed. - " << filename << std::endl;
      return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "TetWild Mesh\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(f, "POINTS %d double\n", mesh.coords.size());
    for (int i = 0; i < mesh.coords.size(); i++) {
      fprintf(f, "%.17g %.17g %.17g\n", mesh.coords[i][0], mesh.coords[i][1],
              mesh.coords[i][2]);
    }
    fprintf(f, "CELLS %d %d\n", mesh.prisms.size()+mesh.pyramids.size()+mesh.tetras.size(), 7*mesh.prisms.size()+6*mesh.pyramids.size()+5*mesh.tetras.size() );
    for (int i = 0; i < mesh.prisms.size(); i++) {
      fprintf(f, "%d ", 6);
      for (int j = 0; j < 6; j++) {
        fprintf(f, "%d ", mesh.prisms[i][j]);
      }
      fprintf(f, "\n");
    }
    for (int i = 0; i < mesh.pyramids.size(); i++) {
      fprintf(f, "%d ", 5);
      for (int j = 0; j < 5; j++) {
        fprintf(f, "%d ", mesh.pyramids[i][j]);
      }
      fprintf(f, "\n");
    }
    for (int i = 0; i < mesh.tetras.size(); i++) {
      fprintf(f, "%d ", 4);
      for (int j = 0; j < 4; j++) {
        fprintf(f, "%d ", mesh.tetras[i][j]);
      }
      fprintf(f, "\n");
    }
    fprintf(f, "CELL_TYPES %d\n",
            mesh.prisms.size() + mesh.pyramids.size() + mesh.tetras.size());
    for (int i = 0; i < mesh.prisms.size(); i++) {
      fprintf(f, "13\n");
    }
    for (int i = 0; i < mesh.pyramids.size(); i++) {
      fprintf(f, "14\n");
    }
    for (int i = 0; i < mesh.tetras.size(); i++) {
      fprintf(f, "10\n");
    }
    fclose(f);
   }
// The eps file format is 
/* eps 0.005          
 * id 1 5 7 9 100
 * eps 0.008
 * id 56 78 79
 * ...
 */
int TiGER::MESHIO::readEPS(std::string filename, int& cou, std::map<int, double>& mpd, std::map<int, vector<int>>& mpi) {
	std::ifstream eps_file;
	eps_file.open(filename);
	if(!eps_file.is_open()){
		std::cout << "No Such file. - " << filename << std::endl;
		return -1;
	}
	string line;
	string _;
	stringstream input;
    // int cou = 0;
    // std::map<int, double> mpd;
    // std::map<int, vector<int>> mpi;
	while(getline(eps_file, line)){
		input.clear();
		input.str(line);
		input >> _;
		if(_[0] == '#') continue;
		else if(_ == "eps" || _ == "EPS"){
            cou++;
			input >> mpd[cou];
		}
		else if(_ == "id" || _ == "ID"){
			int tmpId;
			while (input >> tmpId)
			{
				mpi[cou].push_back(tmpId);
			}
		}
		else{
			continue;
		}
	}
	return 0;
}
int TiGER::MESHIO::readVTK_newer_version(std::string filename, Mesh& mesh, std::string mark_pattern) {
    auto& M = mesh.Masks;
    auto& V = mesh.Vertex;
    auto& T = mesh.Topo;
    int nPoints = 0;
    int nFacets = 0;
    std::ifstream vtk_file;
    vtk_file.open(filename);
    if (!vtk_file.is_open()) {
        std::cout << "No such file. - " << filename << std::endl;
        return -1;
    }
    std::string vtk_type_str = "CONNECTIVITY";
    char buffer[BUFFER_LENGTH];
    while (!vtk_file.eof()) {
        vtk_file.getline(buffer, BUFFER_LENGTH);
        std::string line = (std::string)buffer;
        if (vtk_file.peek() == EOF) break;
        if (line.length() < 2 || buffer[0] == '#')
            continue;

        if (line.find("POINTS ") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            nPoints = stoi(words[1]);
            V.resize(nPoints, 3);
            for (int i = 0; i < nPoints; i++) {
                vtk_file >> V(i, 0) >> V(i, 1) >> V(i, 2);
                
            }
        }
        if (line.find("CELLS") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            nFacets = stoi(words[1])-1;
            T.resize(nFacets, 3);
        }
        if (line.find("CONNECTIVITY") != std::string::npos) {
            
            string trash;
            //vtk_file >> trash;
            for (int i = 0; i < nFacets; i++) {
                vtk_file >> T(i,0) >> T(i, 1) >> T(i, 2);
            }
          //  std::cout << T(0, 0) << std::endl;
        }
    }
  //  std::cout << T(0, 0) << std::endl;
    vtk_file.close();
    return 1;
}


int TiGER::MESHIO::readVTK_WYF_HW(std::string filename,
                           std::vector<std::array<double, 3>>& V,
                           std::vector<std::array<int, 3>>& F,
                           std::vector<std::array<int, 4>>& T,
                           std::vector<int>& F_geo, std::vector<int>& T_geo) {
  int nPoints = 0;
  int nFacets = 0;
  int nTets = 0;
  std::ifstream vtk_file;
  vtk_file.open(filename);
  if (!vtk_file.is_open()) {
    std::cout << "No such file: " << filename << std::endl;
    return 0;
  }
  std::string vtk_type_str = "POLYDATA ";
  char buffer[BUFFER_LENGTH];
  while (!vtk_file.eof()) {
    vtk_file.getline(buffer, BUFFER_LENGTH);
    std::string line = (std::string)buffer;
    if (line.length() < 2 || buffer[0] == '#') continue;
    if (line.find("DATASET") != std::string::npos) {
      std::vector<std::string> words = seperate_string(line);
      if (words[1] == "POLYDATA")
        vtk_type_str = "POLYGONS ";
      else if (words[1] == "UNSTRUCTURED_GRID")
        vtk_type_str = "CELLS ";
      else {
        std::cout
            << "The format of VTK file is illegal, No clear DATASET name. - "
            << filename << std::endl;
      }
    }
    if (line.find("POINTS ") != std::string::npos) {
      std::vector<std::string> words = seperate_string(line);
      nPoints = stoi(words[1]);
      V.resize(nPoints);
      for (int i = 0; i < nPoints; i++) {
        if (!(vtk_file >> V[i][0] >> V[i][1] >> V[i][2])) return 0;
      }
    }
    if (line.find(vtk_type_str) != std::string::npos) {
      std::vector<std::string> words = seperate_string(line);
      int a = stoi(words[1]);
      int b = stoi(words[2]);
      nTets = b - 4 * a;
      nFacets = a - nTets;
      F.resize(nFacets);
      for (int i = 0; i < nFacets; i++) {
        int nindices;
        vtk_file >> nindices;
        if (nindices != 3) return 0;
        vtk_file >> F[i][0] >> F[i][1] >> F[i][2];
      }
      T.resize(nTets);
      for (int i = 0; i < nTets; i++) {
        int nindices;
        vtk_file >> nindices;
        if (nindices != 4) return 0;
        vtk_file >> T[i][0] >> T[i][1] >> T[i][2] >> T[i][3];
      }
    }
    if (line.find("GeoId ") != std::string::npos) {
      vtk_file.getline(buffer, BUFFER_LENGTH);
      F_geo.resize(nFacets);
      T_geo.resize(nTets);
      for (int i = 0; i < nFacets; i++) {
        vtk_file >> F_geo[i];
      }
      for (int i = 0; i < nTets; i++) {
        vtk_file >> T_geo[i];
      }
    }
  }
  vtk_file.close();
  return 0;
}

int TiGER::MESHIO::readVTK(std::string filename, Mesh& mesh, std::string mark_pattern) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    int nPoints = 0;
    int nFacets = 0;
    std::ifstream vtk_file;
    vtk_file.open(filename);
    if(!vtk_file.is_open()) {
        std::cout << "No such file. - " << filename << std::endl;
        return -1;
    }
    std::string vtk_type_str = "POLYDATA ";
    char buffer[BUFFER_LENGTH];
    while(!vtk_file.eof()) {
        vtk_file.getline(buffer, BUFFER_LENGTH);
        std::string line = (std::string)buffer;
        if(line.length() < 2 )
            continue;
        if (line.find("Version") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            if (words[4][0] == '5') {
                vtk_file.close();
                return readVTK_newer_version(filename,mesh,mark_pattern);              
            }
        }
        if(line.find("DATASET") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            if(words[1] == "POLYDATA") 
                vtk_type_str = "POLYGONS ";
            else if(words[1] == "UNSTRUCTURED_GRID")
                vtk_type_str = "CELLS ";
            else {
                std::cout << "The format of VTK file is illegal, No clear DATASET name. - " << filename << std::endl;
            }
        }
        if(line.find("POINTS ") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            nPoints = stoi(words[1]);
            V.resize(nPoints, 3);
            for(int i = 0; i < nPoints; i++) {
                vtk_file.getline(buffer, BUFFER_LENGTH);
                words = seperate_string(std::string(buffer));
                for(int j = 0; j < 3; ++j){
                    V(i, j) = stod(words[j]);
                }
            }
        }
        if(line.find(vtk_type_str) != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            nFacets = stoi(words[1]);
            T.resize(nFacets, stoi(words[2]) / nFacets - 1);
            for(int i = 0; i < nFacets; i++) {
                vtk_file.getline(buffer, BUFFER_LENGTH);
                words = seperate_string(std::string(buffer));
                for(int j = 0; j < stoi(words[0]); j++) 
                    T(i, j) = stoi(words[j + 1]);
            }
            M.resize(nFacets, 1);
			for (int i = 0; i < nFacets; i++)
				M(i, 0) = 0;
        }
        if(line.find("CELL_DATA ") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            if(stoi(words[1]) != nFacets) {
                std::cout << "The number of CELL_DATA is not equal to number of cells. -" << filename;
                std::cout << "Ignore CELL_DATA" << std::endl;
                return 0;
            }
            if (vtk_type_str == "POLYGONS ")
            {
                vtk_file.getline(buffer, BUFFER_LENGTH);
                vtk_file.getline(buffer, BUFFER_LENGTH);
                std::string data_type = seperate_string(std::string(buffer))[0];
                if (data_type != mark_pattern)
                    continue;
                for (int i = 0; i < nFacets; i++){
                    vtk_file.getline(buffer, BUFFER_LENGTH);
                    int surface_id = stoi(std::string(buffer));
                    M(i, 0) = surface_id;
                }
            }
            else
            {
                vtk_file.getline(buffer, BUFFER_LENGTH);
                std::string data_type = seperate_string(std::string(buffer))[1];
                if (data_type != mark_pattern)
                    continue;
                vtk_file.getline(buffer, BUFFER_LENGTH);
                for (int i = 0; i < nFacets; i++)
                {
                    vtk_file.getline(buffer, BUFFER_LENGTH);
                    int surface_id = stoi(std::string(buffer));
                    M(i, 0) = surface_id;
                }
            }
        }
    }
    vtk_file.close();
    return 1;
}

int MESHIO::readVTK_toVolumemesh(std::string filename, TiGER::VolumeMesh &mesh, std::string mark_pattern){
    auto& coords = mesh.coords;
    auto& tetras = mesh.tetras;
    auto& prisms = mesh.prisms;
    auto& pyramids = mesh.pyramids;
    auto& polygons = mesh.polygons;
    auto& regoins = mesh.regions;
    int nPoints = 0;
    int nVolume = 0;
    std::ifstream vtk_file;
    vtk_file.open(filename);
    if(!vtk_file.is_open()) {
        std::cout << "No such file. - " << filename << std::endl;
        return -1;
    }
    std::string vtk_type_str = "POLYDATA ";
    char buffer[BUFFER_LENGTH];
    while(!vtk_file.eof()) {
        vtk_file.getline(buffer,BUFFER_LENGTH);
        std::string line = (std::string)buffer;
        if(line.length() < 2)
            continue;;
        if (line.find("Version") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            if (words[4][0] == '5') {
                vtk_file.close();
                std::cout << "Can't read - " << filename << std::endl;
                return -1;            
            }
        }
        if(line.find("DATASET") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            if(words[1] == "POLYDATA") 
                vtk_type_str = "POLYGONS ";
            else if(words[1] == "UNSTRUCTURED_GRID")
                vtk_type_str = "CELLS ";
            else {
                std::cout << "The format of VTK file is illegal, No clear DATASET name. - " << filename << std::endl;
            }
        }
        if(line.find("POINTS ") != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            nPoints = stoi(words[1]);
            coords.resize(nPoints);
            for(int i = 0; i < nPoints; i++) {
                vtk_file.getline(buffer, BUFFER_LENGTH);
                words = seperate_string(std::string(buffer));
                for(int j = 0; j < 3; ++j){
                    coords[i][j] = stod(words[j]);
                }
            }
        }
        if(line.find(vtk_type_str) != std::string::npos) {
            std::vector<std::string> words = seperate_string(line);
            nVolume = stoi(words[1]);
            for(int i = 0; i < nVolume; i++) {
                vtk_file.getline(buffer, BUFFER_LENGTH);
                words = seperate_string(std::string(buffer));
                if(stoi(words[0]) == 4){
                    std::array<int,4> tet;
                    for(int j = 0; j < stoi(words[0]); j++) 
                        tet[j] = stoi(words[j+1]);
                    tetras.push_back(tet);
                }
                else if(stoi(words[0]) == 5){
                    std::array<int,5> pyramid;
                    for(int j = 0; j < stoi(words[0]); j++) 
                        pyramid[j] = stoi(words[j+1]);
                    pyramids.push_back(pyramid);
                }
                else if(stoi(words[0]) == 6){
                    std::array<int,6> prism;
                    for(int j = 0; j < stoi(words[0]); j++) 
                        prism[j] = stoi(words[j+1]);
                    prisms.push_back(prism);
                }
            }
            regoins.resize(nVolume);
        }
    }
}

int MESHIO::readOBJ(std::string filename, Mesh &mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    ifstream objFile;
    objFile.open(filename);
    if(!objFile.is_open()) {
        std::cout << "No such file. - " << filename << std::endl;
        return -1;
    }
    char buffer[BUFFER_LENGTH];
    vector<vector<double>> plist;
    vector<vector<int>> flist;
    vector<int> mlist;
    int curMark = 0;
    while(!objFile.eof()) {
        objFile.getline(buffer, BUFFER_LENGTH);
        if(buffer[0] == 'v') {
            vector<string> words = seperate_string(string(buffer));
            if(words.size() < 4) {
                continue;
            }
            vector<double> coord;
            for(int ii = 1; ii < 4; ii++) {
                int wordLen = 0;
                while(wordLen < words[ii].length()) {
                    if(words[ii][wordLen] == '/') {
                        break;
                    }
                    wordLen++;
                }
                coord.push_back(stod(words[ii].substr(0, wordLen)));
            }
            plist.push_back(coord);
        }
        if(buffer[0] == 'f') {
            vector<string> words = seperate_string(string(buffer));
            if(words.size() < 4) {
                continue;
            }
            vector<int> facet;
            for(int ii = 1; ii < 4; ii++) {
                int wordLen = 0;
                while(wordLen < words[ii].length()) {
                    if(words[ii][wordLen] == '/') {
                        break;
                    }
                    wordLen++;
                }
                facet.push_back(stoi(words[ii].substr(0, wordLen)) - 1);
            }
            flist.push_back(facet);
            mlist.push_back(curMark);
        }
        if(buffer[0] == 'g') {
            ++curMark;
        }
    }

    V.resize(plist.size(), 3);
    for(int i = 0; i < plist.size(); i++) {
        for(int j = 0; j < 3; j++) {
            V(i, j) = plist[i][j];
        }
    }
    T.resize(flist.size(), 3);
    for(int i = 0; i < flist.size(); i++) {
        for(int j = 0; j < 3; j++) {
            T(i, j) = flist[i][j];
        }
    }
    M.resize(mlist.size(), 1);
    for(int i = 0; i < mlist.size(); i++) {
        M(i, 0) = mlist[i];
    }

    return 1;
}

int TiGER::MESHIO::writeVTK(std::string filename, const Mesh &mesh, std::string mark_pattern) {
    auto& M = mesh.Masks;
    auto& V = mesh.Vertex;
    auto& T = mesh.Topo;
    FILE *f = fopen(filename.c_str(), "w");
    if (!f) {
        std::cout << "Write VTK file failed. - " << filename << std::endl;
        return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "TetWild Mesh\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(f, "POINTS %d double\n", V.rows());
    for (int i = 0; i < V.rows(); i++) {
        fprintf(f, "%.17g %.17g %.17g\n", V(i, 0), V(i, 1), V(i, 2));
    }
    fprintf(f, "CELLS %d %d\n", T.rows(), T.rows() * (T.cols() + 1));
    for (int i = 0; i < T.rows(); i++) {
        fprintf(f, "%d ", T.cols());
        for (int j = 0; j < T.cols(); j++) {
            fprintf(f, "%d ", T(i, j));
        }
        fprintf(f, "\n");
    }
    fprintf(f, "CELL_TYPES %d\n", T.rows());
    int cellType = 0;
    if (T.cols() == 2) {
        cellType = 3;
    } else if (T.cols() == 3) {
        cellType = 5;
    } else if (T.cols() == 4) {
        cellType = 10;
    }
    else if (T.cols() == 6) {
        cellType = 13;
    }

    for (int i = 0; i < T.rows(); i++) {
        fprintf(f, "%d\n", cellType);
    }
    if (M.rows() != T.rows()) {
        fclose(f);
        return 1;
    }
    fprintf(f, "POINT_DATA %d\n", V.rows());
    fprintf(f, "SCALARS point_id int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < V.rows(); i++) {
        fprintf(f, "%d\n", i);
    }
    // fprintf(f, "CELL_DATA %d\n", T.rows());
    // fprintf(f, "SCALARS %s int %d\n", mark_pattern.c_str(), 1);
    // fprintf(f, "LOOKUP_TABLE default\n");
    // for (int i = 0; i < T.rows(); i++) {
    //         fprintf(f, "%d", i);
    //     fprintf(f, "\n");
    // }
    fprintf(f, "CELL_DATA %d\n", M.rows());
    fprintf(f, "SCALARS %s int %d\n", mark_pattern.c_str(), M.cols());
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.cols(); j++) {
            fprintf(f, "%d", M(i, j));
        }
        fprintf(f, "\n");
    }
    if (T.cols() == 3) {
        std::cout << "This have normal information.\n";
        fprintf(f, "NORMALS facet_normals double\n");
        for (int i = 0; i < T.rows(); i++) {
            Eigen::Vector3d ab = V.row(T(i, 1)) - V.row(T(i, 0));
            Eigen::Vector3d ac = V.row(T(i, 2)) - V.row(T(i, 0));
            Eigen::Vector3d t_nor = ab.cross(ac);
            t_nor.normalize();
            fprintf(f, "%.17g %.17g %.17g\n", t_nor.x(), t_nor.y(), t_nor.z());
        }
    }
    fclose(f);
    return 1;
}


// int MESHIO::writeVTK(std::string filename, const Mesh &mesh, std::string mark_pattern) {
// 	auto& M = mesh.Masks;
// 	auto& V = mesh.Vertex;
// 	auto& T = mesh.Topo;
//     std::ofstream f(filename);
//     if(!f.is_open()) {
//         std::cout << "Write VTK file failed. - " << filename << std::endl;
//         return -1;
//     }
//     std::cout << "Writing mesh to - " << filename << std::endl;
//     f.precision(std::numeric_limits<double>::digits10 + 1);
//     f << "# vtk DataFile Version 2.0" << std::endl;
//     f << "TetWild Mesh" << std::endl;
//     f << "ASCII" << std::endl;
//     f << "DATASET UNSTRUCTURED_GRID" << std::endl;
//     f << "POINTS " << V.rows() << " double" << std::endl;
//     for(int i = 0; i < V.rows(); i++)
//         f << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
//     f << "CELLS " << T.rows() << " " << T.rows() * (T.cols() + 1) << std::endl;
//     for(int i = 0; i < T.rows(); i++) {
//         f << T.cols() << " ";
//         for(int j = 0; j < T.cols(); j++)
//             f <<  T(i, j) << " ";
//         f << std::endl;
//     }
//     f << "CELL_TYPES " << T.rows() << std::endl;
//     int cellType = 0;
//     if(T.cols() == 2)
//         cellType = 3;
//     else if(T.cols() == 3)
//         cellType = 5;
//     else if(T.cols() == 4)
//         cellType = 10;
//     for(int i = 0; i < T.rows(); i++)
//         f << cellType << std::endl;
//     if(M.rows() != T.rows()) {
//         f.close();
//         return 1;
//     }
//     f << "POINT_DATA " << V.rows() << std::endl;
//     f << "SCALARS point_id int 1" << std::endl;
//     f << "LOOKUP_TABLE default" << std::endl;
//     for(int i = 0; i < V.rows(); i++) {
//         f << i << endl;
//     }
//     f << "CELL_DATA " << M.rows() << std::endl;
//     f << "SCALARS " << mark_pattern << " int " << M.cols() << std::endl;
//     f << "LOOKUP_TABLE default" << std::endl;
//     for(int i = 0; i < M.rows(); i++) {
//         for(int j = 0; j < M.cols(); j++)
//             f << M(i, j);
//         f << std::endl;
//     }
//     if(T.cols() == 3)
//     {
//         std::cout << "This have normal information.\n";
//         f << "NORMALS facet_normals double\n";
//         for(int i = 0; i < T.rows(); i++)
//         {
//             Eigen::Vector3d ab = V.row(T(i, 1)) - V.row(T(i, 0));
//             Eigen::Vector3d ac = V.row(T(i, 2)) - V.row(T(i, 0));
//             Eigen::Vector3d t_nor = ab.cross(ac);
//             t_nor.normalize();
//             f << t_nor.x() << " " << t_nor.y() << " " << t_nor.z() << std::endl;
//         }
//     }
//     f << std::endl;
//     f.close();
//     return 1;
// }

int TiGER::MESHIO::writeEpsVTK(std::string filename, const Mesh &mesh, int& cou,  std::map<int, double> &mpd, std::map<int, vector<int>> &mpi, std::string mark_pattern) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
	std::ofstream f(filename);
	if(!f.is_open()) {
		std::cout << "Write VTK file failed. - " << filename << std::endl;
		return -1;
	}
	std::cout << "Writing mesh to - " << filename << std::endl;
	f.precision(std::numeric_limits<double>::digits10 + 1);
	f << "# vtk DataFile Version 2.0" << std::endl;
	f << "TetWild Mesh" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;
	f << "POINTS " << V.rows() << " double" << std::endl;
	for(int i = 0; i < V.rows(); i++)
		f << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
	f << "CELLS " << V.rows() << " " << V.rows() * 2 << std::endl;
	for(int i = 0; i < V.rows(); i++) {
		f << 1 <<  " " << i ;
		f << std::endl;
	}
	f << "CELL_TYPES " << V.rows() << std::endl;
	int cellType = 1;
	for(int i = 0; i < V.rows(); i++)
		f << cellType << std::endl;

	f << "POINT_DATA " << V.rows() << std::endl;
	f << "SCALARS local_epsilon double 1" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;

	int curLoc = 0;

    vector<double> vec;
    vec.resize(V.rows(), -1);
    for(int i = 1; i <= cou; i++)
    {
        for(int id : mpi[i]){
            vec[id] = mpd[i];
        }
    }

	for(int i = 0; i < V.rows(); i++) {
		f << vec[i] << endl;
        if(vec[i] != -1){
            std::cout << i << " " << vec[i] << '\n';
        }
	}
	f << std::endl;
	f.close();
	return 1;
}

int TiGER::MESHIO::readMESH(std::string filename, Mesh& mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    std::ifstream mesh_file;
    mesh_file.open(filename);
    if(!mesh_file.is_open()) {
        std::cout << "No such file. - " << filename << std::endl;
        return -1;
    }
    int dimension = 3;
    int nPoints;
    int nFacets;
    char buffer[BUFFER_LENGTH];
    while(!mesh_file.eof()) {
        mesh_file.getline(buffer, BUFFER_LENGTH);
        std::string line = (std::string)buffer;
        if(line.length() < 2)
            continue;
        std::vector<std::string> words = seperate_string(line);
        if(words.empty())
            continue;
        if(line.find("Dimension") != std::string::npos) {
            dimension = stoi(words[1]);
            std::cout << "Reading mesh dimension - " << dimension << std::endl;
        }
        if(line == "Vertices") {
            words = seperate_string(line);
            if(words.size() > 1 || words[0].length() > 8)
                continue;
            line.clear();
            while(line.empty()) {
                mesh_file.getline(buffer, BUFFER_LENGTH);
                line = (std::string)buffer;
            }
            words = seperate_string(line);
            nPoints = std::stoi(words[0]);
            std::cout << "Number of points : " << nPoints << std::endl;
            V.resize(nPoints, dimension);
            int i = 0;
            while(i < nPoints) {
                mesh_file.getline(buffer, BUFFER_LENGTH);
                words = seperate_string(std::string(buffer));
                if(words.size() != (dimension + 1)) continue;
                for(int j = 0; j < dimension; j++)
                    V(i, j) = std::stod(words[j]);
                i++;
            }
        }
        if(line.find("Triangles") != std::string::npos) {
            line.clear();
            while(line.empty()) {
                mesh_file.getline(buffer, BUFFER_LENGTH);
                line = (std::string)buffer;
            }
            words = seperate_string(line);
            nFacets = stoi(words[0]);
            std::cout << "Number of facets : " << nFacets << std::endl;
            T.resize(nFacets, 3);
            M.resize(nFacets, 1);
            int i = 0;
            while(i < nFacets) {
                mesh_file.getline(buffer, BUFFER_LENGTH);
                words = seperate_string(std::string(buffer));
                if(words.size() < 4)
                {
                    std::cout << "Warning : The number of triangles element is not equal 4.\n";
                }
                for(int j = 0; j < 3; j++)
                    T(i, j) = std::stoi(words[j]) - 1;
                M(i, 0) = std::stoi(words[3]) - 1;
                i++;
            }
        }
    }
    mesh_file.close();
    return 1;
}

#include <cstdio>

int TiGER::MESHIO::writeMESH(std::string filename, const Mesh &mesh) {
    auto& M = mesh.Masks;
    auto& V = mesh.Vertex;
    auto& T = mesh.Topo;
    FILE* file = fopen(filename.c_str(), "w");
    if(!file) {
        std::cout << "Write MESH file failed. - " << filename << std::endl;
        return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    fprintf(file, "MeshVersionFormatted 1\n");
    fprintf(file, "Dimension %ld\n", V.cols());
    fprintf(file, "Vertices\n");
    fprintf(file, "%ld\n", V.rows());
    for(int i = 0; i < V.rows(); i++) {
        for(int j = 0; j < V.cols(); j++)
            fprintf(file, "%.15g ", V(i, j));
        fprintf(file, "%d\n", i + 1);
    }
    if(T.cols() == 3)
        fprintf(file, "Triangles\n");
    else if(T.cols() == 4)
        fprintf(file, "Tetrahedra\n");
    else {
        std::cout << "Unsupported format for .mesh file." << std::endl;
        fclose(file);
        return -1;
    }
    fprintf(file, "%ld\n", T.rows());
    for(int i = 0; i < T.rows(); i++) {
        for(int j = 0; j < T.cols(); j++)
            fprintf(file, "%d ", T(i, j) + 1);
        if(mesh.Masks.rows() == T.rows())
            fprintf(file, "%d\n", M(i, 0));
        else
            fprintf(file, "0\n");
    }
    fclose(file);
    return 1;
}


int TiGER::MESHIO::writePLY(std::string filename, const Mesh &mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    if(T.cols() != 3) {
        std::cout << "Unsupported format for .ply file." << std::endl;
        return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    std::ofstream plyfile;

    plyfile.open(filename);
	plyfile.precision(std::numeric_limits<double>::digits10 + 1);
    plyfile << "ply" << std::endl;
    plyfile << "format ascii 1.0" << std::endl;
    plyfile << "comment VTK generated PLY File" << std::endl;
    plyfile << "obj_info vtkPolyData points and polygons: vtk4.0" << std::endl;
    plyfile << "element vertex " << V.rows() << std::endl;
    plyfile << "property float x" << std::endl;
    plyfile << "property float y" << std::endl;
    plyfile << "property float z" << std::endl;
    plyfile << "element face " << T.rows() << std::endl;
    plyfile << "property list uchar int vertex_indices" << std::endl;
    plyfile << "end_header" << std::endl;
    for(int i = 0; i < V.rows(); i++)
        plyfile << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
    for(int i = 0; i < T.rows(); i++)
        plyfile << T.cols() << " " << T(i, 0) << " " << T(i, 1) << " " << T(i, 2) << std::endl;
    plyfile.close();
    return 1;
}

int TiGER::MESHIO::writePLS(std::string filename, const Mesh &mesh)
{
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    if(T.cols() != 3)
    {
        std::cout << "Unsupported format for .pls file." << std::endl;
        return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    std::ofstream plsfile;
    plsfile.open(filename);
    plsfile << T.rows() << " " << V.rows() << " " << "0 0 0 0\n";

	int minmask = std::numeric_limits<int>::max();
    Eigen::MatrixXd M_dup(T.rows(),1);
    if (M.rows() == 0) {
        minmask = 0;
    }
    else

	for (int i = 0; i < M.rows(); i++) {
		minmask = min(M(i, 0), minmask);
	}
	
	int diff = 1 - minmask;// min mask must be one in pls


	plsfile.precision(15);
    for(int i = 0; i < V.rows(); i++)
        plsfile << i + 1 << " " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;

	for (int i = 0; i < T.rows(); i++) {
		plsfile << i + 1 << " " << T(i, 0) + 1 << " " << T(i, 1) + 1 << " " << T(i, 2) + 1;
		if (M.rows()>=T.rows() && M.cols()>=1)
			plsfile << " " << M(i, 0);
		else
			plsfile << " " << 1;
		plsfile << std::endl;
	}

        
    plsfile.close();
    std::cout << "Finish\n";
    return 1;
}

int TiGER::MESHIO::readPLS(std::string filename, Mesh &mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo; 
	std::ifstream pls_file;
	pls_file.open(filename);
	if (!pls_file.is_open()) {
		std::cout << "No such file. - " << filename << std::endl;
		return -1;
	}
	int nPoints = 0;
	int nFacets = 0;

	std::string str;
	std::getline(pls_file, str);
	std::stringstream ss(str);
	ss >> nFacets >> nPoints;
	V.resize(nPoints,3);
	for (int i = 0; i < nPoints; i++) {
		int index;
		double point_coordinate[3];
		pls_file >> index >> point_coordinate[0] >> point_coordinate[1] >> point_coordinate[2];
		for(int k=0;k<3;k++)
		V(index-1,k)=point_coordinate[k];
	}
	M.resize(nFacets,1);
	T.resize(nFacets, 3);
	for (int i = 0; i < nFacets; i++) {
		int index;
		int tri[3];
		pls_file >> index >> tri[0] >> tri[1] >> tri[2] >> M(i,0);
		tri[0]--; tri[1]--; tri[2]--;
		for (int k = 0; k < 3; k++)
			T(i, k) = tri[k];
	}

	return 1;
}

int TiGER::MESHIO::readPLY(std::string filename, Mesh& mesh) {
    auto& M = mesh.Masks;
    auto& V = mesh.Vertex;
    auto& T = mesh.Topo;

    std::ifstream plyfile;
    plyfile.open(filename);
    if(!plyfile.is_open()) {
		std::cout << "No such file. - " << filename << std::endl;
        return -1;
    }

    int nPoints = 0;
    int nFacets = 0;
    
    string line;
    vector<string> words;
    while(!plyfile.eof()) {
        getline(plyfile, line);
        if(line.find("element") != std::string::npos) {
            words = seperate_string(line);
            if(words[1] == "vertex") {
                nPoints = stoi(words[2]);
                V.resize(nPoints, 3);
                cout << nPoints << endl;
            }
            else if(words[1] == "face") {
                nFacets = stoi(words[2]);
                T.resize(nFacets, 3);
                cout << nFacets << endl;
            }
        }
        if(line.find("end_header") != std::string::npos) {
            for(int i = 0; i < nPoints; i++) {
                getline(plyfile, line);
                words = seperate_string(line);
                for(int j = 0; j < 3; j++) {
                    V(i, j) = stod(words[j]);
                }
            }
            for(int i = 0; i < nFacets; i++) {
                getline(plyfile, line);
                stringstream ss(line);
                int _;
                ss >> _ >> T(i, 0) >> T(i, 1) >> T(i, 2);
            }
        }
    }

    plyfile.close();
    return 1;
}

int TiGER::MESHIO::writeFacet(std::string filename, const Mesh &mesh)
{
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
	std::cout << "Writing mesh to - " << filename << std::endl;
	std::ofstream facetfile;
	facetfile.open(filename);
	facetfile.precision(std::numeric_limits<double>::digits10 + 1);
	facetfile << "FACET FILE V3.0  exported from Meshconverter http://10.12.220.71/tools/meshconverter " << std::endl;
	facetfile << 1 << std::endl;
	facetfile << "Grid" << std::endl;
	facetfile << "0, 0.00 0.00 0.00 0.00" << endl;
	facetfile << V.rows() << endl;
	for (int i = 0; i < V.rows(); i++)
		facetfile << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
	facetfile << 1 << endl;
	facetfile << "Triangles" << endl;
	facetfile << T.rows() << " 3" << endl;
	
	for (int i = 0; i < T.rows(); i++) {
		facetfile << " " << T(i, 0) + 1 << " " << T(i, 1) + 1 << " " << T(i, 2) + 1 << " 0 ";
		if (M.rows() >= T.rows())
			facetfile << M(i, 0);
		else
			facetfile << 0;
		facetfile<< " " << i + 1 << std::endl;
	}
	return 0;
}

int TiGER::MESHIO::readFacet(string filename, Mesh &mesh) {
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo; 
	std::ifstream facetfile;
    facetfile.open(filename);
    if(!facetfile.is_open()) {
        return -1;
    }
    int nPoints = 0;
    int nFacets = 0;

    string line;
    while(!facetfile.eof()) {
        getline(facetfile, line);
        if(line[0] == '#') {
            continue;
        }
        if(line.find("Grid") != std::string::npos) {
            getline(facetfile, line);   // ignore line like "0, 0.00, 0.00, 0.00, 0.00"
            getline(facetfile, line);
            nPoints = stoi(line);
            V.resize(nPoints, 3);
            vector<string> words;
            for(int i = 0; i < nPoints; i++) {
                getline(facetfile, line);
                words = seperate_string(line);
                if(words.size() != 3) {
                    break;
                }
                for(int j = 0; j < 3; j++) {
                    V(i, j) = stod(words[j]);
                }
            }
        }
        if(line.find("Triangles") != std::string::npos) {
            getline(facetfile, line);
            vector<string> words = seperate_string(line);
            nFacets = stoi(words[0]);
            M.resize(nFacets, 1);
            T.resize(nFacets, 3);
            for(int i = 0; i < nFacets; i++) {
                getline(facetfile, line);
                words = seperate_string(line);
                for(int j = 0; j < 3; j++) {
                    T(i, j) = stoi(words[j]) - 1;
                }
                M(i, 0) = stoi(words[4]);
            }
        }
    }
    facetfile.close();
    return 1;
}

int TiGER::MESHIO::writeOBJ(string filename, const Mesh &mesh) {
    // Facet group
	auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& F = mesh.Topo; 
	bool doGroup = (M.rows() == F.rows());
    vector<vector<int>> flist;
    for(int i = 0; i < F.rows(); i++) {
        vector<int> facet;
        facet.push_back(doGroup ? M(i, 0) : 0);
        for(int j = 0; j < F.cols(); j++) {
            facet.push_back(F(i, j));
        }
        flist.push_back(facet);
    }
    if(doGroup) {
        sort(flist.begin(), flist.end(), [](vector<int> A, vector<int> B){ return A[0] < B[0]; });
    }

	cout << "Writing mesh to - " << filename << endl;
	ofstream objFile;
	objFile.open(filename);
	objFile.precision(std::numeric_limits<double>::digits10 + 1);

    // Get current time.
    std::string export_time;
    char stime[BUFFER_LENGTH] = {0};
    time_t now_time;
    time(&now_time);
    strftime(stime,sizeof(stime),"%H:%M:%S",localtime(&now_time));
    export_time = stime + '\0';

    // Header 
    objFile << "# TIGER Mesh converter. (c) 2021." << endl;
    objFile << "# Created File: " << export_time << endl;
    objFile << "# " << endl;
    objFile << "# object default" << endl;
    objFile << "# " << endl;
    objFile << endl;

    // Write points
    for(int i = 0; i < V.rows(); i++) {
        objFile << "v";
        for(int j = 0; j < V.cols(); j++) {
            objFile << " " << V(i, j);
        }
        objFile << endl;
    }
    objFile << "# " << V.rows() << " vertices" << endl << endl;

    // Write facets with groups
    int curGroup = INT_MIN;
    for(int i = 0; i < flist.size(); i++) {
        if(flist[i][0] != curGroup) {
            curGroup = flist[i][0];
            objFile << "g " << curGroup << endl;
        }
        objFile << "f";
        for(int j = 1; j < flist[i].size(); j++) {
            objFile << " " << flist[i][j] + 1;
        }
        objFile << endl;
    }

    // Write facets 
    objFile << "# " << F.rows() << " faces" << endl << endl;

    objFile.close();

    return 1;
}

int TiGER::MESHIO::readTetgen(string nodefilename, string elefilename, Mesh &mesh) {
    auto& T = mesh.Topo;
    auto& V = mesh.Vertex;
    string line;
    vector<string> words;

    ifstream nodefile(nodefilename);
    if(!nodefile.is_open()) {
        std::cout << "No such file. - " << nodefilename << std::endl;
        return -1;
    }
    getline(nodefile, line);
    words = seperate_string(line);
    int nPoint = stoi(words[0]);
    int pointType = 3;
    V.resize(nPoint, pointType);
    while(!nodefile.eof()) {
        getline(nodefile, line);
        if(line.length() < 2 || line[0] == '#') {
            continue;
        }
        words = seperate_string(line);
        int pindex = stoi(words[0]) - 1;
        for(int j = 0; j < pointType; j++) {
            V(pindex, j) = stod(words[j + 1]);
        }
    }
    nodefile.close();

    if(elefilename.empty()) {
        return 1;
    }
    ifstream elefile(elefilename);
    if(!elefile.is_open()) {
        std::cout << "No such file - " << elefilename << std::endl;
        return -1;
    }
    getline(elefile, line);
    words = seperate_string(line);
    int nElement = stoi(words[0]);
    int elementType = stoi(words[1]);
    if(elementType == 1) elementType = 3;
    T.resize(nElement, elementType);
    while(!elefile.eof()) {
        getline(elefile, line);
        if(line.length() < 2 || line[0] == '#') {
            continue;
        }
        words = seperate_string(line);
        int eindex = stoi(words[0]) - 1;
        for(int j = 0; j < elementType; j++) {
            T(eindex, j) = stoi(words[j + 1]) - 1;
        }
    }
    elefile.close();

    return 1;
}

int TiGER::MESHIO::readSTL(std::string filename, Mesh &mesh){
    auto& M = mesh.Masks;
    auto& V = mesh.Vertex;
    auto& T = mesh.Topo;
    int nPoints = 0;
    int nFacets = 0;

    stl_reader::StlMesh<> stlmesh;
    int result = stlmesh.read_file(filename);
    if(!result){
        return -1;
    }

    nPoints = stlmesh.num_vrts(); // 点数量
    V.resize(nPoints, 3);

    for(int i = 0; i < nPoints; i++){
        V.row(i) << stlmesh.raw_coords()[i*3 + 0] , stlmesh.raw_coords()[i*3 + 1] ,  stlmesh.raw_coords()[i*3 + 2];
    }

    nFacets = stlmesh.num_tris(); // 三角面片数量
    T.resize(nFacets,3);
    for(int i = 0; i < nFacets; i++){
        for(int j = 0; j < 3 ; j++){
            T(i, j) = stlmesh.tri_corner_ind(i,j);
        }
    }

    M.resize(nFacets,1);  // 面id
    for (int i = 0; i < nFacets; i++)
		M(i, 0) = 0;  // 所有三角面片都在同一个面

    return 1;
}
int TiGER::MESHIO::writeStlIn(std::string filename, const Mesh &mesh)
{
    auto& M = mesh.Masks;
	auto& V = mesh.Vertex;
	auto& T = mesh.Topo;
    std::ofstream f(filename);
    if(!f.is_open()) {
        std::cout << "Write VTK file failed. - " << filename << std::endl;
        return -1;
    }
    std::cout << "Writing mesh to - " << filename << std::endl;
    f.precision(std::numeric_limits<double>::digits10 + 1);
    f << T.rows() << " " << V.rows() << " 0 0 0 0\n";
    for(int i = 0; i < V.rows(); i++)
        f << i + 1 << " " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
    for(int i = 0; i < T.rows(); ++i)
    {
        f << i + 1 << " " << T(i, 0) + 1 << " " <<
        T(i, 1) + 1 << " " << T(i, 2) + 1 << " " << M(i, 0) + 1 << "\n";
    }
    f.close();
}

int TiGER::MESHIO::writeLineVTK(std::string filename,
                 Eigen::MatrixXd &V,
                 std::vector<std::array<int, 2>> &Lines)
{
    FILE *vtkFile = fopen(filename.c_str(), "w");

    if (vtkFile == nullptr)
    {
        perror("Error opening VTK file");
        return 0;
    }

    // 写入VTK文件头部信息
    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "VTK output\n");
    fprintf(vtkFile, "ASCII\n");
    fprintf(vtkFile, "DATASET POLYDATA\n");

    // 写入顶点信息
    fprintf(vtkFile, "POINTS %ld double\n", V.rows());
    for (int i = 0; i < V.rows(); ++i)
    {
        fprintf(vtkFile, "%.6lf %.6lf %.6lf\n", V(i, 0), V(i, 1), V(i, 2));
    }

    // 写入线信息
    fprintf(vtkFile, "LINES %ld %ld\n", Lines.size(), 3 * Lines.size());
    for (const std::array<int, 2> &line : Lines)
    {
        fprintf(vtkFile, "2 %d %d\n", line[0], line[1]);
    }

    fclose(vtkFile);

    return 1;
}

}