#include <vector>
#include<set>
#include<array>
#include<algorithm>
#include"../include/mesh_data_struct.h"
#include <Eigen/Dense>
#ifndef _SURFACE_FROM_VOLUME_H_
#define  _SURFACE_FROM_VOLUME_H_

namespace TiGER {

	bool surfaceFromVolume(
		const VolumeMesh& volume,
		SurfaceMesh& surface)
	{
		std::vector<std::array<int, 3>> tri_surface;
		std::vector<std::array<int, 4>> qua_surface;
		std::vector<std::set<std::array<int, 3>>> tri_hash(volume.coords.size());
		std::vector<std::set<std::array<int, 4>>> qua_hash(volume.coords.size());

		for (int i = 0; i < volume.tetras.size(); i++) {
			std::vector<std::vector<int>> tri_id_lst{ {0, 2, 1}, { 0,1,3 }, { 0,3,2 }, { 1,2,3 } };
			for (int j = 0; j < tri_id_lst.size();j++) {
				std::array<int, 3> temp;
				for (int k = 0; k < tri_id_lst[j].size();k++) {
					temp[k] = volume.tetras[i][tri_id_lst[j][k]];
				}
				tri_surface.push_back(temp);
				std::sort(temp.begin(), temp.end());
				if (tri_hash[temp[0]].find(temp) == tri_hash[temp[0]].end()) {
					tri_hash[temp[0]].insert(temp);
				}
				else {
					tri_hash[temp[0]].erase(temp);
				}
			}
		}
		for (int i = 0; i < volume.prisms.size(); i++) {
			std::vector<std::vector<int>> tri_id_lst{ {0,1,2},{3,5,4} };
			for (int j = 0; j < tri_id_lst.size(); j++) {
				std::array<int, 3> temp;
				for (int k = 0; k < tri_id_lst[j].size(); k++) {
					temp[k] = volume.prisms[i][tri_id_lst[j][k]];
				}
				tri_surface.push_back(temp);
				std::sort(temp.begin(), temp.end());
				if (tri_hash[temp[0]].find(temp) == tri_hash[temp[0]].end()) {
					tri_hash[temp[0]].insert(temp);
				}
				else {
					tri_hash[temp[0]].erase(temp);
				}
			}
			std::vector<std::vector<int>> qua_id_lst{ {0,3,4,1},{1,4,5,2},{0,2,5,3} };
			for (int j = 0; j < qua_id_lst.size(); j++) {
				std::array<int, 4> temp;
				for (int k = 0; k < qua_id_lst[j].size(); k++) {
					temp[k] = volume.prisms[i][qua_id_lst[j][k]];
				}
				qua_surface.push_back(temp);
				std::sort(temp.begin(), temp.end());
				if (qua_hash[temp[0]].find(temp) == qua_hash[temp[0]].end()) {
					qua_hash[temp[0]].insert(temp);
				}
				else {
					qua_hash[temp[0]].erase(temp);
				}
			}
		}
		for (int i = 0; i < volume.pyramids.size(); i++) {
			std::vector<std::vector<int>> tri_id_lst{ {0,1,4},{1,2,4},{2,3,4},{3,0,4} };
			for (int j = 0; j < tri_id_lst.size(); j++) {
				std::array<int, 3> temp;
				for (int k = 0; k < tri_id_lst[j].size(); k++) {
					temp[k] = volume.pyramids[i][tri_id_lst[j][k]];
				}
				tri_surface.push_back(temp);
				std::sort(temp.begin(), temp.end());
				if (tri_hash[temp[0]].find(temp) == tri_hash[temp[0]].end()) {
					tri_hash[temp[0]].insert(temp);
				}
				else {
					tri_hash[temp[0]].erase(temp);
				}
			}
			std::vector<std::vector<int>> qua_id_lst{ {0,3,2,1} };
			for (int j = 0; j < qua_id_lst.size(); j++) {
				std::array<int, 4> temp;
				for (int k = 0; k < qua_id_lst[j].size(); k++) {
					temp[k] = volume.pyramids[i][qua_id_lst[j][k]];
				}
				qua_surface.push_back(temp);
				std::sort(temp.begin(), temp.end());
				if (qua_hash[temp[0]].find(temp) == qua_hash[temp[0]].end()) {
					qua_hash[temp[0]].insert(temp);
				}
				else {
					qua_hash[temp[0]].erase(temp);
				}
			}
		}

		for (int i = 0; i < tri_surface.size(); i++) {
			std::array<int, 3> temp = tri_surface[i];
			std::sort(temp.begin(), temp.end());
			if (tri_hash[temp[0]].find(temp) != tri_hash[temp[0]].end()) {
				surface.tris.push_back(temp);
			}
		}
		for (int i = 0; i < qua_surface.size(); i++) {
			std::array<int, 4> temp = qua_surface[i];
			std::sort(temp.begin(), temp.end());
			if (qua_hash[temp[0]].find(temp) != qua_hash[temp[0]].end()) {
				surface.quads.push_back(temp);
			}
		}
		surface.coords = volume.coords;
		return 0;
	}
}  // namespace TiGER

#endif