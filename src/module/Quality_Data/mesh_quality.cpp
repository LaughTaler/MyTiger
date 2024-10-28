#include <mesh_quality.h>
#include <verdict.h>
#include <cmath>
#include <map>
#include <iostream>
#include <algorithm>
#include "alias.h"
using namespace std;


namespace TiGER
{

	TriMetricVals** pTriMetric;
	TetMetricVals** pTetMetric;
	HexMetricVals** pHexMetric;
	PrismMetricVals** pPrismMetric;
	PyramidMetricVals** pPyraMetric;
	std::vector<std::pair<int, int>> adjancent_s;
	std::vector<std::pair<int, int>> adjancent_v;

	static void connectivity_s(const TiGER::SurfaceMesh& surface_mesh)
	{
		int nstri = surface_mesh.tris.size();
		std::map<int, std::map<std::vector<int>, int>> hash_s;
		int stritris[3][2] = { {0,1},{0,2},{1,2} };


		for (int i = 0; i < nstri; i++) {
			for (int j = 0; j < 3; j++) {
				vector<int> index;
				for (int k = 0; k < 2; k++) {
					index.push_back(surface_mesh.tris[i][stritris[j][k]]);
				}
				sort(index.begin(), index.end());
				if (hash_s[index[0]].find(index) == hash_s[index[0]].end()) {
					hash_s[index[0]][index] = i;
				}
				else {
					adjancent_s.push_back({ hash_s[index[0]][index] ,i });
				}
			}
		}
	}

	static void connectivity_v(const TiGER::VolumeMesh& volume_mesh)
	{
		int ntetra = volume_mesh.tetras.size();
		//int nhex = volume_mesh.hexs.size();
		int nprism = volume_mesh.prisms.size();
		int npyra = volume_mesh.pyramids.size();
		int numelem = ntetra /*+ nhex*/ + nprism + npyra;
		vector<vector<int>> ej(numelem);
		std::map<int, std::map<std::vector<int>, int>>hash_e;

		int tetratris[4][3] = { {0,1,2},{0,2,3},{1,2,3},{0,1,3} };
		int prismtris[5][4] = { {0,1,2,-1},{3,4,5,-1},{1,2,4,5},{0,1,3,4},{2,0,5,3} };
		int hextris[6][4] = { {0,1,2,3},{0,1,4,5},{1,2,6,5},{2,3,7,6},{0,3,4,7},{4,5,6,7} };
		int pyratris[5][5] = { {0,1,2,3},{0,3,4,-1},{0,4,1,-1},{1,2,4,-1},{2,3,4,-1} };


		for (int i = 0; i < ntetra; i++) {
			for (int j = 0; j < 4; j++) {
				vector<int> index;
				for (int k = 0; k < 3; k++) {
					index.push_back(volume_mesh.tetras[i][tetratris[j][k]]);
				}
				sort(index.begin(), index.end());
				if (hash_e[index[0]].find(index) == hash_e[index[0]].end()) {
					hash_e[index[0]][index] = i;
				}
				else {
					adjancent_v.push_back({ hash_e[index[0]][index] ,i });
				}
			}
		}
		
		for (int i = 0; i < nprism; i++) {
			for (int j = 0; j < 5; j++) {
				vector<int> index;
				for (int k = 0; k < 4; k++) {
					if (prismtris[j][k] >= 0)
						index.push_back(volume_mesh.prisms[i][prismtris[j][k]]);
				}
				sort(index.begin(), index.end());
				if (hash_e[index[0]].find(index) == hash_e[index[0]].end()) {
					hash_e[index[0]][index] = i + ntetra;
				}
				else {
					adjancent_v.push_back({ hash_e[index[0]][index] ,i + ntetra });
				}
			}
		}

		//for (int i = 0; i < nhex; i++) {
		//	for (int j = 0; j < 6; j++) {
		//		vector<int> index;
		//		for (int k = 0; k < 4; k++) {
		//			index.push_back(volume_mesh.hexs[i][hextris[j][k]]);
		//		}
		//		sort(index.begin(), index.end());
		//		if (hash_e[index[0]].find(index) == hash_e[index[0]].end()) {
		//			hash_e[index[0]][index] = i + ntetra + nprism;
		//		}
		//		else {
		//			adjancent_v.push_back({ hash_e[index[0]][index] ,i + ntetra + nprism });
		//		}
		//	}
		//}

		for (int i = 0; i < npyra; i++) {
			for (int j = 0; j < 5; j++) {
				vector<int> index;
				for (int k = 0; k < 5; k++) {
					if (pyratris[j][k] >= 0)
						index.push_back(volume_mesh.pyramids[i][pyratris[j][k]]);
				}
				sort(index.begin(), index.end());
				if (hash_e[index[0]].find(index) == hash_e[index[0]].end()) {
					hash_e[index[0]][index] = i + ntetra + nprism/* + nhex*/;
				}
				else {
					adjancent_v.push_back({ hash_e[index[0]][index] ,i + ntetra + nprism /*+ nhex */});
				}
			}
		}

	}

	static void CalculateSurfRatio(QualityResult& result) {
		for (int i = 0; i < adjancent_s.size(); i++) {
			int sur1 = adjancent_s[i].first;
			int sur2 = adjancent_s[i].second;

			double s1 = pTriMetric[sur1]->area;
			double s2 = pTriMetric[sur2]->area;

			double arear = std::max(s1 / s2, s2 / s1);

			if (i == 0) {
				result.min_value = arear;
				result.max_value = arear;
				result.ave_value = 0.0;
				result.max_index = sur1;
				result.min_index = sur2;
				result.max_index_2 = sur2;
				result.min_index_2 = sur1;
			}

			if (result.min_value > arear) {
				result.min_value = arear;
				result.min_index = sur1;
				result.min_index_2 = sur2;
			}

			if (result.max_value < arear) {
				result.max_value = arear;
				result.max_index = sur1;
				result.max_index_2 = sur2;
			}

			result.ave_value += arear;
		}

		result.ave_value /= adjancent_s.size();

	}

	static void CalculateVolumeRatio(const TiGER::VolumeMesh& volume_mesh,QualityResult& result) {
		int ntetra = volume_mesh.tetras.size();
		//int nhex = volume_mesh.hexs.size();
		int nprism = volume_mesh.prisms.size();
		int npyra = volume_mesh.pyramids.size();
		int numelem = ntetra /*+ nhex*/ + nprism + npyra;
		for (int i = 0; i < adjancent_v.size(); i++) {
			int vol1 = adjancent_v[i].first;
			int vol2 = adjancent_v[i].second;

			double v1 = 0.0;
			double v2 = 0.0;

			if (vol1 < ntetra) {
				v1 = pTetMetric[vol1]->volume;
			}
			else if (vol1 >= ntetra && vol1 < ntetra + nprism) {
				vol1 -= ntetra;
				v1 = pPrismMetric[vol1]->volume;
			}
			else if (vol1 >= ntetra + nprism && vol1 < ntetra + nprism /*+ nhex*/) {
				vol1 = vol1 - ntetra - nprism;
				v1 = pHexMetric[vol1]->volume;
			}
			else {
				vol1 = vol1 - ntetra - nprism/* - nhex*/;
				v1 = pPyraMetric[vol1]->volume;
			}

			if (vol2 < ntetra) {
				v2 = pTetMetric[vol2]->volume;
			}
			else if (vol2 >= ntetra && vol2 < ntetra + nprism) {
				vol2 -= ntetra;
				v2 = pPrismMetric[vol2]->volume;
			}
			else if (vol2 >= ntetra + nprism && vol2 < ntetra + nprism /*+ nhex*/) {
				vol2 = vol2 - ntetra - nprism;
				v2 = pHexMetric[vol2]->volume;
			}
			else {
				vol2 = vol2 - ntetra - nprism /*- nhex*/;
				v2 = pPyraMetric[vol2]->volume;
			}

			double volumer = std::max(v1 / v2, v2 / v1);

			if (i == 0) {
				result.min_value = volumer;
				result.max_value = volumer;
				result.ave_value = 0.0;
				result.max_index = vol1;
				result.max_index_2 = vol2;
				result.min_index = vol1;
				result.min_index_2 = vol2;
			}

			if (result.min_value > volumer) {
				result.min_value = volumer; 
				result.min_index = vol1;
				result.min_index_2 = vol2;
			}

			if (result.max_value < volumer) {
				result.max_value = volumer;
				result.max_index = vol1;
				result.max_index_2 = vol2;
			}

			result.ave_value += volumer;

		}

		result.ave_value /= adjancent_v.size();

	}

	void Tool_SurfQuality(
		const TiGER::SurfaceMesh& surface_mesh,
		const SurfMeshQualityType& type,
		std::vector<QualityResult>& quality)
	{
		connectivity_s(surface_mesh);
		int nstri = surface_mesh.tris.size();
		unsigned long metrics_flag;
		QualityResult result;

		if (nstri > 0)
		{
			pTriMetric = new TriMetricVals * [nstri];
			for (int i = 0; i < nstri; i++)
				pTriMetric[i] = new TriMetricVals;
		}

		if (nstri > 0)
		{
			metrics_flag = 0;
			metrics_flag += V_TRI_ALL;
			double coord[3][3];

			for (int i = 0; i < nstri; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					coord[j][0] = surface_mesh.coords[surface_mesh.tris[i][j]][0];
					coord[j][1] = surface_mesh.coords[surface_mesh.tris[i][j]][1];
					coord[j][2] = surface_mesh.coords[surface_mesh.tris[i][j]][2];
				}
				v_tri_quality(3, coord, metrics_flag, pTriMetric[i]);
			}
			switch (type)
			{
				//TRI_ASPECT
			case TRI_ASPECT:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->aspect;
						result.max_value = pTriMetric[i]->aspect;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->aspect < result.min_value)
					{
						result.min_value = pTriMetric[i]->aspect;
						result.min_index = i;
					}

					if (pTriMetric[i]->aspect > result.max_value)
					{
						result.max_value = pTriMetric[i]->aspect;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->aspect / nstri;
				}
				quality.push_back(result);
				break;
				//TRI_AREA
			case TRI_AREA:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->area;
						result.max_value = pTriMetric[i]->area;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->area < result.min_value)
					{
						result.min_value = pTriMetric[i]->area;
						result.min_index = i;
					}

					if (pTriMetric[i]->area > result.max_value)
					{
						result.max_value = pTriMetric[i]->area;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->area / nstri;
				}
				quality.push_back(result);
				break;
				//TRI_ANGLE
			case TRI_ANGLE:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->minimum_angle;
						result.max_value = pTriMetric[i]->maximum_angle;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->minimum_angle < result.min_value)
					{
						result.min_value = pTriMetric[i]->minimum_angle;
						result.min_index = i;
					}

					if (pTriMetric[i]->maximum_angle > result.max_value)
					{
						result.max_value = pTriMetric[i]->maximum_angle;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->minimum_angle / nstri + pTriMetric[i]->maximum_angle / nstri;
					result.ave_value /= 2;
				}
				quality.push_back(result);
				break;
				//TRI_CONDITION
			case TRI_CONDITION:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->condition;
						result.max_value = pTriMetric[i]->condition;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->condition < result.min_value)
					{
						result.min_value = pTriMetric[i]->condition;
						result.min_index = i;
					}

					if (pTriMetric[i]->condition > result.max_value)
					{
						result.max_value = pTriMetric[i]->condition;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->condition / nstri;
				}
				quality.push_back(result);
				break;
				//TRI_SCALED_JACOBIAN
			case TRI_SCALED_JACOBIAN:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->scaled_jacobian;
						result.max_value = pTriMetric[i]->scaled_jacobian;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->scaled_jacobian < result.min_value)
					{
						result.min_value = pTriMetric[i]->scaled_jacobian;
						result.min_index = i;
					}

					if (pTriMetric[i]->scaled_jacobian > result.max_value)
					{
						result.max_value = pTriMetric[i]->scaled_jacobian;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->scaled_jacobian / nstri;
				}
				quality.push_back(result);
				break;
				//TRI_SHAPE
			case TRI_SHAPE:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->shape;
						result.max_value = pTriMetric[i]->shape;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->shape < result.min_value)
					{
						result.min_value = pTriMetric[i]->shape;
						result.min_index = i;
					}

					if (pTriMetric[i]->shape > result.max_value)
					{
						result.max_value = pTriMetric[i]->shape;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->shape / nstri;
				}
				quality.push_back(result);
				break;
				//TRI_RELATIVE_SIZE_SQUARED
			case TRI_RELATIVE_SIZE_SQUARED:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->relative_size_squared;
						result.max_value = pTriMetric[i]->relative_size_squared;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->relative_size_squared < result.min_value)
					{
						result.min_value = pTriMetric[i]->relative_size_squared;
						result.min_index = i;
					}

					if (pTriMetric[i]->relative_size_squared > result.max_value)
					{
						result.max_value = pTriMetric[i]->relative_size_squared;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->relative_size_squared / nstri;
				}
				quality.push_back(result);
				break;
				//TRI_SHAPE_AND_SIZE
			case TRI_SHAPE_AND_SIZE:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->shape_and_size;
						result.max_value = pTriMetric[i]->shape_and_size;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->shape_and_size < result.min_value)
					{
						result.min_value = pTriMetric[i]->shape_and_size;
						result.min_index = i;
					}

					if (pTriMetric[i]->shape_and_size > result.max_value)
					{
						result.max_value = pTriMetric[i]->shape_and_size;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->shape_and_size / nstri;
				}
				quality.push_back(result);
				break;
				//TRI_DISTORTION
			case TRI_DISTORTION:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->distortion;
						result.max_value = pTriMetric[i]->distortion;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->distortion < result.min_value)
					{
						result.min_value = pTriMetric[i]->distortion;
						result.min_index = i;
					}

					if (pTriMetric[i]->distortion > result.max_value)
					{
						result.max_value = pTriMetric[i]->distortion;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->distortion / nstri;
				}
				quality.push_back(result);
				break;
				//TRI_EDGE_RATIO
			case TRI_EDGE_RATIO:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->edge_ratio;
						result.max_value = pTriMetric[i]->edge_ratio;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->edge_ratio < result.min_value)
					{
						result.min_value = pTriMetric[i]->edge_ratio;
						result.min_index = i;
					}

					if (pTriMetric[i]->edge_ratio > result.max_value)
					{
						result.max_value = pTriMetric[i]->edge_ratio;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->edge_ratio / nstri;
				}
				quality.push_back(result);
				break;
				//TRI_SIZE_CHANGE
			case TRI_SIZE_CHANGE:
				break;
				//TRI_EQUIANGLE_SKEWNESS
			case TRI_EQUIANGLE_SKEWNESS:
				for (int i = 0; i < nstri; i++)
				{
					if (i == 0)
					{
						result.min_value = pTriMetric[i]->eqangle_skew;
						result.max_value = pTriMetric[i]->eqangle_skew;
						result.ave_value = 0.0;
						result.min_index = result.max_index = 0;
					}

					if (pTriMetric[i]->eqangle_skew < result.min_value)
					{
						result.min_value = pTriMetric[i]->eqangle_skew;
						result.min_index = i;
					}

					if (pTriMetric[i]->eqangle_skew > result.max_value)
					{
						result.max_value = pTriMetric[i]->eqangle_skew;
						result.max_index = i;
					}
					result.ave_value += pTriMetric[i]->eqangle_skew / nstri;
				}
				quality.push_back(result);
				break;
			case TRI_AREA_RATIO:
				CalculateSurfRatio(result);
				quality.push_back(result);
				break;
			}
		}

		if (pTriMetric)
		{
			for (int i = 0; i < nstri; i++)
			{
				delete[] pTriMetric[i];
				pTriMetric[i] = nullptr;
			}
			delete[] pTriMetric;
			pTriMetric = nullptr;
		}
	};

	void Tool_VolQuality(
		const TiGER::VolumeMesh& volume_mesh,
		const VolMeshQualityType& type,
		std::vector<QualityResult>& quality)
	{
		connectivity_v(volume_mesh);
		int ntetra = volume_mesh.tetras.size();
		//int nhex = volume_mesh.hexs.size();
		int nprism = volume_mesh.prisms.size();
		int npyra = volume_mesh.pyramids.size();

		unsigned long metrics_flag;
		QualityResult result;

		if (ntetra > 0)
		{
			pTetMetric = new TetMetricVals * [ntetra];
			for (int i = 0; i < ntetra; i++) {
				pTetMetric[i] = new TetMetricVals;
			}
		}

		//if (nhex > 0)
		//{
		//	pHexMetric = new HexMetricVals * [nhex];
		//	for (int i = 0; i < nhex; i++) {
		//		pHexMetric[i] = new HexMetricVals;
		//	}
		//}

		if (nprism > 0)
		{
			pPrismMetric = new PrismMetricVals * [nprism];
			for (int i = 0; i < nprism; i++) {
				pPrismMetric[i] = new PrismMetricVals;
			}
		}

		if (npyra > 0)
		{
			pPyraMetric = new PyramidMetricVals * [npyra];
			for (int i = 0; i < npyra; i++) {
				pPyraMetric[i] = new PyramidMetricVals;
			}
		}

		if (ntetra > 0)
		{
			metrics_flag = 0;
			metrics_flag += V_TET_ALL;
			double coord[4][3];

			for (int i = 0; i < ntetra; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					coord[j][0] = volume_mesh.coords[volume_mesh.tetras[i][j]][0];
					coord[j][1] = volume_mesh.coords[volume_mesh.tetras[i][j]][1];
					coord[j][2] = volume_mesh.coords[volume_mesh.tetras[i][j]][2];
				}
				v_tet_quality(4, coord, metrics_flag, pTetMetric[i]);
			}
		}

		//if (nhex > 0)
		//{
		//	metrics_flag = 0;
		//	metrics_flag += V_HEX_ALL;
		//	double coord[8][3];

		//	for (int i = 0; i < nhex; i++)
		//	{
		//		for (int j = 0; j < 8; j++)
		//		{
		//			coord[j][0] = volume_mesh.coords[volume_mesh.hexs[i][j]][0];
		//			coord[j][1] = volume_mesh.coords[volume_mesh.hexs[i][j]][1];
		//			coord[j][2] = volume_mesh.coords[volume_mesh.hexs[i][j]][2];
		//		}
		//		v_hex_quality(8, coord, metrics_flag, pHexMetric[i]);
		//	}
		//}

		if (nprism > 0)
		{
			metrics_flag = 0;
			metrics_flag += V_WEDGE_ALL;
			double coord[6][3];

			for (int i = 0; i < nprism; i++)
			{
				for (int j = 0; j < 6; j++)
				{
					coord[j][0] = volume_mesh.coords[volume_mesh.prisms[i][j]][0];
					coord[j][1] = volume_mesh.coords[volume_mesh.prisms[i][j]][1];
					coord[j][2] = volume_mesh.coords[volume_mesh.prisms[i][j]][2];
				}
				v_prism_quality(6, coord, metrics_flag, pPrismMetric[i]);
			}
		}

		if (npyra > 0)
		{
			metrics_flag = 0;
			metrics_flag += V_PYRAMID_ALL;
			double coord[5][3];

			for (int i = 0; i < npyra; i++)
			{
				for (int j = 0; j < 5; j++)
				{
					coord[j][0] = volume_mesh.coords[volume_mesh.pyramids[i][j]][0];
					coord[j][1] = volume_mesh.coords[volume_mesh.pyramids[i][j]][1];
					coord[j][2] = volume_mesh.coords[volume_mesh.pyramids[i][j]][2];
				}
				v_pyramid_quality(5, coord, metrics_flag, pPyraMetric[i]);
			}
		}

		switch (type)
		{
			//TET_ASPECT_BETA
		case TET_ASPECT_BETA:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->aspect_beta;
					result.max_value = pTetMetric[i]->aspect_beta;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->aspect_beta < result.min_value)
				{
					result.min_value = pTetMetric[i]->aspect_beta;
					result.min_index = i;
				}

				if (pTetMetric[i]->aspect_beta > result.max_value)
				{
					result.max_value = pTetMetric[i]->aspect_beta;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->aspect_beta / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_ASPECT_GAMMA
		case TET_ASPECT_GAMMA:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->aspect_gamma;
					result.max_value = pTetMetric[i]->aspect_gamma;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->aspect_gamma < result.min_value)
				{
					result.min_value = pTetMetric[i]->aspect_gamma;
					result.min_index = i;
				}

				if (pTetMetric[i]->aspect_gamma > result.max_value)
				{
					result.max_value = pTetMetric[i]->aspect_gamma;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->aspect_gamma / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_VOLUME
		case TET_VOLUME:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->volume;
					result.max_value = pTetMetric[i]->volume;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->volume < result.min_value)
				{
					result.min_value = pTetMetric[i]->volume;
					result.min_index = i;
				}

				if (pTetMetric[i]->volume > result.max_value)
				{
					result.max_value = pTetMetric[i]->volume;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->volume / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_CONDITION
		case TET_CONDITION:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->condition;
					result.max_value = pTetMetric[i]->condition;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->condition < result.min_value)
				{
					result.min_value = pTetMetric[i]->condition;
					result.min_index = i;
				}

				if (pTetMetric[i]->condition > result.max_value)
				{
					result.max_value = pTetMetric[i]->condition;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->condition / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_JACOBIAN
		case TET_JACOBIAN:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->jacobian;
					result.max_value = pTetMetric[i]->jacobian;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->jacobian < result.min_value)
				{
					result.min_value = pTetMetric[i]->jacobian;
					result.min_index = i;
				}

				if (pTetMetric[i]->jacobian > result.max_value)
				{
					result.max_value = pTetMetric[i]->jacobian;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->jacobian / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_SCALED_JACOBIAN
		case TET_SCALED_JACOBIAN:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->scaled_jacobian;
					result.max_value = pTetMetric[i]->scaled_jacobian;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->scaled_jacobian < result.min_value)
				{
					result.min_value = pTetMetric[i]->scaled_jacobian;
					result.min_index = i;
				}

				if (pTetMetric[i]->scaled_jacobian > result.max_value)
				{
					result.max_value = pTetMetric[i]->scaled_jacobian;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->scaled_jacobian / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_SHAPE
		case TET_SHAPE:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->shape;
					result.max_value = pTetMetric[i]->shape;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->shape < result.min_value)
				{
					result.min_value = pTetMetric[i]->shape;
					result.min_index = i;
				}

				if (pTetMetric[i]->shape > result.max_value)
				{
					result.max_value = pTetMetric[i]->shape;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->shape / ntetra;
			}
			quality.push_back(result);
			break;

			//TET_RELATIVE_SIZE
		case TET_RELATIVE_SIZE:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->relative_size_squared;
					result.max_value = pTetMetric[i]->relative_size_squared;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->relative_size_squared < result.min_value)
				{
					result.min_value = pTetMetric[i]->relative_size_squared;
					result.min_index = i;
				}

				if (pTetMetric[i]->relative_size_squared > result.max_value)
				{
					result.max_value = pTetMetric[i]->relative_size_squared;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->relative_size_squared / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_SHAPE_AND_SIZE
		case TET_SHAPE_AND_SIZE:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->shape_and_size;
					result.max_value = pTetMetric[i]->shape_and_size;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->shape_and_size < result.min_value)
				{
					result.min_value = pTetMetric[i]->shape_and_size;
					result.min_index = i;
				}

				if (pTetMetric[i]->shape_and_size > result.max_value)
				{
					result.max_value = pTetMetric[i]->shape_and_size;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->shape_and_size / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_DISTORTION
		case TET_DISTORTION:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->distortion;
					result.max_value = pTetMetric[i]->distortion;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->distortion < result.min_value)
				{
					result.min_value = pTetMetric[i]->distortion;
					result.min_index = i;
				}

				if (pTetMetric[i]->distortion > result.max_value)
				{
					result.max_value = pTetMetric[i]->distortion;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->distortion / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_EDGE_RATIO
		case TET_EDGE_RATIO:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->edge_ratio;
					result.max_value = pTetMetric[i]->edge_ratio;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->edge_ratio < result.min_value)
				{
					result.min_value = pTetMetric[i]->edge_ratio;
					result.min_index = i;
				}

				if (pTetMetric[i]->edge_ratio > result.max_value)
				{
					result.max_value = pTetMetric[i]->edge_ratio;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->edge_ratio / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_MINIMUM_DIHEDRAL
		case TET_MINIMUM_DIHEDRAL:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->minimum_dihedral;
					result.max_value = pTetMetric[i]->minimum_dihedral;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->minimum_dihedral < result.min_value)
				{
					result.min_value = pTetMetric[i]->minimum_dihedral;
					result.min_index = i;
				}

				if (pTetMetric[i]->minimum_dihedral > result.max_value)
				{
					result.max_value = pTetMetric[i]->minimum_dihedral;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->minimum_dihedral / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_MAXIMUM_DIHEDRAL
		case TET_MAXIMUM_DIHEDRAL:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->maximum_dihedral;
					result.max_value = pTetMetric[i]->maximum_dihedral;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->maximum_dihedral < result.min_value)
				{
					result.min_value = pTetMetric[i]->maximum_dihedral;
					result.min_index = i;
				}

				if (pTetMetric[i]->maximum_dihedral > result.max_value)
				{
					result.max_value = pTetMetric[i]->maximum_dihedral;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->maximum_dihedral / ntetra;
			}
			quality.push_back(result);
			break;
			//TET_SIZE_CHANGE
		case TET_SIZE_CHANGE:
			break;
			//TET_EQUIANGLE_SKEWNESS
		case TET_EQUIANGLE_SKEWNESS:
			for (int i = 0; i < ntetra; i++)
			{
				if (i == 0)
				{
					result.min_value = pTetMetric[i]->eqangle_skew;
					result.max_value = pTetMetric[i]->eqangle_skew;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pTetMetric[i]->eqangle_skew < result.min_value)
				{
					result.min_value = pTetMetric[i]->eqangle_skew;
					result.min_index = i;
				}

				if (pTetMetric[i]->eqangle_skew > result.max_value)
				{
					result.max_value = pTetMetric[i]->eqangle_skew;
					result.max_index = i;
				}
				result.ave_value += pTetMetric[i]->eqangle_skew / ntetra;
			}
			quality.push_back(result);
			break;
			//HEX_VOLUME
		case HEX_VOLUME:
			break;
			//for (int i = 0; i < nhex; i++)
			//{
			//	if (i == 0)
			//	{
			//		result.min_value = pHexMetric[i]->volume;
			//		result.max_value = pHexMetric[i]->volume;
			//		result.ave_value = 0.0;
			//		result.min_index = result.max_index = 0;
			//	}

			//	if (pHexMetric[i]->volume < result.min_value)
			//	{
			//		result.min_value = pHexMetric[i]->volume;
			//		result.min_index = i;
			//	}

			//	if (pHexMetric[i]->volume > result.max_value)
			//	{
			//		result.max_value = pHexMetric[i]->volume;
			//		result.max_index = i;
			//	}
			//	result.ave_value += pHexMetric[i]->volume / nhex;
			//}
			//quality.push_back(result);
			//PRISM_VOLUME
		case PRISM_VOLUME:
			for (int i = 0; i < nprism; i++)
			{
				if (i == 0)
				{
					result.min_value = pPrismMetric[i]->volume;
					result.max_value = pPrismMetric[i]->volume;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pPrismMetric[i]->volume < result.min_value)
				{
					result.min_value = pPrismMetric[i]->volume;
					result.min_index = i;
				}

				if (pPrismMetric[i]->volume > result.max_value)
				{
					result.max_value = pPrismMetric[i]->volume;
					result.max_index = i;
				}
				result.ave_value += pPrismMetric[i]->volume / nprism;
			}
			quality.push_back(result);
			break;
			//PRISM_EQUIANGLE_SKEWNESS
		case PRISM_EQUIANGLE_SKEWNESS:
			for (int i = 0; i < nprism; i++)
			{
				if (i == 0)
				{
					result.min_value = pPrismMetric[i]->eqangle_skew;
					result.max_value = pPrismMetric[i]->eqangle_skew;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pPrismMetric[i]->eqangle_skew < result.min_value)
				{
					result.min_value = pPrismMetric[i]->eqangle_skew;
					result.min_index = i;
				}

				if (pPrismMetric[i]->eqangle_skew > result.max_value)
				{
					result.max_value = pPrismMetric[i]->eqangle_skew;
					result.max_index = i;
				}
				result.ave_value += pPrismMetric[i]->eqangle_skew / nprism;
			}
			quality.push_back(result);
			break;
			//PRISM_CENTROID_SKEWNESS
		case PRISM_CENTROID_SKEWNESS:
			for (int i = 0; i < nprism; i++)
			{
				if (i == 0)
				{
					result.min_value = pPrismMetric[i]->centroid_skew;
					result.max_value = pPrismMetric[i]->centroid_skew;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pPrismMetric[i]->centroid_skew < result.min_value)
				{
					result.min_value = pPrismMetric[i]->centroid_skew;
					result.min_index = i;
				}

				if (pPrismMetric[i]->centroid_skew > result.max_value)
				{
					result.max_value = pPrismMetric[i]->centroid_skew;
					result.max_index = i;
				}
				result.ave_value += pPrismMetric[i]->centroid_skew / nprism;
			}
			quality.push_back(result);
			break;
			//PYRAMID_VOLUME
		case PYRAMID_VOLUME:
			for (int i = 0; i < npyra; i++)
			{
				if (i == 0)
				{
					result.min_value = pPyraMetric[i]->volume;
					result.max_value = pPyraMetric[i]->volume;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pPyraMetric[i]->volume < result.min_value)
				{
					result.min_value = pPyraMetric[i]->volume;
					result.min_index = i;
				}

				if (pPyraMetric[i]->volume > result.max_value)
				{
					result.max_value = pPyraMetric[i]->volume;
					result.max_index = i;
				}
				result.ave_value += pPyraMetric[i]->volume / npyra;
			}
			quality.push_back(result);
			break;
			//PYRAMID_EQUIANGLE_SKEWNESS
		case PYRAMID_EQUIANGLE_SKEWNESS:
			for (int i = 0; i < npyra; i++)
			{
				if (i == 0)
				{
					result.min_value = pPyraMetric[i]->eqangle_skew;
					result.max_value = pPyraMetric[i]->eqangle_skew;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pPyraMetric[i]->eqangle_skew < result.min_value)
				{
					result.min_value = pPyraMetric[i]->eqangle_skew;
					result.min_index = i;
				}

				if (pPyraMetric[i]->eqangle_skew > result.max_value)
				{
					result.max_value = pPyraMetric[i]->eqangle_skew;
					result.max_index = i;
				}
				result.ave_value += pPyraMetric[i]->eqangle_skew / npyra;
			}
			quality.push_back(result);
			break;
			//PYRAMID_CENTROID_SKEWNESS
		case PYRAMID_CENTROID_SKEWNESS:
			for (int i = 0; i < npyra; i++)
			{
				if (i == 0)
				{
					result.min_value = pPyraMetric[i]->centroid_skew;
					result.max_value = pPyraMetric[i]->centroid_skew;
					result.ave_value = 0.0;
					result.min_index = result.max_index = 0;
				}

				if (pPyraMetric[i]->centroid_skew < result.min_value)
				{
					result.min_value = pPyraMetric[i]->centroid_skew;
					result.min_index = i;
				}

				if (pPyraMetric[i]->centroid_skew > result.max_value)
				{
					result.max_value = pPyraMetric[i]->centroid_skew;
					result.max_index = i;
				}
				result.ave_value += pPyraMetric[i]->centroid_skew / npyra;
			}
			quality.push_back(result);
			break;
		case VOLUME_RATIO:
			CalculateVolumeRatio(volume_mesh, result);
			quality.push_back(result);
			break;
		}

		if (pTetMetric)
		{
			for (int i = 0; i < ntetra; i++)
			{
				delete[] pTetMetric[i];
				pTetMetric[i] = nullptr;
			}
			delete[] pTetMetric;
			pTetMetric = nullptr;
		}

		//if (pHexMetric)
		//{
		//	for (int i = 0; i < nhex; i++)
		//	{
		//		delete[] pHexMetric[i];
		//		pHexMetric[i] = nullptr;
		//	}
		//	delete[] pHexMetric;
		//	pHexMetric = nullptr;
		//}

		if (pPrismMetric)
		{
			for (int i = 0; i < nprism; i++)
			{
				delete[] pPrismMetric[i];
				pPrismMetric[i] = nullptr;
			}
			delete[] pPrismMetric;
			pPrismMetric = nullptr;
		}

		if (pPyraMetric)
		{
			for (int i = 0; i < npyra; i++)
			{
				delete[] pPyraMetric[i];
				pPyraMetric[i] = nullptr;
			}
			delete[] pPyraMetric;
			pPyraMetric = nullptr;
		}
	};

}; // namespace TiGER
