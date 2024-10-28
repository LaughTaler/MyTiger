/**
 * @author :wangyifei
*/
#pragma once

#include <vector>
#include <array>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include "pair_hash.h"
#include <Eigen/Dense>



//考虑了三角形四边形混合网格及非二边流形的情况,faceid从1开始
void facet_classfication(TiGER::SurfaceMesh& SurfaceMesh,double angle = 30){
	SurfaceMesh.regions.assign(SurfaceMesh.tris.size()+SurfaceMesh.quads.size(),-1);
	std::vector<std::array<int,4>> elements;
	for(int i=0;i<SurfaceMesh.tris.size();i++){
		std::array<int,4> temp_array={SurfaceMesh.tris[i][0],SurfaceMesh.tris[i][1],SurfaceMesh.tris[i][2],-1};
		elements.emplace_back(temp_array);
	}
	for(int i=0;i<SurfaceMesh.quads.size();i++){
		elements.emplace_back(SurfaceMesh.quads[i]);
	}
	std::unordered_map<std::pair<int,int>,std::vector<int>,pair_hash> edge_elements;
	for(int i=0;i<elements.size();i++){
		std::array<int,4> id = elements[i];
		if(id[3]==-1){
			edge_elements[make_pair_in_order(id[0],id[1])].emplace_back(i);
			edge_elements[make_pair_in_order(id[1],id[2])].emplace_back(i);
			edge_elements[make_pair_in_order(id[2],id[0])].emplace_back(i);
		}
		else{
			edge_elements[make_pair_in_order(id[0],id[1])].emplace_back(i);
			edge_elements[make_pair_in_order(id[1],id[2])].emplace_back(i);
			edge_elements[make_pair_in_order(id[2],id[3])].emplace_back(i);
			edge_elements[make_pair_in_order(id[3],id[0])].emplace_back(i);
		}
	}
	int current_faceid = 1;
	for(int i=0;i<SurfaceMesh.regions.size();i++){
		if(SurfaceMesh.regions[i]<0){
			std::queue<int> Q;
			Q.push(i);
			SurfaceMesh.regions[i]=current_faceid++;
			while(!Q.empty()){
				int current_element = Q.front();
				Q.pop();
				std::array<int,4> vertex_id = elements[current_element];
				std::vector<std::pair<int,int>> temp_edges;
				if(vertex_id[3]==-1){
					temp_edges.emplace_back(make_pair_in_order(vertex_id[0],vertex_id[1]));
					temp_edges.emplace_back(make_pair_in_order(vertex_id[1],vertex_id[2]));
					temp_edges.emplace_back(make_pair_in_order(vertex_id[2],vertex_id[0]));
				}
				else{
					temp_edges.emplace_back(make_pair_in_order(vertex_id[0],vertex_id[1]));
					temp_edges.emplace_back(make_pair_in_order(vertex_id[1],vertex_id[2]));
					temp_edges.emplace_back(make_pair_in_order(vertex_id[2],vertex_id[3]));
					temp_edges.emplace_back(make_pair_in_order(vertex_id[3],vertex_id[0]));
				}
				for(auto temp_edge : temp_edges){
					if(edge_elements[temp_edge].size()==2){
						int temp_element;
						if(current_element!=edge_elements[temp_edge][0]) temp_element=edge_elements[temp_edge][0];
						else temp_element=edge_elements[temp_edge][1];
						if(SurfaceMesh.regions[temp_element]>=0) continue;
						int id_v_edge[2]={temp_edge.first,temp_edge.second};
						int id_v_ele[2];
						int id_ele[2]={current_element,temp_element};
						for(int j=0;j<2;j++){
							for(int k=0;k<4;k++){
								int idk=elements[id_ele[j]][k];
								if(idk!=-1&&idk!=id_v_edge[0]&&idk!=id_v_edge[1]){
									id_v_ele[j]=idk;
									break;
								}
							}
						}
						Eigen::Vector3d v0,v1,ve0,ve1;
						int idnow=id_v_ele[0];
						for(int i=0;i<3;i++) v0(i)=SurfaceMesh.coords[idnow][i];
						idnow=id_v_ele[1];
						for(int i=0;i<3;i++) v1(i)=SurfaceMesh.coords[idnow][i];
						idnow=id_v_edge[0];
						for(int i=0;i<3;i++) ve0(i)=SurfaceMesh.coords[idnow][i];
						idnow=id_v_edge[1];
						for(int i=0;i<3;i++) ve1(i)=SurfaceMesh.coords[idnow][i];
						Eigen::Vector3d n0,n1;
						n0=((v0-ve0).cross(ve1-ve0)).normalized();
						n1=((ve1-ve0).cross(v1-ve0)).normalized();
						if(n0.dot(n1)>=cos(angle*3.14159265358979632/180.0)){
							Q.push(temp_element);
							SurfaceMesh.regions[temp_element]=SurfaceMesh.regions[current_element];
						}
					}
				}
			}
		}
	}

}

using namespace TiGER
{
    void findFacetMarkID(TiGER::SurfaceMesh & mesh_in_without_label, TiGER::SurfaceMesh & mesh_in_with_label);
} // namespace TiGER
