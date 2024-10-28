#include <vector>
#include <unordered_map>
#include "../include/mesh_data_struct.h"
#include "mesh_repair.h"
#include "alias.h"
namespace TiGER {
    void merge_volume_mesh(
        const std::vector<VolumeMesh>& volume_in,
        VolumeMesh& volume_out){
        // 将所有mesh合并
        std::vector<std::array<double,3>> point_global;
        for(int i=0;i<volume_in.size();i++){
            for(auto point : volume_in[i].coords){
                point_global.push_back(point);
                // volume_out.coords.push_back(point);
            }
        }
        int count = 0;
        std::vector<std::array<int,4>> tet_global;
        std::vector<std::array<int,6>> prism_global;
        std::vector<std::array<int,5>> pyramid_global;
        std::vector<std::array<int,8>> hex_global;
        for(int i=0;i<volume_in.size();i++){
            for(auto tet : volume_in[i].tetras){
                std::array<int,4> new_tet = {
                    tet[0]+count,
                    tet[1]+count,
                    tet[2]+count,
                    tet[3]+count};
                tet_global.push_back(new_tet);
                // volume_out.tetras.push_back(new_tet);
            }

             for (auto tet : volume_in[i].prisms) {
                std::array<int, 6> new_tet = {tet[0] + count, tet[1] + count,
                                              tet[2] + count, tet[3] + count, tet[4] + count, tet[5] + count
                };
                prism_global.push_back(new_tet);
                // volume_out.tetras.push_back(new_tet);
            }
             for (auto tet : volume_in[i].pyramids) {
                std::array<int, 5> new_tet = {tet[0] + count, tet[1] + count,
                                              tet[2] + count, tet[3] + count,
                                              tet[4] + count};
                pyramid_global.push_back(new_tet);
                // volume_out.tetras.push_back(new_tet);
             }

               for (auto tet : volume_in[i].hexs) {
                std::array<int, 8> new_hex = {tet[0] + count, tet[1] + count,
                                              tet[2] + count, tet[3] + count,
                                              tet[4] + count,tet[5] + count,
                                              tet[6] + count,tet[7] + count};
                hex_global.push_back(new_hex);
                // volume_out.tetras.push_back(new_tet);
             }


            count += volume_in[i].coords.size();
        }
        // return;
        // 去重
        double eps = 1e-4;
        std::vector<int> bcj_pre(point_global.size(),-1);
        TiGER::repair::Tool_detectPointDuplicates(point_global,eps,bcj_pre);
        std::unordered_map<int,int> point_index;
        int index = 0;
        for(int i=0;i<bcj_pre.size();i++){
            if(bcj_pre[i] == -1){ // 将非重复点加入
                volume_out.coords.push_back(point_global[i]);
                point_index[i] = index;
                index++;
            }
        }
        for (int i = 0; i < tet_global.size(); i++) {
            std::array<int, 4> tet;
            for (int j = 0; j < 4; j++) {
                if (bcj_pre[tet_global[i][j]] != -1)
                  tet[j] = point_index[bcj_pre[tet_global[i][j]]];
                else
                  tet[j] = point_index[tet_global[i][j]];
            }
            volume_out.tetras.push_back(tet);
        }
        for (int i = 0; i < prism_global.size(); i++) {
            std::array<int, 6> prism;
            for (int j = 0; j < 6; j++) {
                if (bcj_pre[prism_global[i][j]] != -1)
                  prism[j] = point_index[bcj_pre[prism_global[i][j]]];
                else
                  prism[j] = point_index[prism_global[i][j]];
            }
            volume_out.prisms.push_back(prism);
        }
        for (int i = 0; i < pyramid_global.size(); i++) {

            std::array<int, 5> pyramid;
            for (int j = 0; j < 5; j++) {
                if (bcj_pre[pyramid_global[i][j]] != -1)
                    pyramid[j] = point_index[bcj_pre[pyramid_global[i][j]]];
                else
                    pyramid[j] = point_index[pyramid_global[i][j]];
            }
            volume_out.pyramids.push_back(pyramid);

        }
         for (int i = 0; i < hex_global.size(); i++) {

            std::array<int, 8> hex;
            for (int j = 0; j < 5; j++) {
                if (bcj_pre[hex_global[i][j]] != -1)
                    hex[j] = point_index[bcj_pre[hex_global[i][j]]];
                else
                    hex[j] = point_index[hex_global[i][j]];
            }
            volume_out.hexs.push_back(hex);

        }
        return;
        // 
        // double ave_eps = 0;
        // int loop_count = 0;
        // for(const auto& volume : volume_in){
        //     Eigen::MatrixXi topo;
        //     TiGER::list_to_matrix<4>(volume.tetras, topo);
        //     Eigen::MatrixXi loop;
        //     igl::boundary_loop(topo, loop);

        //     for (int k = 0; k < loop.rows() - 1; k++) {
        //         ave_eps += (Eigen::RowVector3d(volume.coords[loop(k, 0)].data()) -
        //                     Eigen::RowVector3d(volume.coords[loop(k + 1, 0)].data()))
        //                         .norm();
        //     }
        //     if (loop.rows() > 1) {
        //         ave_eps /= loop.rows();
        //         loop_count += loop.rows();
        //     }
        // }
        // if (loop_count > 0) {
        //     ave_eps /= loop_count;
        // }

        // std::map<std::array<double,3>,int> vertexIndexMap;
        // int newindex = 0;
        // for(const auto& mesh : volume_in){
        //     for(const auto& coord : mesh.coords){
        //         bool found = false;
        //         int existingIndex = -1;
        //         for(const auto& item : vertexIndexMap){
        //             if (std::sqrt(std::pow(item.first[0] - coord[0], 2) +
        //                 std::pow(item.first[1] - coord[1], 2) +
        //                 std::pow(item.first[2] - coord[2], 2)) < eps) {
        //                 found = true;
        //                 existingIndex = item.second;
        //                 break;
        //             }
        //         }
        //         if (!found) {
        //             vertexIndexMap[coord] = newIndex;
        //             output.coords.push_back(coord);
        //             existingIndex = newIndex++;
        //         }
        //     }
        // }
        // int volumeid = 0;
        // for(const auto& mesh : volume_in){
        //     for(const auto& tet : mesh.tetras){
        //         std::array<int,4> newtet;
        //         for(int i=0;i<4;i++){
        //             newtet[i] = vertexIndexMap[mesh.coords[tet[i]]];
        //         }
        //         volume_out.tetras.push_back(newtet);
        //         volume_out.regions.push_back(volumeid);
        //     }
        //     volumeid++;
        // }
        // std::cout << "Total unique vertices after merging: " << output.coords.size()
        //             << std::endl;
    }
}