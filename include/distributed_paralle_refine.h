#include "mesh_data_struct.h"
#include "alias.h"
namespace TiGER
{
    /**
     * @brief 并行细分网格生成（待完善） 
     * 
     * @param volmesh_in 
     * @param surfmesh 
     * @param args 
     * @param volmesh_out 
     */
    void HexGrid_ParallelRefine(
        const VolumeMesh& volmesh_in,
        const SurfaceMesh& surfmesh,
        const std::string args,
        ParallelVolumeMesh& volmesh_out
    );
} // namespace TiGER
