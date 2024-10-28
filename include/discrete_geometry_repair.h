#include <vector>
#include <array>
#include <set>
#include "alias.h"
#include "./mesh_data_struct.h"
#include "./discrete_geometry_struct.h"
#include "./discrete_geometry_repair_parameters.h"

#if defined(_WIN32) || defined(_WIN64)
    #ifdef Tiger_EXPORTS
    #define Tiger_API __declspec(dllexport)
    #else
    #define Tiger_API __declspec(dllimport)
    #endif
#else
    #ifdef Tiger_EXPORTS
    #define Tiger_API __attribute__((visibility("default")))
    #else
    #define Tiger_API
    #endif
#endif

namespace TiGER{
    namespace discrete_geometry_repair{

        struct EdgeTriPair{
            std::array<int,2> edge;
            int tri;  // 
            int type;
            std::array<double, 3> point;
        };

        // 交线
        struct InterEdge{
            int startid;  //
            int endid;  
            int shareface;  // edge相连的面
            std::set<int> faceid;
        };
        

        /**
         * @brief 初始化 数据结构
         * 
         * @param[in] surf_in 输入三角形表面网格
         */
        Tiger_API void Tool_initializeData(const SurfaceMesh& surf_in);
        
        /**
         * @brief 面网格转离散几何
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[out] args 输入参数
         * @param[out] discrete_geometry_out 输出离散几何
         */
        Tiger_API void Tool_convertToDiscreteGeometry(
            const SurfaceMesh& surf_in,
            const RepairParameters& args,
            DiscreteGeometry& discrete_geometry_out
        );

        //************* 检测 *****************//

        /**
         * @brief 穿刺面检测，检测与另一个面的一条或多条边相交的面
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[out] edges_out 输出相交的边
         * @param[out] tris_out 输出被穿过的三角形表面网格
         */
        Tiger_API void TMeshSurf_detectSurfacePunctures(
            const SurfaceMesh& surf_in,
            // std::vector<DiscreteGeometry::Edge>& edges_out,  // 两个点id定一条边
            std::vector<EdgeTriPair>& edge_tri
            // std::vector<DiscreteGeometry::TriangleIndex>& tris_out   // 三角形面id 
        );

        /**
         * @brief 低质量面检测
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] args 输入参数 包括:
         *      operator 输入运算符  (低于 高于 介于 不介于)
         *      tolerance 输入阈值
         *      type 质量类型
         * @param[out] tris_out 输出低质量的三角形表面网格
         */
        Tiger_API void detectSufaceQuality(
            const SurfaceMesh& surf_in,
            const RepairParameters& args,
            std::vector<DiscreteGeometry::TriangleIndex>& tris_out
        );

        /**
         * @brief 接近面检测 
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] args 输入参数 包括:
         *      operator 输入运算符  (低于 高于 介于 不介于)
         *      tolerance 输入阈值
         * @param[out] tri_out 输出面邻近值小于阈值的面对
         */
        Tiger_API void TMeshSurf_identifyMeshAdjacency(
            const SurfaceMesh& surf_in,
            const RepairParameters& args,
            std::vector<DiscreteGeometry::Triangle_Triangle>& tri_out  // 两个面id 一对
        );

        /**
         * @brief 自由边检测，检测仅与一个面接合的边
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[out] edges_out 输出自由边
         */
        Tiger_API void TMeshSurf_identifyFreeEdges(
            const SurfaceMesh& surf_in,
            std::vector<DiscreteGeometry::Edge>& edges_out  // 两个点id定一条边
        );

        /**
         * @brief 多重边检测，检测三个或更多其他面共享的边
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[out] edges_out 输出多重边
         */
        Tiger_API void TMeshSurf_identifyMultipleEdges(
            const SurfaceMesh& surf_in,
            std::vector<DiscreteGeometry::Edge>& edges_out  // 两个点id定一条边
        );

        /**
         * @brief 多重点检测，检测非流形点
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[out] point_out 输出多重点
         */
        Tiger_API void TMeshSurf_identifyMultiplePoints(
            const SurfaceMesh& surf_in,
            std::vector<DiscreteGeometry::PointIndex>& points_out
        );

        //************* 修复 *****************//
        // 可能会改变拓扑信息

        /**
         * @brief 删除选定面
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] tris_delete 输入删掉三角面片的id
         * @param[out] surf_out 输出删除指定三角面的面网格
         */
        void TMeshSurf_removeTriangles(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::TriangleIndex>& tris_delete,  // 三角片id
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 选中三个节点创建三角面片
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] tri_points 输入三个节点
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_createTrianglesFromPoints(
            const DiscreteGeometry& surf_in,
            const std::array<DiscreteGeometry::PointIndex,3>& tri_points,
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 重叠节点（重叠到第一个节点）
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] points_in 输入一系列节点
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_identifyOverlappingPoints(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::PointIndex>& points_in,
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 分割选定边
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] edges_in 输入一系列边
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_splitMeshEdges(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::Edge>& edges_in ,
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 分割选定三角片
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] tris_in 输入一系列三角形面
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_splitMeshTriangles(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::TriangleIndex>& tris_in,
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 交换选定边
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] edge_in 输入一条边
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_flipMeshEdges(
            const DiscreteGeometry& surf_in,
            const DiscreteGeometry::Edge& edge_in ,
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 用三角片填充多边形
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] points_coords 输入点集合
         * @param[in] args 输入字符串控制参数 
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_fillPolygonWithTriangles(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::Point>& points_coords,
            const RepairParameters& args,
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 自动修复表面错误 （面网格质量、面接近值、自相交）
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_automateSurfaceRepair(
            const DiscreteGeometry& surf_in,
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 使用选定边填充孔
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] edges_in 输入一系列边
         * @param[in] args 输入字符串控制参数 （填充模式 开环选项 填充选项）
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_createTrianglesFromEdges(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::Edge>& edges_in,
            const RepairParameters& args,
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 缝合选定边
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] edges_in 输入一系列边
         * @param[in] args tolerance 容差
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_sutureEdges(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::Edge>& edges_in,
            const RepairParameters& args,
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 重构选定面
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] tris_in 输入一系列三角形面
         * @param[out] surf_out 输出三角形表面网格
         */
        void remesh(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::TriangleIndex>& tris_in, 
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 平滑选定/关联节点
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] points_in 输入一系列节点
         * @param[in] edges_in 输入一系列边
         * @param[in] tris_in 输入一系列三角形面
         * @param[in] args 输入字符串控制参数(保留特征边  约束平滑  迭代次数)
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_smoothMeshPoints(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::PointIndex>& points_in,
            const std::vector<DiscreteGeometry::Edge>& edges_in,
            const std::vector<DiscreteGeometry::TriangleIndex>& tris_in,
            const RepairParameters& args,
            DiscreteGeometry& surf_out
        );

        /**
         * @brief 合并相邻节点
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] points_in 输入一系列节点
         * @param[in] args tolerance 输入容差
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_mergeAdjacentMeshPoints(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::PointIndex>& points_in,
            const RepairParameters& args,
            DiscreteGeometry& surf_out
        );
        
        /**
         * @brief 合并重复面
         * 
         * @param[in] surf_in 输入三角形表面网格
         * @param[in] tris_in 输入一系列三角形面
         * @param[out] surf_out 输出三角形表面网格
         */
        void TMeshSurf_mergeDuplicateSurfaces(
            const DiscreteGeometry& surf_in,
            const std::vector<DiscreteGeometry::TriangleIndex>& tris_in,
            DiscreteGeometry& surf_out
        );


        //************* 修复 *****************//

        /**
         * @brief 相交面-面集合
         *     两个面集合求交  （需要提供特征边）
         *     两个属性
         *     各自面id
         * 
         * @param[in] surf_in_A 输入三角形表面网格集合
         * @param[in] surf_in_B 输入三角形表面网格集合
         * @param[in] args eps_bool_option 是否开启容差布尔  0 
         *                 bool_eps  容差 1e-6
         * @param[in] surf_out 输出三角形表面网格
         * @param[out] edges_out 输出交线
         */
        Tiger_API void intersectionFaceSets(
            const SurfaceMesh& surf_in_A,
            const SurfaceMesh& surf_in_B,
            const RepairParameters& args,
            SurfaceMesh& surf_out,
            std::vector<InterEdge>& edges_out  // 抽取交线
        );

        /**
         * @brief 相交面
         *   两个domain求交  两个domain给定初始面id
         *   返回surface面id的结构为     （surf_A id）(surf_B id) (new id)  各8位
         * 
         * @param[in] surf_in_A 输入三角形表面网格集合
         * @param[in] surf_in_B 输入三角形表面网格集合
         * @param[in] args eps_bool_option 是否开启容差布尔  0 
         *                  bool_eps  容差 1e-6
         * @param[in] surf_out 输出三角形表面网格
         * @param[out] edges_out 输出交线
         */
        // Tiger_API void intersectionFaces(
        //     const SurfaceMesh& surf_in_A,
        //     const SurfaceMesh& surf_in_B,
        //     const RepairParameters& args,
        //     SurfaceMesh& surf_out,
        //     std::vector<DiscreteGeometry::Edge>& edges_out  // 抽取交线
        // );

        /**
         * @brief 投影
         * 
         * @param[in] surf_in 输入网格
         * @param[in] edges_in  指定边
         * @param[in] tris_id 指定投影面三角形
         * @param[in] args  指定print_eps
         * @param[out] surf_out 投影之后重新重网格化
         */
        Tiger_API void project(
            const SurfaceMesh& surf_in,
            const std::vector<DiscreteGeometry::Edge>& edges_in,
            const std::vector<int>& tris_id,
            const RepairParameters& args,
            SurfaceMesh& surf_out);

        /**
         * @brief 压印
         * 
         * @param[in] edges_in
         * @param[in] surf_in
         * @param[out] surf_out
         */
        // Tiger_API void print(
        //     const std::vector<DiscreteGeometry::Edge>& edges_in,
        //     const SurfaceMesh& surf_in,
        //     SurfaceMesh& surf_out);
    }
}