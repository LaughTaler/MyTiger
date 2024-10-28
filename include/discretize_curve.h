#include "geometry_data_struct.h"
#include "./sizing_field_interface.h"
#include "./mesh_data_struct.h"
#include "./definitions.h"
#include "./discretize_curve_parameters.h"
#include "alias.h"
#include <memory>
#include <string>
/**
 * @brief 曲线离散函数集 
 *  
 */
namespace TiGER{
namespace discretize_curve{


    Tiger_API void Grid1D_configureLineMesh(std::shared_ptr<TiGER::GeometryCurve>curve,const std::vector<std::array<double,3>>& points);

    /**
     * @brief 根据尺寸场离散虚曲线 
     * 
     * @param curves 一组收尾相连的几何曲线 
     * @param sizefield 尺寸场 
     * @param discrete_curve 离散曲线结果 
     */

    Tiger_API void Grid1D_Legacy(
        const std::shared_ptr<TiGER::GeometryCurve>curve,
        const SizingFunction& sizefield, 
        const TiGER::CurveParameters& args,
        TiGER::LineMesh& segments);

    /**
     * @brief 根据尺寸场离散曲线的部分曲线段
     *
     * @param curve 一个收尾相连的几何曲线
     * @param StartPt 曲线段首端点
     * @param EndPt 曲线段末端点
     * @param sizefield 尺寸场
     * @param discrete_curve 离散曲线结果
     */

    Tiger_API void discretePartCurve(const std::shared_ptr<TiGER::GeometryCurve> curve,
                                     const std::array<double, 3> StartPt, const std::array<double, 3> EndPt,
                                     const SizingFunction& sizefield, const TiGER::CurveParameters& args,
                                     TiGER::LineMesh& segments);


    /**
     * @brief 以Dimension方式离散一组曲线
     * 类型包括 
     * 
     * @param curves 一组曲线
     * @param args Dimension参数
     * @param segements 离散结果
     */
    Tiger_API void Grid1D_Dimension(
        const std::vector<std::shared_ptr<TiGER::GeometryCurve>>& curves, 
        const std::vector<std::vector<std::shared_ptr<TiGER::GeometrySurface>>>&surfaces,
        const TiGER::CurveParametersDimension& args,
        std::vector<TiGER::LineMesh>& segements
    );
    /**
     * @brief 以Shape等特殊方式离散一组虚曲线
     * 类型包括 
     * 
     * @param curves 一组虚曲线
     * @param sizefield 尺寸场
     * @param discrete_curve 离散结果
     */
    Tiger_API void Grid1D_Tanh(
        const std::vector<std::shared_ptr<GeometryCurve>>& curves, 
        const TiGER::CurveParametersTanh& args, 
        std::vector<TiGER::LineMesh>& segements);

    /**
     * @brief 以Shape等特殊方式离散一组虚曲线
     * 类型包括 
     * 
     * @param curves 一组虚曲线
     * @param sizefield 尺寸场
     * @param discrete_curve 离散结果
     */
    Tiger_API void Grid1D_Geometric(
        const std::vector<std::shared_ptr<GeometryCurve>>& curves, 
        const TiGER::CurveParametersGeometric& args, 
        std::vector<TiGER::LineMesh>& segements);
    /**
     * @brief 以Growth等特殊方式离散一组虚曲线
     * 类型包括 
     * 
     * @param curves 一组虚曲线
     * @param sizefield 尺寸场
     * @param discrete_curve 离散结果
     */
    Tiger_API void Grid1D_Shape(
        const std::vector<TiGER::GeometryCurve>& curves, 
        const TiGER::CurveParametersGrowth& args, 
        std::vector<TiGER::LineMesh>& segements);


    /**
     * @brief 以Automatic等特殊方式离散一组虚曲线
     * 类型包括 
     * 
     * @param curves 一组虚曲线
     * @param sizefield 尺寸场
     * @param discrete_curve 离散结果
     */
    Tiger_API void Grid1D_Automate(
        const std::vector<TiGER::GeometryCurve>& curves, 
        const TiGER::CurveParametersAutomatic& args, 
        std::vector<TiGER::LineMesh>& segements);


} // discretize_curve
} // TiGER