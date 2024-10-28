#pragma once
#ifndef CREATE_GEOMETRY_H_
#define CREATE_GEOMETRY_H_
#include "geometry_data_struct.h"
#include "mesh_data_struct.h"
#include "./definitions.h"
#include "alias.h"
#include "./geometry_interface_parameters.h"

namespace TiGER
{
/**
 * @brief 读取几何到内存中
 *
 * @param[in] filename 读入文件名，支持 stp, igs, step, iges, sat, brep, x_t ...
 **/
Tiger_API void Read_Model(const std::string& filename, TiGER::Geometry& geometry, const readModelParameters& args);


// Tiger_API void readModel_Ferguson(const std::string& filename, TiGER::Geometry& geometry,
//            const readModelParameters &args);

/**
 * @brief 生成STL（离散背景网格）
 *
 * @param[in] geometry 几何文件
 * @param[out] surfaceOut 以面网格的形式生成的STL
 **/

Tiger_API void CAD_Tessellation(TiGER::Geometry& geometry, TiGER::SurfaceMesh& surfaceOut,
                               TiGER::generateStlParameters& arg_stl);

Tiger_API void VCAD_Tessellation(TiGER::Geometry& geometry, TiGER::SurfaceMesh& surfaceOut,
                                  TiGER::generateStlDTParameters& arg_stldt);

Tiger_API void mapReal2Normal_Curve(std::shared_ptr<GeometryCurve> curve,double& real,double& normal);

Tiger_API void mapNormal2Real_Curve(std::shared_ptr<GeometryCurve> curve,double& real,double& normal);

Tiger_API std::shared_ptr<GeometryCurve> trimmedCurve(std::shared_ptr<GeometryCurve> curve,double& real1,double& real2);

Tiger_API void dealcomplexcurve(TiGER::Geometry& Geo);

} // namespace TiGER
#endif