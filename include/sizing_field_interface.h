#pragma once
#include <vector>
#include <array>
#include <memory>

#include "./geometry_data_struct.h"
#include "./definitions.h"
#include "./binary_tree.hpp"
#include "sizing_field_parameters.h"
// #include <igl/AABB.h>
#include <Eigen/Dense>

namespace TiGER{
  typedef std::function<double(const double&,const double&,const double&)> SizingFunction;


  class sizeobj{
    public:
    // SizingFunction getsize;
    virtual double getsize(const double& x, const double& y, const double& z) = 0;
  };

  class constsizefield: public sizeobj{
    public:
    double size;

    double getsize(const double& x, const double& y, const double& z);
  };

  class sourse: public sizeobj{
    public:
    int spacingtype;
  };

  class generalsource: public sourse{
    public:
    SourceParameters args;
    SurfaceMesh backgroundmesh;

    double getsize(const double& x, const double& y, const double& z);
  };

  class pointSource: public sourse{
    public:
    PointSourceParameters args;
    std::array<double,3> xyz;

    double getsize(const double& x, const double& y, const double& z);
  };

  class lineSource: public sourse{
    public:
    LineSourceParameters args;
    std::array<std::array<double,3>,2> xyz;

    double getsize(const double& x, const double& y, const double& z);
  };

  
  class triangleSource: public sourse{
    public:
    TriangleSourceParameters args;
    std::array<std::array<double,3>,3> xyz;

    double getsize(const double& x, const double& y, const double& z);
  };



  /*
             0---------1
            /|        /|
           / |       / |
          /  |      /  |
         3---------2   |
         |   |     |   |
         |   4-----|---5
         |  /      |  /
         | /       | /
         7---------6
*/

  class cubicSource: public sourse{
    public:
    CubicSourceParameters args;
    std::array<std::array<double,3>,8> xyz;

    double getsize(const double& x, const double& y, const double& z);
  };

  class cylinderSource: public sourse{
    public:
    CylinderSourceParameters args;
    std::array<std::array<double,3>,2> xyz;
    std::array<double,2> radius;
    double getsize(const double& x, const double& y, const double& z);
  };

    class sphereSource: public sourse{
    public:
    SphereSourceParameters args;
    std::array<double,3> xyz;
    std::array<double,1> radius;
    double getsize(const double& x, const double& y, const double& z);
  };





  class backgroundsizefield: public sizeobj{
    public:
    backgroundsizefield(){}
    BackGroundParameters args;
    SurfaceMesh backgroundmesh;
    
    // common::BinaryTree<Eigen::Matrix3d> octree_;
    TiGER::common::BinaryTree<Eigen::Matrix3d> octree_;

    
    double getsize(const double& x, const double& y, const double& z);
    // void exportSizeFieldToVTK(const std::string &filePath);

};


  struct SizingManager{
      public:
      int combinetype = 0;
      std::vector<std::shared_ptr<sizeobj>> sf;
  };



} // TiGER


namespace TiGER{
  namespace size_field{
   Tiger_API void setSizeField(const SurfaceMesh& SurfaceMesh,const BackGroundParameters& args,SizingManager& manager);

    // Tiger_API void SizingFunction_setSizeFromVolume(const VolumeMesh& VolumeMesh,SizingManager& manager,double scale = 1);

    Tiger_API void SizingFunction_setUniformSize(double size,SizingManager& manager);

    Tiger_API void SizingFunction_addPointSource(const std::array<double,3>& xyz,const PointSourceParameters& args,SizingManager& manager);

    Tiger_API void SizingFunction_addLineSource(const std::array<std::array<double,3>,2>& xyz,const LineSourceParameters& args,SizingManager& manager);

    Tiger_API void SizingFunction_addTriangleSource(const std::array<std::array<double,3>,3>& xyz,const TriangleSourceParameters& args,SizingManager& manager);

    Tiger_API void SizingFunction_addCubicSource(const std::array<std::array<double,3>,8>& points,const CubicSourceParameters& args,SizingManager& manager);

    Tiger_API void addSphereSource(const std::array<double,3>& xyz, const std::array<double,1>& radius, SphereSourceParameters& args, SizingManager& manager);

    Tiger_API void addCylinderSource(const std::array<std::array<double,3>,2>& xyz, const std::array<double,2>& radius, CylinderSourceParameters& args, SizingManager& manager);

    Tiger_API void setAttributeSize(const std::string& attribute_name,const std::vector<int>& attribute_index,const double& sizing_value,backgroundsizefield& sf);

    Tiger_API double SizingFunction_getSizeAtPoint(const double& x, const double& y, const double& z, SizingManager& manager);

    Tiger_API void API_Enable_Source(int obj);

    Tiger_API void API_Disable_Source(int obj);

    Tiger_API void API_Delete_Source(int obj);

    
  }



    void setGeometrySizeField(
      const Geometry& geometry,
      const SourceParameters& args,
      SizingManager& sf
      ) ;
    




    void addPolygonSource(const std::vector<std::array<double,3>>& xyz,
        const std::vector<double>& inner_radius,
        const std::vector<double>& outter_radius,
        const double& sizing_value,
        SizingManager& manager);
        

    void addEllipsoidSource(const std::array<double,3>& xyz1,
    const std::array<double,3>& xyz2,
        const double& inner_radius,
        const double& outter_radius,
        const double& sizing_value,
        SizingManager& manager);



    void addFrustumSource(const std::array<double,3>& upper_xyz,
        const std::array<double,3>& lower_xyz,
        const double& upper_radius,
        const double& lower_radius,
        const double& outter_radius,
        const double& sizing_value,
        SizingManager& manager);

    void addCubicSource(
        const std::array<std::array<double,3>,8>& points,
        const double& outter_radius,
        const double& sizing_value,
        SizingManager& manager);
        

    void addUserDefineSize(
        const std::function<std::array<double,9>(std::array<double,3>)> function,
        SizingManager& manager);
 

    
}