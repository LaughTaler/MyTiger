#pragma once
#ifndef _TIGER_EXCEPTION_H_
#define _TIGER_EXCEPTION_H_
#include <exception>
#include <string>
namespace TiGER{


class CustomException : public std::exception {
    std::string message;
public:
    CustomException(const std::string& msg) : message(msg) {}
    virtual const char* what() const noexcept override {
        return message.c_str();
    }
};

#define DEFINE_SPECIFIC_EXCEPTION(NAME, MESSAGE_PREFIX,MESSAGE_HINT,MESSAGE_EXTRA) \
class NAME : public CustomException { \
public: \
    NAME(const std::string& message) : CustomException(MESSAGE_PREFIX+'(' MESSAGE_HINT + ")--"+ MESSAGE_EXTRA + " :"+message) {} \
};

DEFINE_SPECIFIC_EXCEPTION(BoundaryLayerMesh_InputMeshHasOpenEdge, "The input mesh is not closed", "The algorithm only supports closed 2-manifold triangular meshes", "IDs of the two points of the first open edge");
DEFINE_SPECIFIC_EXCEPTION(BoundaryLayerMesh_InputMeshIsNonTwoManifold, "The input mesh has an edge shared by more than two triangles", "The algorithm only supports closed 2-manifold triangular meshes", "IDs of the two points of the first non-manifold edge");
DEFINE_SPECIFIC_EXCEPTION(BoundaryLayerMesh_InputMeshHasSelfIntersections, "The input mesh has self-intersections", "The algorithm only supports non-self-intersecting meshes", "ID of the first intersecting triangle");
DEFINE_SPECIFIC_EXCEPTION(BoundaryLayerMesh_InsufficientMemory, "Insufficient memory for mesh generation", "Consider reducing the size of the surface mesh", "");
DEFINE_SPECIFIC_EXCEPTION(BoundaryLayerMesh_InvalidInputParameters, "Invalid input parameters", "Please refer to the documentation for correct parameters", "");

DEFINE_SPECIFIC_EXCEPTION(DelaunayTriangulation_PointCollinearOnBoundaryEdge, "The input mesh has points collinear with boundary edges", "The algorithm does not support collinear points and edges", "Collinearity occurred, {} located in {} {}");
DEFINE_SPECIFIC_EXCEPTION(DelaunayTriangulation_InternalError, "Neighbor error due to topological transformation error", "Internal error", "");
DEFINE_SPECIFIC_EXCEPTION(DelaunayTriangulation_InputMeshHasSelfIntersections, "The input mesh has self-intersections", "The algorithm only supports non-self-intersecting meshes", "ID of the first intersecting triangle");

DEFINE_SPECIFIC_EXCEPTION(AutoGrid_InputMeshHasPoorQuality, "The input triangular mesh has poor quality", "Consider improving the quality of the triangular mesh", "Energy value of the worst triangle");
DEFINE_SPECIFIC_EXCEPTION(AutoGrid_InsufficientMemory, "Insufficient memory for mesh generation", "Consider reducing the size of the surface mesh", "");
DEFINE_SPECIFIC_EXCEPTION(AutoGrid_InvalidInputParameters, "Invalid input parameters", "Please refer to the documentation for correct parameters", "");

DEFINE_SPECIFIC_EXCEPTION(Remesh_InputMeshHasSelfIntersections, "The input mesh has self-intersections", "The algorithm only supports non-self-intersecting meshes", "Number of intersecting triangles");
DEFINE_SPECIFIC_EXCEPTION(Remesh_InvalidInputParameters, "Invalid input parameters", "Please refer to the documentation for correct parameters", "");

DEFINE_SPECIFIC_EXCEPTION(SMesh_GeometryHasSelfIntersections, "The input geometry has self-intersections", "Only non-self-intersecting geometries are supported", "");
DEFINE_SPECIFIC_EXCEPTION(SMesh_GeometryHasOpenLoops, "The input geometry has open loops", "Only closed-loop geometries are supported", "");
DEFINE_SPECIFIC_EXCEPTION(SMesh_InputSizeFieldHasZeroValues, "Invalid size field values", "Only size fields with values greater than or equal to zero are supported", "");



}

#endif //!_TIGER_EXCEPTION_H_