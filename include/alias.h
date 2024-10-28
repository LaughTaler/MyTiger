#ifndef _ALIAS_DEFINE_H_
#define _ALIAS_DEFINE_H_

#define TMeshSurf_processTriangleLoop triangleLoop
#define TMeshSurf_Repair repairDiscreteModel
#define CADSurf_exportToSTL discretizeSurfaceTri_STL

#define Grid1D_configureLineMesh setLineMesh
#define Grid1D_Legacy discreteCurve
#define Grid1D_Dimension discreteCurveDimension
#define Grid1D_Tanh discreteCurvesTanh
#define Grid1D_Geometric discreteCurvesGeometric
#define Grid1D_Shape discreteCurvesShape
#define Grid1D_Automate discreteCurvesAutomatic

#define HexGrid_createPrismaticMesh generatePrismaticMesh
#define HexGrid_createHybridMesh generateHybridMesh
#define HexGrid_Legacy generateHexMesh

#define TMeshSurf_Opt triangleRemesh
#define CADSurf_TGrid_AFT discretizeSurfaceTri
#define CADSurf_TGrid_AFLR discretizeSurfaceTriAFLR
#define CADSurf_QGrid_AFLR discretizeSurfaceQuadAFLR
#define VCADSurf_TGrid_AFT discretizeSurfacesTri
#define VSurf_TMesh_meshTriangulately discretizeVirtualSurfaceTri
#define CADSurf_QGrid_Paving discretizeSurfaceQuad
#define VSurf_QMesh_Paving discretizeVirtualSurfaceQuad
#define CADSurf_QGrid_Submap discretizeSurfaceStrucQuad
#define VSurf_QGrid_Submap discretizeVirtualSurfaceStrucQuad

#define Tool_initializeData init_data
#define Tool_convertToDiscreteGeometry convertTODiscreteGeometry

#define TMeshSurf_detectSurfacePunctures detectPunctureSurface
#define TMeshSurf_assessSurfaceQuality detectSurfaceQuality
#define TMeshSurf_identifyMeshAdjacency detectAdjacence
#define TMeshSurf_identifyFreeEdges detectFreeEdge
#define TMeshSurf_identifyMultipleEdges detectMultiedge
#define TMeshSurf_identifyMultiplePoints detectMultipoint
#define TMeshSurf_removeTriangles deletedTris
#define TMeshSurf_createTrianglesFromPoints createTriByPoints
#define TMeshSurf_identifyOverlappingPoints overlapPoints
#define TMeshSurf_splitMeshEdges splitEdges
#define TMeshSurf_splitMeshTriangles splitTris
#define TMeshSurf_flipMeshEdges flipEdges
#define TMeshSurf_fillPolygonWithTriangles fillPolygonByTris
#define TMeshSurf_automateSurfaceRepair autoRepairSurface
#define TMeshSurf_createTrianglesFromEdges CreateTrisByEdges
#define TMeshSurf_sutureEdges SutureEdges
#define TMeshSurf_smoothMeshPoints smoothPoints
#define TMeshSurf_mergeAdjacentMeshPoints mergeAdjacentPoints
#define TMeshSurf_mergeDuplicateSurfaces mergeDuplicateSurfaces

#define HexGrid_Cartesian generateCartesianMesh
#define HexGrid_ParallelRefine distributeParalleRefine
#define PolyGrid_Dual generatePolygonalMesh
#define HexGrid_Sweep generateSweepHexMesh
#define HexGrid_Map generateMappingHexMesh
#define HexGrid_Extrude2D extrude2d
#define HexGrid_Extrude3D extrude3d

// #define Read_Model readModel_OCC
// #define CAD_Tessellation generateStl_OCC
#define VCAD_Tessellation generateStl_Binary
#define Tool_SurfQuality calculateSurfaceMeshQuality
#define Tool_VolQuality calculateVolMeshQuality

#define TetGrid_Constrained constrainDelaunay
#define TetGrid_BoxConforming boxConformingDelaunay
#define TetGrid_Conforming conformingDelaunay
#define TetGrid_opt tetOptimization

#define SizingFunction_setSizeFromVolume setSizeFieldwithVolumeMesh
#define SizingFunction_setUniformSize setConstSize
#define SizingFunction_addPointSource addPointSource
#define SizingFunction_addLineSource addLineSource
#define SizingFunction_addTriangleSource addTriangleSource
#define SizingFunction_addCubicSource addCubicSource
#define SizingFunction_getSizeAtPoint getSizeValue
#define SizingFunction_configureGeometryField setGeometrySizeField
#define SizingFunction_addPolygonSource addPolygonSource
#define SizingFunction_addEllipsoidSource addEllipsoidSource
#define SizingFunction_addFrustumSource addFrustumSource
#define SizingFunction_defineCustomSize addUserDefineSize

#define Tool_writeMeshToCGNS writeCGNS
#define Grid1D_createMeshLine creatLine
#define Grid2D_createMeshPlane creatPlane
#define Tool_detectPointDuplicates  detectDuplicatePoints
#define Tool_nativelyDetectPointDuplicates native_detectDuplicatePoints
#define Tool_projectPointToSegment project_PointOntoSegment
#define Tool_projectPointToTriangle project_PointOntoTriangle
#define Tool_intersectLines line_line_intersection
#define Tool_intersectLines2D line_line_2dintersection
#define Tool_intersectLineWithTriangle line_tri_intersection
#define Tool_intersectLineWithTriangle2D line_tri_2dintersection
#define Tool_intersectPointWithTriangle point_tri_intersection

#endif //!_ALIAS_DEFINE_H_
