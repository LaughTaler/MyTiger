#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <Eigen/Dense>
#include<memory>

#include "../../include/geometry_data_struct.h"
#include "../../include/geometry_interface_occ.h"
#include "./geometry_occ.h"
#include <GCPnts_AbscissaPoint.hxx>

using namespace std;

namespace TiGER
{
struct curve_hash_fn
{
    size_t operator()(const TopoDS_Edge& edge) const
    {
        gp_Pnt p1, p2, p3;
        BRepAdaptor_Curve curve(edge);
        curve.D0(curve.FirstParameter(), p1);
        curve.D0(curve.LastParameter(), p2);
        curve.D0((curve.FirstParameter() + curve.LastParameter()) / 2, p3);
        std::size_t h1 = std::hash<double>()(std::round(p1.X() * 10000) / 10000);
        std::size_t h2 = std::hash<double>()(std::round(p1.Y() * 10000) / 10000);
        std::size_t h3 = std::hash<double>()(std::round(p1.Z() * 10000) / 10000);
        std::size_t h4 = std::hash<double>()(std::round(p2.X() * 10000) / 10000);
        std::size_t h5 = std::hash<double>()(std::round(p2.Y() * 10000) / 10000);
        std::size_t h6 = std::hash<double>()(std::round(p2.Z() * 10000) / 10000);
        std::size_t h7 = std::hash<double>()(std::round(p3.X() * 10000) / 10000);
        std::size_t h8 = std::hash<double>()(std::round(p3.Y() * 10000) / 10000);
        std::size_t h9 = std::hash<double>()(std::round(p3.Z() * 10000) / 10000);
        return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4) ^ (h5 << 1) ^ (h6 << 2) ^ (h7 << 3) ^ (h8 << 4) ^ (h9 << 5);
    }
};

struct curve_equal_fn
{
    bool operator()(const TopoDS_Edge& e0, const TopoDS_Edge& e1) const
    {
        return e0.IsSame(e1);
    }
};
TopoDS_Solid makeSolidFromFaces(const std::vector<TopoDS_Face>& faces)
{
    // 创建一个空的Shell
    TopoDS_Shell shell;
    BRep_Builder builder;
    builder.MakeShell(shell);

    // 将每个Face添加到Shell中
    for( const TopoDS_Face& face : faces )
    {
        builder.Add(shell, face);
    }

    // 检查是否成功生成Shell
    if( shell.IsNull() )
    {
        std::cerr << "Shell creation failed!" << std::endl;
        return TopoDS_Solid(); // 返回一个空的Solid
    }

    // 使用Shell创建Solid
    BRepBuilderAPI_MakeSolid solidMaker(shell);

    // 检查是否成功生成Solid
    if( !solidMaker.IsDone() )
    {
        std::cerr << "Solid creation failed!" << std::endl;
        return TopoDS_Solid(); // 返回一个空的Solid
    }

    // 返回生成的Solid
    return solidMaker.Solid();
}
void ClearGeometry(Geometry& geometry)
{
    // 清空几何曲线和几何面的智能指针容器，自动释放内存
    geometry.curves_.clear();
    geometry.surfaces_.clear();

    // 清空拓扑信息
    geometry.topo_surf_to_curves_.clear();
    geometry.topo_body_to_surfs_.clear();

    // 清空面名字信息
    geometry.mp_index_name.clear();
}

int performBoolean(const TopTools_ListOfShape& shapeList, TopoDS_Shape& resultShape)
{
    if( shapeList.IsEmpty() )
    {
        return -1;
    }
    BOPAlgo_MakerVolume aMV;
    aMV.SetArguments(shapeList);
    // aMV.SetParallelMode(true);
    aMV.SetRunParallel(true);
    // Additional option of the algorithm
    Standard_Boolean bAvoidInternalShapes = Standard_True;
    aMV.SetAvoidInternalShapes(bAvoidInternalShapes);

    // Perform the operation
    aMV.Perform();
    if( aMV.HasErrors() )
    { // check error status
        return -1;
    }
    resultShape = aMV.Shape();

    if( resultShape.NbChildren() == 0 )
    {
        resultShape = shapeList.First();
        return -1;
    }
    return 0;
}

int cutEnclosedFaces(TopoDS_Shape& shape, int splitPoints)
{
    ShapeUpgrade_ShapeDivideClosed shapeDivideTool(shape);
    shapeDivideTool.SetNbSplitPoints(splitPoints);
    Standard_Boolean res = shapeDivideTool.Perform();
    if( res )
    {
        shape = shapeDivideTool.Result();
    }
    return 0;
}

int fixDegeneratedGeo(TopoDS_Shape& shape)
{
    TopExp_Explorer edgeExplorer;
    ShapeBuild_ReShape rebuild;
    for( edgeExplorer.Init(shape, TopAbs_EDGE); edgeExplorer.More(); edgeExplorer.Next() )
    {
        TopoDS_Edge edge = TopoDS::Edge(edgeExplorer.Current());
        if( BRep_Tool::Degenerated(edge) )
            rebuild.Remove(edge);
    }
    shape = rebuild.Apply(shape);
    TopExp_Explorer faceExplorer;
    for( faceExplorer.Init(shape, TopAbs_FACE); faceExplorer.More(); faceExplorer.Next() )
    {
        TopoDS_Face face = TopoDS::Face(faceExplorer.Current());

        ShapeFix_Face sff(face);
        sff.FixAddNaturalBoundMode() = Standard_True;
        sff.FixSmallAreaWireMode() = Standard_True;
        sff.Perform();

        if( sff.Status(ShapeExtend_DONE1) || sff.Status(ShapeExtend_DONE2) || sff.Status(ShapeExtend_DONE3) ||
            sff.Status(ShapeExtend_DONE4) || sff.Status(ShapeExtend_DONE5) )
        {
            TopoDS_Face newface = sff.Face();
            rebuild.Replace(face, newface);
        }
    }
    shape = rebuild.Apply(shape);

    TopExp_Explorer edgeExplorer1;
    for( edgeExplorer1.Init(shape, TopAbs_EDGE); edgeExplorer1.More(); edgeExplorer1.Next() )
    {
        TopoDS_Edge edge = TopoDS::Edge(edgeExplorer1.Current());
        if( BRep_Tool::Degenerated(edge) )
            rebuild.Remove(edge);
    }
    shape = rebuild.Apply(shape);

    return 0;
}

int fixSmallEdges(TopoDS_Shape& shape, double tolerance)
{
    std::cout << ("Shape Healing: Fixing small edges") << std::endl;

    TopExp_Explorer faceExp;
    ShapeBuild_ReShape rebuild;

    for( faceExp.Init(shape, TopAbs_FACE); faceExp.More(); faceExp.Next() )
    {
        TopoDS_Face face = TopoDS::Face(faceExp.Current());
        TopExp_Explorer wireExp;
        for( wireExp.Init(face, TopAbs_WIRE); wireExp.More(); wireExp.Next() )
        {
            TopoDS_Wire oldwire = TopoDS::Wire(wireExp.Current());
            ShapeFix_Wire sfw(oldwire, face, tolerance);
            sfw.ModifyTopologyMode() = Standard_True;
            sfw.ClosedWireMode() = Standard_True;
            bool replace = false;
            replace = sfw.FixReorder() || replace;
            replace = sfw.FixConnected() || replace;

            if( sfw.FixSmall(Standard_False, tolerance) &&
                !(sfw.StatusSmall(ShapeExtend_FAIL1) || sfw.StatusSmall(ShapeExtend_FAIL2) ||
                  sfw.StatusSmall(ShapeExtend_FAIL3)) )
            {
                std::cout << (" . Fixed small edge in wire");
                replace = true;
            }
            else if( sfw.StatusSmall(ShapeExtend_FAIL1) )
                std::cout << "Failed to fix small edge in wire, edge cannot be checked (no 3d curve and no pcurve)"
                          << std::endl;
            else if( sfw.StatusSmall(ShapeExtend_FAIL2) )
                std::cout << "Failed to fix small edge in wire, edge is null-"
                          << "length and has different vertives at begin and end, "
                          << "and lockvtx is True or ModifiyTopologyMode is False" << std::endl;
            else if( sfw.StatusSmall(ShapeExtend_FAIL3) )
                std::cout << "Failed to fix small edge in wire, CheckConnected has failed" << std::endl;

            replace = sfw.FixEdgeCurves() || replace;
            replace = sfw.FixDegenerated() || replace;
            replace = sfw.FixSelfIntersection() || replace;
            replace = sfw.FixLacking(Standard_True) || replace;
            if( replace )
            {
                TopoDS_Wire newwire = sfw.Wire();
                rebuild.Replace(oldwire, newwire);
            }
        }
    }
    shape = rebuild.Apply(shape);

    TopExp_Explorer edgeExp;
    for( edgeExp.Init(shape, TopAbs_EDGE); edgeExp.More(); edgeExp.Next() )
    {
        TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());
        GProp_GProps system;
        BRepGProp::LinearProperties(edge, system);
        if( system.Mass() < tolerance )
        {
            std::cout << (" Removing degenerated edge") << std::endl;
            rebuild.Remove(edge);
        }
    }
    shape = rebuild.Apply(shape);

    TopExp_Explorer edgeExp1;
    for( edgeExp1.Init(shape, TopAbs_EDGE); edgeExp1.More(); edgeExp1.Next() )
    {
        TopoDS_Edge edge = TopoDS::Edge(edgeExp1.Current());
        if( BRep_Tool::Degenerated(edge) )
            rebuild.Remove(edge);
    }
    shape = rebuild.Apply(shape);

    ShapeFix_Wireframe sfwf;
    sfwf.SetPrecision(tolerance);
    sfwf.Load(shape);
    sfwf.ModeDropSmallEdges() = Standard_True;

    if( sfwf.FixWireGaps() )
    {
        std::cout << (" - Fixing wire gaps") << std::endl;
        if( sfwf.StatusWireGaps(ShapeExtend_OK) )
            std::cout << ("  no gaps found") << std::endl;
        if( sfwf.StatusWireGaps(ShapeExtend_DONE1) )
            std::cout << (" . Some 2D gaps fixed") << std::endl;
        if( sfwf.StatusWireGaps(ShapeExtend_DONE2) )
            std::cout << (" . Some 3D gaps fixed") << std::endl;
        if( sfwf.StatusWireGaps(ShapeExtend_FAIL1) )
            std::cout << (" . Failed to fix some 2D gaps") << std::endl;
        if( sfwf.StatusWireGaps(ShapeExtend_FAIL2) )
            std::cout << (" . Failed to fix some 3D gaps") << std::endl;
    }

    sfwf.SetPrecision(tolerance);

    if( sfwf.FixSmallEdges() )
    {
        std::cout << ("Fixing wire frames") << std::endl;
        if( sfwf.StatusSmallEdges(ShapeExtend_OK) )
            std::cout << (" . No small edges found") << std::endl;
        if( sfwf.StatusSmallEdges(ShapeExtend_DONE1) )
            std::cout << (" . Some small edges fixed") << std::endl;
        if( sfwf.StatusSmallEdges(ShapeExtend_FAIL1) )
            std::cout << (" . Failed to fix some small edges") << std::endl;
    }

    shape = sfwf.Shape();
    return 0;
}

int fixSmallFaces(TopoDS_Shape& shape, double tolerance)
{
    std::cout << ("Shape Healing: Fixing spot and strip edges") << std::endl;
    ShapeFix_FixSmallFace sffsm;
    sffsm.Init(shape);
    sffsm.SetPrecision(tolerance);
    sffsm.Perform();
    shape = sffsm.FixShape();
    return 0;
}

int sewFaces(TopoDS_Shape& shape, double tolerance)
{
    BRepOffsetAPI_Sewing sewedObj(tolerance);
    TopExp_Explorer faceExp;
    for( faceExp.Init(shape, TopAbs_FACE); faceExp.More(); faceExp.Next() )
    {
        TopoDS_Face face = TopoDS::Face(faceExp.Current());
        sewedObj.Add(face);
    }

    sewedObj.Perform();

    if( !sewedObj.SewedShape().IsNull() )
    {
        shape = sewedObj.SewedShape();
    }
    else
    {
        return -1;
    }
    return 0;
}

int makeSolids(TopoDS_Shape& shape, double tolerance)
{
    std::cout << ("Shape Healing: Making solids") << std::endl;

    BRepBuilderAPI_MakeSolid ms;
    int count = 0;
    TopExp_Explorer shellExp;
    for( shellExp.Init(shape, TopAbs_SHELL); shellExp.More(); shellExp.Next() )
    {
        count++;
        ms.Add(TopoDS::Shell(shellExp.Current()));
    }

    if( count == 0 )
    {
        return -1;
    }
    else
    {
        BRepCheck_Analyzer ba(ms);
        if( ba.IsValid() )
        {
            ShapeFix_Shape sfs;
            sfs.Init(ms);
            sfs.SetPrecision(tolerance);
            sfs.SetMaxTolerance(tolerance);
            sfs.Perform();
            shape = sfs.Shape();
            TopExp_Explorer solidExp;
            for( solidExp.Init(shape, TopAbs_SOLID); solidExp.More(); solidExp.Next() )
            {
                TopoDS_Solid solid = TopoDS::Solid(solidExp.Current());
                TopoDS_Solid newsolid = solid;
                BRepLib::OrientClosedSolid(newsolid);
                ShapeBuild_ReShape rebuild;
                rebuild.Replace(solid, newsolid);
                shape = rebuild.Apply(shape, TopAbs_COMPSOLID);
            }
        }
        else
        {
            return -1;
        }
    }
}

int travelShape(vector<TopoDS_Shape> shapes, TiGER::Geometry& geometry)
{
    TopExp_Explorer exp_model, exp_quilt, exp_wire, exp_curve;
    std::unordered_map<TopoDS_Edge, int, curve_hash_fn, curve_equal_fn> curve_id_map;
    std::unordered_map<TopoDS_Face, std::string, quilt_hash_fn, quilt_equal_fn> mp;
    for( auto shape : shapes )
    {
        std::vector<int> temp_surf_on_body;
        for( exp_quilt.Init(shape, TopAbs_FACE); exp_quilt.More(); exp_quilt.Next() )
        {
            TopoDS_Face topo_face = TopoDS::Face(exp_quilt.Current());
            BRepAdaptor_Surface surface(topo_face);
            std::shared_ptr<GeometrySurface_OCC> surface_temp =
                std::make_shared<GeometrySurface_OCC>(topo_face, surface);
            int surfaceid = geometry.surfaces_.size();

            int idd = geometry.surfaces_.size();
            geometry.surfaces_.push_back(surface_temp);
            if( mp.count(surface_temp->Ref_face) != 0 )
            {
                geometry.mp_index_name[idd] = mp[surface_temp->Ref_face];
            }

            temp_surf_on_body.push_back(surfaceid);
            std::vector<int> temp_curve_on_face;

            for( exp_curve.Init(topo_face, TopAbs_EDGE); exp_curve.More(); exp_curve.Next() )
            {
                TopoDS_Edge topo_edge = TopoDS::Edge(exp_curve.Current());
                int curveid;
                if( curve_id_map.count(topo_edge) )
                {
                    curveid = curve_id_map[topo_edge];
                }
                else
                {
                    BRepAdaptor_Curve curve(topo_edge);
                    double first = curve.FirstParameter();
                    double last = curve.LastParameter();
                    double length = GCPnts_AbscissaPoint::Length(curve, first, last);
                    if (length < 1e-6) {
                        continue;
                    }
                    std::shared_ptr<GeometryCurve_OCC> curve_temp =
                        std::make_shared<GeometryCurve_OCC>(topo_edge, curve);
                    curveid = geometry.curves_.size();
                    curve_id_map[topo_edge] = curveid;
                    geometry.curves_.push_back(curve_temp);
                }

                temp_curve_on_face.push_back(curveid);
            }
            geometry.topo_surf_to_curves_.push_back(temp_curve_on_face);
        }
        geometry.topo_body_to_surfs_.push_back(temp_surf_on_body);
    }
    return 0;
}

int repair(TiGER::Geometry& geometry, int heal)
{
    // heal fix model
    bool fixDegenerated_fg, fixSmallEdges_fg, fixSmallFaces_fg, sewFaces_fg, makeSolids_fg, mydeletesmalledge_fg,
        PreserveFeatureFaces_fg;
    double tolerance = 1e-8;
    fixDegenerated_fg = (heal & 0b00001);
    fixSmallEdges_fg = (heal & 0b00010) >> 1;
    fixSmallFaces_fg = (heal & 0b00100) >> 2;
    sewFaces_fg = (heal & 0b01000) >> 3;
    makeSolids_fg = (heal & 0b10000) >> 4;
    mydeletesmalledge_fg = (heal & 0b100000) >> 5;
    PreserveFeatureFaces_fg = (heal & 0b1000000) >> 6;

    TopTools_ListOfShape solids;
    for( auto model_topo : geometry.topo_body_to_surfs_ )
    {
        std::vector<TopoDS_Face> faces;
        for( auto face_topo : model_topo )
        {
            std::shared_ptr<GeometrySurface> face_geo = geometry.surfaces_[face_topo];
            GeometrySurface_OCC* surface = static_cast<GeometrySurface_OCC*>(face_geo->getptr());

            faces.push_back(surface->Ref_face);
            TopoDS_Solid solid = makeSolidFromFaces(faces);
            if( solid.IsNull() )
            {
                return -1;
            }
            solids.Append(solid);
        }
    }

    if( fixDegenerated_fg )
    {
        for( auto it = TopTools_ListIteratorOfListOfShape(solids); it.More(); it.Next() )
        {
            auto& shape = it.Value();
            fixDegeneratedGeo(shape);
        }
    }

    if( fixSmallEdges_fg )
    {
        for( auto it = TopTools_ListIteratorOfListOfShape(solids); it.More(); it.Next() )
        {
            auto& shape = it.Value();
            fixSmallEdges(shape, tolerance);
        }
    }

    if( fixSmallFaces_fg )
    {
        for( auto it = TopTools_ListIteratorOfListOfShape(solids); it.More(); it.Next() )
        {
            auto& shape = it.Value();
            fixSmallFaces(shape, tolerance);
        }
    }

    if( sewFaces_fg )
    {
        for( auto it = TopTools_ListIteratorOfListOfShape(solids); it.More(); it.Next() )
        {
            auto& shape = it.Value();
            sewFaces(shape, tolerance);
        }
    }

    if( makeSolids_fg )
    {
        for( auto it = TopTools_ListIteratorOfListOfShape(solids); it.More(); it.Next() )
        {
            auto& shape = it.Value();
            makeSolids(shape, tolerance);
        }
    }

    vector<TopoDS_Shape> input;
    for( auto it = TopTools_ListIteratorOfListOfShape(solids); it.More(); it.Next() )
    {
        auto& shape = it.Value();
        input.push_back(shape);
    }
    ClearGeometry(geometry);
    travelShape(input, geometry);
}


int repair_ori(TopoDS_Shape& shape,int heal)
{
    // heal fix model
    bool fixDegenerated_fg, fixSmallEdges_fg, fixSmallFaces_fg, sewFaces_fg, makeSolids_fg, mydeletesmalledge_fg,
        PreserveFeatureFaces_fg;
    double tolerance = 1e-8;
    fixDegenerated_fg = (heal & 0b00001);
    fixSmallEdges_fg = (heal & 0b00010) >> 1;
    fixSmallFaces_fg = (heal & 0b00100) >> 2;
    sewFaces_fg = (heal & 0b01000) >> 3;
    makeSolids_fg = (heal & 0b10000) >> 4;
    mydeletesmalledge_fg = (heal & 0b100000) >> 5;
    PreserveFeatureFaces_fg = (heal & 0b1000000) >> 6;

    if( fixDegenerated_fg )
    {
        fixDegeneratedGeo(shape);
    }

    if( fixSmallEdges_fg )
    {
        fixSmallEdges(shape, tolerance);
    }

    if( fixSmallFaces_fg )
    {
        fixSmallFaces(shape, tolerance);
    }

    if( sewFaces_fg )
    {
        sewFaces(shape, tolerance);
    }

    if( makeSolids_fg )
    {
        makeSolids(shape, tolerance);
    }
    return 0;
}

void Read_Model(const std::string& filename, TiGER::Geometry& geometry, const readModelParameters& args)
{
    bool getname_flag = args.getgetname_flag();
    bool divideclosed_flag = args.getdivideclosed_flag();
    int heal = args.getheal();
    bool bool_flag = args.getbool_flag();
    TopoDS_Shape shape;
    TopoDS_Shape r;
    size_t dot_pos = filename.find_last_of('.');
    std::string postfix = filename.substr(dot_pos + 1, filename.length() - dot_pos - 1);
    std::unordered_map<TopoDS_Face, std::string, quilt_hash_fn, quilt_equal_fn>mp;
    int ret = 0;
    int ft = 0;
    if( postfix == "igs" || postfix == "iges" || postfix == "IGS" || postfix == "IGES" )
    {
        ft = 1;
        ret = readIgesFromFile(filename, r);
        if( ret == -1 )
        {
            return;
        }
    }
    else if( postfix == "stp" || postfix == "step" || postfix == "STP" || postfix == "STEP" )
    {
        ft = 2;
        if( bool_flag )
        {
            ret = readStepFromFile_bool(filename, r, getname_flag, mp, args.getbool_para());
        }
        else
        {
            ret = readStepFromFile(filename, r);
        }

        if( ret == -1 )
        {
            return;
        }
    }
    else
    {
        return;
    }

    if( divideclosed_flag )
    {
        cutEnclosedFaces(r, 2);
    }
    repair_ori(r,heal);
    vector<TopoDS_Shape> shapes = {r};
    travelShape(shapes, geometry);
}

} // namespace TiGER
