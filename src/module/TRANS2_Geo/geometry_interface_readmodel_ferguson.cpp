// #include <algorithm>
// #include <cmath>
// #include <iostream>
// #include <map>
// #include <Eigen/Dense>
// #include<memory>

// #include "../../include/geometry_data_struct.h"
// #include "../../include/geometry_interface_ferguson.h"
// #include "./geometry_occ.h"
// #include <GCPnts_AbscissaPoint.hxx>


// void readModel_Ferguson(const std::string& filename, TiGER::Geometry& geometry,
//                         const readModelParameters& args) {
//   bool getname_flag = args.getgetname_flag();
//   bool divideclosed_flag = args.getdivideclosed_flag();
//   int heal = args.getheal();
//   TopoDS_Shape shape;
//   TopoDS_Shape r;
//   size_t dot_pos = filename.find_last_of('.');
//   std::string postfix =
//       filename.substr(dot_pos + 1, filename.length() - dot_pos - 1);
//   int ret = 0;
//   int ft = 0;
//   std::unordered_map<TopoDS_Face, std::string, quilt_hash_fn, quilt_equal_fn>
//       mp;
//   if (postfix == "igs" || postfix == "iges" || postfix == "IGS" ||
//       postfix == "IGES") {
//     ft = 1;
//     ret = readIgesFromFile(filename, r);
//     if (ret == -1) {
//       return;
//     }
//   }
//   else if (postfix == "stp" || postfix == "step" || postfix == "STP" ||
//              postfix == "STEP") {
//     ft = 2;
//     ret = readStepFromFile_bool(filename, r, getname_flag, mp);

//     if (ret == -1) {
//       return;
//     }
//   }
//   else {
//     return;
//   }

//    if (divideclosed_flag) {
//     ShapeUpgrade_ShapeDivideClosed shapeDivideTool(r);
//     shapeDivideTool.SetNbSplitPoints(2);
//     shapeDivideTool.Perform();
//     shape = shapeDivideTool.Result();
//   } else {
//     shape = r;
//   }
//   TopExp_Explorer exp_model, exp_quilt, exp_wire, exp_curve;
//   for (exp_model.Init(shape, TopAbs_SOLID); exp_model.More();
//        exp_model.Next()) {
//     std::vector<int> temp_surf_on_body;
//     TopoDS_Solid topo_body = TopoDS::Solid(exp_model.Current());
//     for (exp_quilt.Init(topo_body, TopAbs_FACE); exp_quilt.More();
//          exp_quilt.Next()) {
//       TopoDS_Face topo_face = TopoDS::Face(exp_quilt.Current());
//       BRepAdaptor_Surface surface(topo_face);
//       std::vector<int> m_uv(2);
//       int intervalu = 43;  // round(8 * ((surface.LastUParameter() -
//                            // surface.FirstUParameter()) / M_PI));
//       int intervalv = 43;  // round(8 * ((surface.LastVParameter() -
//                            // surface.FirstVParameter()) / M_PI));
//       std::vector<vec3d> m_uv_pts;
//       m_uv[0] = intervalu + 1;
//       m_uv[1] = intervalv + 1;
//       for (int u = 0; u <= intervalu; u++) {
//         for (int v = 0; v <= intervalv; v++) {
//           gp_Pnt p;
//           surface.D0(
//               u * 1.0 / intervalu *
//                       (surface.LastUParameter() - surface.FirstUParameter()) +
//                   surface.FirstUParameter(),
//               v * 1.0 / intervalv *
//                       (surface.LastVParameter() - surface.FirstVParameter()) +
//                   surface.FirstVParameter(),
//               p);
//           vec3d temp = {p.X(), p.Y(), p.Z()};
//           m_uv_pts.push_back(temp);
//         }
//       }

//       std::shared_ptr<TiGER::GeometrySurface_Ferguson> surface_temp =
//           std::make_shared<TiGER::GeometrySurface_Ferguson>(m_uv_pts, m_uv[0],
//                                                      m_uv[1]);
//       int surfaceid = geometry.surfaces_.size();

//       int idd = geometry.surfaces_.size();
//       geometry.surfaces_.push_back(surface_temp);
//       if (mp.count(topo_face) != 0) {
//         geometry.mp_index_name[idd] = mp[topo_face];
//       }

//       temp_surf_on_body.push_back(surfaceid);
//       std::vector<int> temp_curve_on_face;

//       for (exp_curve.Init(topo_face, TopAbs_EDGE); exp_curve.More();
//            exp_curve.Next()) {
//         TopoDS_Edge topo_edge = TopoDS::Edge(exp_curve.Current());
//         BRepAdaptor_Curve curve(topo_edge);
//         int interval = 43;
//         std::vector<vec3d> curve_pts;
//         for (int para = 0; para <= interval; para++) {
//           gp_Pnt p;
//           curve.D0(para * 1.0 *
//                            (curve.LastParameter() - curve.FirstParameter()) /
//                            interval +
//                        curve.FirstParameter(),
//                    p);
//           vec3d temp = {p.X(), p.Y(), p.Z()};
//           curve_pts.push_back(temp);
//         }

//         std::shared_ptr<GeometryCurve_Ferguson> curve_temp =
//             std::make_shared<GeometryCurve_Ferguson>(curve_pts);
//         int curveid = geometry.curves_.size();
//         geometry.curves_.push_back(curve_temp);
//         temp_curve_on_face.push_back(curveid);
//       }
//       geometry.topo_surf_to_curves_.push_back(temp_curve_on_face);
//     }
//     geometry.topo_body_to_surfs_.push_back(temp_surf_on_body);
//   }

// }
