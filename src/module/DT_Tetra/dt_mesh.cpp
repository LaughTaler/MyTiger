#include "./dt_API.h"
#include "./tetrahedral_mesh.h"
#include "./Pyramid_surface.h"
using namespace std;
namespace TiGER
{
namespace tetrahedral_mesh
{
void constrainDelaunay(const SurfaceMesh& surface_mesh, const TetrahedraParameters& args, VolumeMesh& volume_mesh,
                       const SizingFunction& func)
{
    dt::Args dtargs;
    dtargs.infolevel = args.getinfolevel();
    dtargs.constrain = 1;
    dtargs.refine = args.getrefine();
    dtargs.optlevel = args.getoptlevel();
    dtargs.optloop = args.getoptloop();
    dtargs.nthread = args.getnthread();
    dtargs.size = args.getsize();
    dtargs.decay = args.getdecay();
    dtargs.minEdge = args.getminEdge();
    dtargs.maxEdge = args.getmaxEdge();
    dtargs.minSize = args.getminSize();
    dtargs.maxSize = args.getmaxSize();
    dtargs.backgroundSpace = args.getbackgroundSpace();
    dtargs.growsize = args.getgrowsize();
    dtargs.optangle = args.getoptangle();
    dtargs.optratio = args.getoptratio();
    dtargs.hole = args.gethole();
    dtargs.layer = args.getlayer();
    if( func )
        dtargs.sizingFunc = func;

    dt::Mesh mesh;
    mesh.V.assign(surface_mesh.coords.begin(), surface_mesh.coords.end());
    mesh.F.assign(surface_mesh.tris.begin(), surface_mesh.tris.end());
    if( surface_mesh.point_attribute_double.size() != 0 )
    {
        if( surface_mesh.point_attribute_double.size() == surface_mesh.coords.size() )
        {
            mesh.faceID.resize(surface_mesh.point_attribute_double.size());
            for( int i = 0; i < surface_mesh.point_attribute_double.size(); i++ )
                mesh.pointSize[i] = surface_mesh.point_attribute_double[i];
        }
        else
        {
            std::printf("Error before DT, error input point attribute information!\n");
        }
    }
    if( API_Tetrahedralize(mesh, dtargs) )
    {
        volume_mesh.coords.assign(mesh.V.begin(), mesh.V.end());
        volume_mesh.tetras.assign(mesh.T.begin(), mesh.T.end());
        volume_mesh.regions.resize(mesh.tetID.size());
        for( int i = 0; i < mesh.tetID.size(); i++ )
        {
            volume_mesh.regions[i] = mesh.tetID[i];
        }
    }
    else
        std::printf("tetrahedralize failed!\n");
    return;
}

void TetGrid_Delaunay_QuadInputs(const SurfaceMesh& surface_mesh, const TetrahedraParameters& args,
                                 VolumeMesh& volume_mesh, const SizingFunction& func)
{
    dt::Args dtargs;
    dtargs.infolevel = args.getinfolevel();
    dtargs.constrain = 1;
    dtargs.refine = args.getrefine();
    dtargs.optlevel = args.getoptlevel();
    dtargs.optloop = args.getoptloop();
    dtargs.nthread = args.getnthread();
    dtargs.size = args.getsize();
    dtargs.decay = args.getdecay();
    dtargs.minEdge = args.getminEdge();
    dtargs.maxEdge = args.getmaxEdge();
    dtargs.minSize = args.getminSize();
    dtargs.maxSize = args.getmaxSize();
    dtargs.backgroundSpace = args.getbackgroundSpace();
    dtargs.growsize = args.getgrowsize();
    dtargs.optangle = args.getoptangle();
    dtargs.optratio = args.getoptratio();
    dtargs.hole = args.gethole();
    dtargs.layer = args.getlayer();
    if( func )
        dtargs.sizingFunc = func;

    SurfaceMesh surface_mesh_out;
    pyramid_surface(surface_mesh, volume_mesh, surface_mesh_out);

    dt::Mesh mesh;
    mesh.V.assign(surface_mesh_out.coords.begin(), surface_mesh_out.coords.end());
    mesh.F.assign(surface_mesh_out.tris.begin(), surface_mesh_out.tris.end());

    if( API_Tetrahedralize(mesh, dtargs) )
    {
        volume_mesh.coords.assign(mesh.V.begin(), mesh.V.end());
        volume_mesh.tetras.assign(mesh.T.begin(), mesh.T.end());
        volume_mesh.regions.resize(mesh.tetID.size());
        for( int i = 0; i < mesh.tetID.size(); i++ )
        {
            volume_mesh.regions[i] = mesh.tetID[i];
        }
    }
    else
        std::printf("tetrahedralize failed!\n");
    return;
}

void TetGrid_Opt(const VolumeMesh& volume_mesh_in, const TetrahedraParameters& args, VolumeMesh& volume_mesh_out,
                 const std::function<std::array<double, 9>(std::array<double, 3>)>& func)
{
    dt::Args dtargs;
    dtargs.infolevel = args.getinfolevel();
    dtargs.constrain = 1;
    dtargs.refine = args.getrefine();
    dtargs.optlevel = args.getoptlevel();
    dtargs.optloop = args.getoptloop();
    dtargs.nthread = args.getnthread();
    dtargs.size = args.getsize();
    dtargs.decay = args.getdecay();
    dtargs.minEdge = args.getminEdge();
    dtargs.maxEdge = args.getmaxEdge();
    dtargs.minSize = args.getminSize();
    dtargs.maxSize = args.getmaxSize();
    dtargs.backgroundSpace = args.getbackgroundSpace();
    dtargs.growsize = args.getgrowsize();
    dtargs.optangle = args.getoptangle();
    dtargs.optratio = args.getoptratio();
    dtargs.hole = args.gethole();
    dt::Mesh mesh;
    // mesh.V.assign(surface_mesh.coords.begin(), surface_mesh.coords.end());
    // mesh.F.assign(surface_mesh.tris.begin(), surface_mesh.tris.end());
    // if (surface_mesh.point_attribute_double.size() != 0) {
    //	if (surface_mesh.point_attribute_double.size() ==
    // surface_mesh.coords.size()) {
    //		mesh.I.resize(surface_mesh.point_attribute_double.size());
    //		for (int i = 0; i < surface_mesh.point_attribute_double.size();
    // i++) 			mesh.I[i].pointSize =
    // surface_mesh.point_attribute_double[i];
    //	}
    //	else {
    //		std::printf("Error before DT, error input point attribute
    // information!\n");
    //	}
    // }
    // if (API_Tetrahedralize(mesh, dtargs)) {
    //	volume_mesh.coords.assign(mesh.V.begin(), mesh.V.end());
    //	volume_mesh.tetras.assign(mesh.T.begin(), mesh.T.end());
    // }
    // else
    //	std::printf("tetrahedralize failed!\n");
    return;
}
} // namespace tetrahedral_mesh
} // namespace TiGER