#include "./discretize_surface_delaunay.h"
#include "./dt2D_API.h"
#include "alias.h"
#include <vector>
using namespace std;

namespace TiGER
{
	namespace discretize_surface_delaunay
	{
		void TMeshSurf_processTriangleLoop(
			const LineMesh& Loop_in,
			const TriangleDealunay& args,
			std::function<double(const std::array<double, 3>)>& func,
			SurfaceMesh& surfmesh_out) {
			dt2d::Args dtargs;
			if (func)
				dtargs.sizingFunc = func;
			dt2d::Mesh mesh;
			mesh.V.assign(Loop_in.coord.begin(), Loop_in.coord.end());
			mesh.L.assign(Loop_in.segments.begin(), Loop_in.segments.end());
			if (API_Triangulation(mesh, dtargs)) {
				surfmesh_out.coords.assign(mesh.V.begin(), mesh.V.end());
				surfmesh_out.tris.assign(mesh.F.begin(), mesh.F.end());
			}
			else
				std::printf("Triangulation failed!\n");
			return;
		}
	}
} // namespace TiGER