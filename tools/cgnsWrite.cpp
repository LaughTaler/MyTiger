
#include "cgnslib.h"
#include "cgnsWrite.h"
#include <cstdio>
#include <cstring>
#include <sstream>
#include <iostream>
#include "alias.h"

const char* EnumToString(BoundaryCondition condition) {
	switch (condition) {
	case WALL:
		return "WALL";
	case SYMMETRY:
		return "SYMMETRY";
	case FARFIELD:
		return "FARFIELD";
	case INLET:
		return "INLET";
	case OUTLET:
		return "OUTLET";
	case PERIODIC_SYMMETRY:
		return "PERIODIC_SYMMETRY";
	case ADJACENT_GRID:
		return "ADJACENT_GRID";
	default:
		return "UNKNOWN";
	}
}

int Tool_writeMeshToCGNS(const char *filename, SurfaceMesh &sm, VolumeMesh &vm)
{
	int rcode = 1;
	int fileID, baseID, zoneID, coordID, eleID, bcID;
	char baseName[256], zoneName[256];

	/* open CGNS file for write */
	int status =  cg_open(filename, CG_MODE_WRITE, &fileID);
	if (status != 0)
	{
		//error
		std::cout << "Failed to open CGNS file." << std::endl;
		return -1;
	}

	/* create base */
	strcpy(baseName, "Base");
	strcpy(zoneName, "Zone");
	int icelldim = 3, iphydim = 3;
	status = cg_base_write(fileID, baseName, icelldim, iphydim, &baseID);
	if (status != 0)
	{
		std::cout << "Failed to write base." << std::endl;
	}

	/* create zone */ 
	cgsize_t isize[3][1];
	int npt = vm.coords.size();
	isize[0][0] = npt; // vertex size
	isize[1][0] = vm.tetras.size() + vm.pyramids.size() + vm.prisms.size(); // cell size
	isize[2][0] = 0; // boundary vertex size(zero if elements not sorted)
	cg_zone_write(fileID, baseID, zoneName, isize[0], Unstructured, &zoneID);

	/* write grid coordinates */
	double* x = nullptr, * y = nullptr, * z = nullptr;
	
	x = new double[npt];
	if (!x)
	{
		std::cout << "ERROR: No enough space (x)" << std::endl;
		rcode = -1;
	}

	y = new double[npt];
	if (!y)
	{
		std::cout << "ERROR: No enough space (y)" << std::endl;
		rcode = -1;
	}

	z = new double[npt];
	if (!z)
	{
		std::cout << "ERROR: No enough space (z)" << std::endl;
		rcode = -1;
	}

	int idx = 0;
	for (int i = 0; i < npt; i++)
	{
		x[i] = vm.coords[i][0];
		y[i] = vm.coords[i][1];
		z[i] = vm.coords[i][2];
		++idx;
	}

	cg_coord_write(fileID, baseID, zoneID, CGNS_ENUMV(RealDouble), "CoordinateX", x, &coordID);
	cg_coord_write(fileID, baseID, zoneID, CGNS_ENUMV(RealDouble), "CoordinateY", y, &coordID);
	cg_coord_write(fileID, baseID, zoneID, CGNS_ENUMV(RealDouble), "CoordinateZ", z, &coordID);

	CGNS_ENUMT(DataType_t) data_type = RealDouble;
	int data_dimension = 2;
	cgsize_t dimension_vectorc[] = { sm.point_normal.size(), 3 };
	cg_array_write("PointNormal", data_type, data_dimension, dimension_vectorc, sm.point_normal.data());

	/* write element connectivity */
	int mshtype[] = { 4, 5, 6 ,8}; // tetrahedron, prism, pyramid
    ElementType_t eletype[] = {TETRA_4, PYRA_5, PENTA_6, HEXA_8};
    const char* typenames[] = {"tetra", "pyramid", "prism", "hexa"};
	int type, cnt = 0;
	cgsize_t *ele = NULL;
	int nelem_start = 0, nelem_end = 0, nbdyelem = 0;
    ele = new cgsize_t[4 * vm.tetras.size() + 5 * vm.pyramids.size() + 6 * vm.prisms.size() + 8 * vm.hexs.size()];
	for (int i = 0; i < 3; i++)
	{
		nelem_start = nelem_end + 1;
		type = mshtype[i];
		cnt = 0;
		int eleSize;
		if (type == 4) eleSize = vm.tetras.size();
		else if (type == 5) eleSize = vm.pyramids.size();
		else if (type == 6) eleSize = vm.prisms.size();
        else if( type == 8 )
            eleSize = vm.hexs.size();
		for (int j = 0; j < eleSize; j++)
		{
            if( type == 6 ) // prism
            {
                for( int k = 0; k < type; k++ )
                {
                    ele[type * cnt + k] = vm.hexs[j][k] + 1;
                    if( ele[type * cnt + k] > vm.coords.size() )
                        printf("type: %d eidx: %d idx:%d\n", type, j, k);
                }
            }
			if (type == 6)// prism
			{
				ele[type * cnt + 0] = vm.prisms[j][0] + 1;
				ele[type * cnt + 1] = vm.prisms[j][1] + 1;
				ele[type * cnt + 2] = vm.prisms[j][2] + 1;
				ele[type * cnt + 3] = vm.prisms[j][3] + 1;
				ele[type * cnt + 4] = vm.prisms[j][4] + 1;
				ele[type * cnt + 5] = vm.prisms[j][5] + 1;
			}
			else if(type == 5)
			{
				for (int k = 0; k < type; k ++)
				{
					ele[type * cnt + k] = vm.pyramids[j][k] + 1;
					if (ele[type * cnt + k] > vm.coords.size())
						printf("type: %d eidx: %d idx:%d\n", type, j, k);
				}
			}
			else if (type == 4)
			{
				for (int k = 0; k < type; k++)
				{
					ele[type * cnt + k] = vm.tetras[j][k] + 1;
					if (ele[type * cnt + k] > vm.coords.size())
						printf("type: %d eidx: %d idx:%d\n", type, j, k);
				}
			}
			cnt++;
		}
		nelem_end = nelem_start + cnt - 1;

		//printf("s: %d, e: %d, cnt: %d\n", nelem_start, nelem_end, cnt);

		cg_section_write(fileID, baseID, zoneID, typenames[i], eletype[i], nelem_start, nelem_end, nbdyelem, ele, &eleID);
	}

	/* write surface elements: triangle */
	int* bcrange = NULL;
	cgsize_t* bc_node = NULL;
	int bc_start, bc_end;
	bc_node = new cgsize_t[3 * sm.tris.size() + 4 * sm.quads.size()];
	bcrange = new int[2 * sm.bc.size() + 2]; 
	char* prefix = "_tri_";
	if (!bc_node)
	{
		std::cout << "ERROR: No enough space (bc_node)" << std::endl;
		rcode = -1;
	}
	bc_start = bc_end = nelem_end;
	bcrange[0] = bc_end;
	for (int i = 0; i < sm.bc.size(); i++)
	{
		char strbcname[256];
		std::ostringstream oss;
		oss << EnumToString(sm.bc[i]) << prefix << i;
		std::string str = oss.str();
		if (str.length() < 256) {
			strcpy(strbcname, str.c_str());
		}
		else {
			std::cerr << "strbcname too long" << std::endl;
			break;
		}
		bc_start = bc_end + 1;
		cnt = 0;
		for (int j = 0; j < sm.tris.size(); j++)
		{
			if (sm.attribute_int[j] == i)
			{
				bc_node[cnt * 3 + 0] = sm.tris[j][0] + 1;
				bc_node[cnt * 3 + 2] = sm.tris[j][1] + 1;
				bc_node[cnt * 3 + 1] = sm.tris[j][2] + 1;
				cnt++;
			}
		}
		bc_end = bc_start + cnt - 1;
		bcrange[i + 1] = bc_end;
		if(bc_end >= bc_start)
			cg_section_write(fileID, baseID, zoneID, strbcname, TRI_3, bc_start, bc_end, nbdyelem, bc_node, &eleID);
		if (cnt > 0)
		{
			std::cout << "strbcname: " << strbcname << std::endl;
			printf("s: %d e:%d c:%d\n", bc_start, bc_end, cnt);
			printf("indexe: %d\n", eleID);
		}
	}

	//write surface elements: quad
	prefix = "_quad_";
	for (int i = 0; i < sm.bc.size(); i++)
	{
		char strbcname[256];
		std::ostringstream oss;
		oss << EnumToString(sm.bc[i]) << prefix << i;
		std::string str = oss.str();
		if (str.length() < 256) {
			strcpy(strbcname, str.c_str());
		}
		else {
			std::cerr << "strbcname too long" << std::endl;
			break;
		}
		bc_start = bc_end + 1;
		cnt = 0;
		for (int j = 0; j < sm.quads.size(); j++)
		{
			if (sm.attribute_int[sm.tris.size() + j] == i)
			{
				bc_node[cnt * 4 + 0] = sm.quads[j][0] + 1;
				bc_node[cnt * 4 + 1] = sm.quads[j][1] + 1;
				bc_node[cnt * 4 + 2] = sm.quads[j][2] + 1;
				bc_node[cnt * 4 + 3] = sm.quads[j][3] + 1;
				cnt++;
			}
		}
		bc_end = bc_start + cnt - 1;
		bcrange[sm.bc.size() + i + 1] = bc_end;
		if(bc_end >= bc_start)
			cg_section_write(fileID, baseID, zoneID, strbcname, QUAD_4, bc_start, bc_end, nbdyelem, bc_node, &eleID);

		if (cnt > 0)
		{
			std::cout << "strbcname: " << strbcname << std::endl;
			printf("s: %d e:%d c:%d\n", bc_start, bc_end, cnt);
			printf("indexe: %d\n", eleID);
		}
	}

	cgsize_t dimension_vectorf[] = { sm.facet_normal.size(), 3 };
	cg_array_write("FacetNormal", data_type, data_dimension, dimension_vectorf, sm.facet_normal.data());

	//write boundary
	cgsize_t* bc_elem = NULL;
	bc_elem = new cgsize_t[sm.attribute_int.size()];
	int bccnt = 0;
	for (int i = 0; i < sm.bc.size(); i++)
	{
		char strbcname[256];
		memset(strbcname, 0, sizeof(strbcname));
        strcpy(strbcname, EnumToString(sm.bc[i]));
        strcat(strbcname, "_tri");
		cnt = 0;
		int istart, iend;
		istart = bcrange[i] + 1;
		iend = bcrange[i + 1];
		for (int j = istart; j <= iend; j++)
		{
			bc_elem[cnt] = j;
			cnt++;
		}

		if (iend >= istart)
		{
            if( sm.bc[i] == 0 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCWall, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 1 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCSymmetryPlane, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 2 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCFarfield, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 3 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCInflow, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 4 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCOutflow, PointList, cnt, bc_elem, &bcID);
			bccnt++;
		}
	}

	int ibc = sm.bc.size();
	for (int i = 0; i < sm.bc.size(); i++)
	{
		char strbcname[256];
		memset(strbcname, 0, sizeof(strbcname));
        strcpy(strbcname, EnumToString(sm.bc[i]));
        strcat(strbcname, "_quad");
		cnt = 0;
		int istart, iend;
		istart = bcrange[ibc] + 1;
		iend = bcrange[ibc + 1];
        ibc ++;
		for (int j = istart; j <= iend; j++)
		{
			bc_elem[cnt] = j;
			cnt++;
		}

		if (iend >= istart)
		{
            if( sm.bc[i] == 0 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCWall, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 1 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCSymmetryPlane, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 2 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCFarfield, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 3 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCInflow, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 4 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCOutflow, PointList, cnt, bc_elem, &bcID);
			bccnt++;
		}
	}

	for (int i = 0; i < bccnt; i++)
	{
		cg_goto(fileID, baseID, "Zone_t", zoneID, "ZoneBC_t", 1, "BC_t", i + 1, "end");
		cg_gridlocation_write(FaceCenter);
	}

	cg_close(fileID);
	if (x) delete []x;
	if (y) delete []y;
	if (z) delete []z;
	if (ele) delete[]ele;
	if (bc_node) delete []bc_node;
	if (bcrange) delete []bcrange;
	if (bc_elem) delete []bc_elem;

	return 0;
}

int Tool_writeMeshToCGNS(const char* filename, SurfaceMesh& sm)
{
	int rcode = 1;
	int fileID, baseID, zoneID, coordID, eleID, bcID;
	char baseName[256], zoneName[256];

	/* open CGNS file for write */
	cg_open(filename, CG_MODE_WRITE, &fileID);

	/* create base */
	strcpy(baseName, "Base");
	strcpy(zoneName, "Zone");
	int icelldim = 2, iphydim = 3;
	cg_base_write(fileID, baseName, icelldim, iphydim, &baseID);

	/* create zone */
	cgsize_t isize[3][1];
	int npt = sm.coords.size();
	isize[0][0] = npt; // vertex size
	isize[1][0] = sm.attribute_int.size(); // cell size
	isize[2][0] = 0; // boundary vertex size(zero if elements not sorted)
	cg_zone_write(fileID, baseID, zoneName, isize[0], Unstructured, &zoneID);

	/* write grid coordinates */
	double* x = nullptr, * y = nullptr, * z = nullptr;

	x = new double[npt];
	if (!x)
	{
		std::cout << "ERROR: No enough space (x)" << std::endl;
		rcode = -1;
	}

	y = new double[npt];
	if (!y)
	{
		std::cout << "ERROR: No enough space (y)" << std::endl;
		rcode = -1;
	}

	z = new double[npt];
	if (!z)
	{
		std::cout << "ERROR: No enough space (z)" << std::endl;
		rcode = -1;
	}
	
	int idx = 0;
	for (int i = 0; i < npt; i++)
	{
		x[i] = sm.coords[i][0];
		y[i] = sm.coords[i][1];
		z[i] = sm.coords[i][2];
		++idx;
	}

	cg_coord_write(fileID, baseID, zoneID, CGNS_ENUMV(RealDouble), "CoordinateX", x, &coordID);
	cg_coord_write(fileID, baseID, zoneID, CGNS_ENUMV(RealDouble), "CoordinateY", y, &coordID);
	cg_coord_write(fileID, baseID, zoneID, CGNS_ENUMV(RealDouble), "CoordinateZ", z, &coordID);

	CGNS_ENUMT(DataType_t) data_type = RealDouble;
	int data_dimension = 2;
	cgsize_t dimension_vectorc[] = { sm.point_normal.size(), 3 };
	cg_array_write("PointNormal", data_type, data_dimension, dimension_vectorc, sm.point_normal.data());

	/* write surface element: triangle*/
	int* bcrange = NULL;
	cgsize_t* bc_node = NULL;
	int bc_start, bc_end;
	bc_node = new cgsize_t[3 * sm.tris.size() + 4 * sm.quads.size()];
	bcrange = new int[2 * sm.bc.size() + 2];
	bc_start = bc_end = 0;
	bcrange[0] = bc_end;
	int cnt = 0, nbdyelem = 0;
	char* prefix = "_tri_";
	for (int i = 0; i < sm.bc.size(); i++)
	{
		char strbcname[256];
		std::ostringstream oss;
        oss << EnumToString(sm.bc[i]) << prefix << i;
		std::string str = oss.str();
		if (str.length() < 256) {
			strcpy(strbcname, str.c_str());
		}
		else {
			std::cerr << "strbcname too long" << std::endl;
			break;
		}
		bc_start = bc_end + 1;
		cnt = 0;
		for (int j = 0; j < sm.tris.size(); j++)
		{
			if (sm.attribute_int[j] == i)
			{
				bc_node[cnt * 3 + 0] = sm.tris[j][0] + 1;
				bc_node[cnt * 3 + 2] = sm.tris[j][1] + 1;
				bc_node[cnt * 3 + 1] = sm.tris[j][2] + 1;
				cnt++;
			}
		}
		bc_end = bc_start + cnt - 1;
		bcrange[i + 1] = bc_end;
		if (bc_end >= bc_start)
			cg_section_write(fileID, baseID, zoneID, strbcname, TRI_3, bc_start, bc_end, nbdyelem, bc_node, &eleID);
		if (cnt > 0)
		{
			std::cout << "strbcname: " << strbcname << std::endl;
			printf("s: %d e:%d c:%d\n", bc_start, bc_end, cnt);
			printf("indexe: %d\n", eleID);
		}
	}

	/* write surface element: quad */
	prefix = "_quad_";
	for (int i = 0; i < sm.bc.size(); i++)
	{
		char strbcname[256];
		std::ostringstream oss;
        oss << EnumToString(sm.bc[i]) << prefix << i;
		std::string str = oss.str();
		if (str.length() < 256) {
			strcpy(strbcname, str.c_str());
		}
		else {
			std::cerr << "strbcname too long" << std::endl;
			break;
		}
		bc_start = bc_end + 1;
		cnt = 0;
		for (int j = 0; j < sm.quads.size(); j++)
		{
			if (sm.attribute_int[j] == i)
			{
				bc_node[cnt * 4 + 0] = sm.quads[j][0] + 1;
				bc_node[cnt * 4 + 1] = sm.quads[j][1] + 1;
				bc_node[cnt * 4 + 2] = sm.quads[j][2] + 1;
				bc_node[cnt * 4 + 3] = sm.quads[j][3] + 1;
				cnt++;
			}
		}
		bc_end = bc_start + cnt - 1;
		bcrange[sm.bc.size() + i + 1] = bc_end;
		if (bc_end >= bc_start)
			cg_section_write(fileID, baseID, zoneID, strbcname, QUAD_4, bc_start, bc_end, nbdyelem, bc_node, &eleID);
		if (cnt > 0)
			printf("  bc[%d]: %-16s global elem_id: [%d, %d], count: %d\n", eleID, strbcname, bc_start, bc_end, cnt);
	}

	/* write normal */
	cgsize_t dimension_vectorf[] = { sm.facet_normal.size(), 3 };
	cg_array_write("FacetNormal", data_type, data_dimension, dimension_vectorf, sm.facet_normal.data());

	/* write boundary */
	cgsize_t* bc_elem = NULL;
	bc_elem = new cgsize_t[sm.attribute_int.size()];
	int bccnt = 0;
	for (int i = 0; i < sm.bc.size(); i++)
	{
		char strbcname[256];
		memset(strbcname, 0, sizeof(strbcname));
        strcpy(strbcname, EnumToString(sm.bc[i]));
        strcat(strbcname, "_tri");
		cnt = 0;
		int istart, iend;
		istart = bcrange[i] + 1;
		iend = bcrange[i + 1];
		for (int j = istart; j <= iend; j++)
		{
			bc_elem[cnt] = j;
			cnt++;
		}

		if (iend >= istart)
		{
            if( sm.bc[i] == 0 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCWall, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 1 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCSymmetryPlane, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 2 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCFarfield, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 3 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCInflow, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 4 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCOutflow, PointList, cnt, bc_elem, &bcID);
			bccnt++;
		}
	}

	int ibc = sm.bc.size();
	for (int i = 0; i < sm.bc.size(); i++)
	{
		char strbcname[256];
		memset(strbcname, 0, sizeof(strbcname));
        strcpy(strbcname, EnumToString(sm.bc[i]));
		strcat(strbcname, "_quad");
		cnt = 0;
		int istart, iend;
		istart = bcrange[ibc] + 1;
		iend = bcrange[ibc + 1];
        ibc++;
		for (int j = istart; j <= iend; j++)
		{
			bc_elem[cnt] = j;
			cnt++;
		}

		if (iend >= istart)
		{
            if( sm.bc[i] == 0 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCWall, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 1 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCSymmetryPlane, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 2 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCFarfield, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 3 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCInflow, PointList, cnt, bc_elem, &bcID);
            if( sm.bc[i] == 4 )
				cg_boco_write(fileID, baseID, zoneID, strbcname, BCOutflow, PointList, cnt, bc_elem, &bcID);
			bccnt++;
		}
	}

	for (int i = 0; i < bccnt; i++)
	{
		cg_goto(fileID, baseID, "Zone_t", zoneID, "ZoneBC_t", 1, "BC_t", i + 1, "end");
		cg_gridlocation_write(FaceCenter);
	}

	cg_close(fileID);
	if (x) delete []x;
	if (y) delete []y;
	if (z) delete []z;
	if (bc_elem) delete[]bc_elem;
	if (bc_node) delete []bc_node;
	if (bcrange) delete []bcrange;

	return 0;
}