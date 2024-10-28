#include "Pyramid_surface.h"

namespace TiGER
{
void pyramid_surface(const SurfaceMesh& surfacemeshin, VolumeMesh& volumemesh, SurfaceMesh& surfacemeshout)
{
    for( auto i : surfacemeshin.coords )
    {
        surfacemeshout.coords.push_back(i);
    }

    for( auto i : surfacemeshin.tris )
    {
        surfacemeshout.tris.push_back(i);
    }

    for( auto i : surfacemeshin.quads )
    {
        double u1[3], u2[3], u3[3], u4[3], normal[3], cross[3];
        for( int j = 0; j < 3; j++ )
        {
            u1[j] = surfacemeshin.coords[i[1]][j] - surfacemeshin.coords[i[0]][j];
            u2[j] = surfacemeshin.coords[i[2]][j] - surfacemeshin.coords[i[1]][j];
            u3[j] = surfacemeshin.coords[i[3]][j] - surfacemeshin.coords[i[2]][j];
            u4[j] = surfacemeshin.coords[i[0]][j] - surfacemeshin.coords[i[3]][j];
            cross[j] = (surfacemeshin.coords[i[0]][j] + surfacemeshin.coords[i[1]][j] + surfacemeshin.coords[i[2]][j] +
                        surfacemeshin.coords[i[3]][j]) /
                       4;
        }
        normal[0] = u1[1] * u2[2] - u1[2] * u2[1];
        normal[1] = u1[2] * u2[0] - u1[0] * u2[2];
        normal[2] = u1[0] * u2[1] - u1[1] * u2[0];
        double norm, height;
        double length = 0;
        length =
            sqrt(u1[0] * u1[0] + u1[1] * u1[1] + u1[2] * u1[2]) + sqrt(u2[0] * u2[0] + u2[1] * u2[1] + u2[2] * u2[2]) +
            sqrt(u3[0] * u3[0] + u3[1] * u3[1] + u3[2] * u3[2]) + sqrt(u4[0] * u4[0] + u4[1] * u4[1] + u4[2] * u4[2]);
        height = length / 4 / 2;

        norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
        for( int j = 0; j < 3; j++ )
        {
            normal[j] = normal[j] / norm * height;
        }

        std::array<double, 3> newpoint;
        for( int j = 0; j < 3; j++ )
        {
            newpoint[j] = cross[j] + normal[j];
        }
        surfacemeshout.coords.push_back(newpoint);

        std::array<int, 5> newpyramid;
        for( int j = 0; j < 4; j++ )
        {
            newpyramid[j] = i[j];
        }
        newpyramid[4] = surfacemeshout.coords.size() - 1;
        volumemesh.pyramids.push_back(newpyramid);
    }

    for( auto i : volumemesh.pyramids )
    {
        std::array<int, 3> temptris0;
        std::array<int, 3> temptris1;
        std::array<int, 3> temptris2;
        std::array<int, 3> temptris3;

        temptris0[0] = i[0];
        temptris0[1] = i[1];
        temptris0[2] = i[4];

        temptris1[0] = i[1];
        temptris1[1] = i[2];
        temptris1[2] = i[4];

        temptris2[0] = i[2];
        temptris2[1] = i[3];
        temptris2[2] = i[4];

        temptris3[0] = i[3];
        temptris3[1] = i[0];
        temptris3[2] = i[4];

        surfacemeshout.tris.push_back(temptris0);
        surfacemeshout.tris.push_back(temptris1);
        surfacemeshout.tris.push_back(temptris2);
        surfacemeshout.tris.push_back(temptris3);
    }
    return;
}
} // namespace TiGER