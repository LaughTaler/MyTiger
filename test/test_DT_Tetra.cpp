#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do
                          // this in one cpp file
#include "alias.h"
#include "test_api.h"

#ifdef _MSC_VER
std::vector<std::string> list_files(const std::string& directory_path)
{
    std::vector<std::string> aaa;
    // 构建 PowerShell 命令
    std::string command = "powershell.exe -Command \"Get-ChildItem -Name '" + directory_path + "'\"";

    std::array<char, 128> buffer;
    std::string result;

    // 打开管道以读取命令输出
    FILE* pipe = _popen(command.c_str(), "r");
    if( !pipe )
    {
        std::cerr << "error" << command << std::endl;
        return aaa;
    }

    try
    {
        while( fgets(buffer.data(), buffer.size(), pipe) != nullptr )
        {
            aaa.push_back(buffer.data());
        }
    }
    catch( ... )
    {
        _pclose(pipe);
        throw;
    }
    _pclose(pipe);

    return aaa;
}
#else
std::vector<std::string> list_files(const std::string& directory_path)
{
    std::vector<std::string> aaa;
    std::string command = "ls " + directory_path;

    std::array<char, 128> buffer;
    std::string result;

    // 打开管道以读取命令输出
    FILE* pipe = popen(command.c_str(), "r");
    if( !pipe )
    {
        std::cerr << "error" << command << std::endl;
        return aaa;
    }

    try
    {
        while( fgets(buffer.data(), buffer.size(), pipe) != nullptr )
        {
            // 将每一行输出添加到结果字符串中
            aaa.push_back(buffer.data());
        }
    }
    catch( ... )
    {
        pclose(pipe);
        throw;
    }
    pclose(pipe);
    return aaa;
}

#endif
TEST_CASE("tiger 1.5 base", "HexGrid_createPrismaticMesh")
{
#ifdef _MSC_VER
    std::string filename = "Z:/home/sftp1/models/1000_sur";
#else
    std::string filename = "/home/sftp1/models/1000_sur";
#endif
    std::vector<std::string> filenames = list_files(filename); // Assuming list_files() returns a vector of filenames

    for( const auto& filename : filenames )
    {
        TiGER::Mesh m;
        test_load_vtk("1000_sur/" + filename, m);
        TiGER::SurfaceMesh surf;
        TiGER::matrix_to_list<3>(m.Vertex, surf.coords);
        TiGER::matrix_to_list<3>(m.Topo, surf.tris);
        TiGER::matrix_to_list<1>(m.Masks, surf.attribute_int);
        TiGER::TetrahedraParameters args;
        TiGER::VolumeMesh vol;
        // args.setrefiene(0);
        // args.setoptlevel(0);
        TiGER::tetrahedral_mesh::TetGrid_Constrained(surf, args, vol, NULL);

        REQUIRE(vol.coords.size() > 0);
        REQUIRE(vol.tetras.size() > 0);
    }
}
