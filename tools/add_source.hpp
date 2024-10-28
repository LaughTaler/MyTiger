/**
 * @author :wangyifei
*/
#include <sstream>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include "tiger.h"
#include "tiger_parameters.h"
#include "alias.h"

void add_source(std::string filename,TiGER::SizingManager &manager){
	int nsource;  
    std::ifstream file(filename);    
    if (file.is_open()){
        while(file){
            std::string line;
            std::getline(file, line);
            std::stringstream lineStream;
            lineStream << line;
            std::string keyString;
            lineStream >> keyString;
            if (keyString=="POINT"){
                lineStream >> nsource;
                for(int i=0; i<nsource; i++){
                    std::string subLine;
                    std::getline(file, subLine);
					std::stringstream sublineStream;
					sublineStream << subLine;
					std::array<double,3> xyz;
                    TiGER::PointSourceParameters args;
					double tmp_hmin,tmp_beta;
					sublineStream >> tmp_hmin >> tmp_beta >> xyz[0] >> xyz[1] >> xyz[2];
                    args.setsizing_value(tmp_hmin);
                    args.setbeta(tmp_beta);
					TiGER::size_field::SizingFunction_addPointSource(xyz,args,manager);
                }
            }
            else if (keyString == "LINE"){
                lineStream >> nsource;
                for(int i=0; i<nsource; i++){
                    std::string subLine;
                    std::getline(file, subLine);
					std::stringstream sublineStream;
                    std::array<std::array<double,3>,2> xyz;
					sublineStream << subLine;
                    TiGER::LineSourceParameters args;
					double tmp_hmin,tmp_beta;
					sublineStream >> tmp_hmin >> tmp_beta >> xyz[0][0] >> xyz[0][1] >> xyz[0][2] >> xyz[1][0] >> xyz[1][1] >> xyz[1][2];
					args.setsizing_value(tmp_hmin);
                    args.setbeta(tmp_beta);
					TiGER::size_field::SizingFunction_addLineSource(xyz,args,manager);
                }
            }
            else if (keyString == "TRIANGLE"){
                lineStream >> nsource;
                for(int i=0; i<nsource; i++){
                    std::string subLine;
                    std::getline(file, subLine);
					std::stringstream sublineStream;
                    std::array<std::array<double,3>,3> xyz;
					sublineStream << subLine;
                    TiGER::TriangleSourceParameters args;
					double tmp_hmin,tmp_beta;
					sublineStream >> tmp_hmin >> tmp_beta >> xyz[0][0] >> xyz[0][1] >> xyz[0][2] >> xyz[1][0] >> xyz[1][1] >> xyz[1][2] >> xyz[2][0] >> xyz[2][1] >> xyz[2][2];
					args.setsizing_value(tmp_hmin);
                    args.setbeta(tmp_beta);
					TiGER::size_field::SizingFunction_addTriangleSource(xyz,args,manager);
                }
            }
			else if (keyString == "CUBOID"){
                lineStream >> nsource;
                for(int i=0; i<nsource; i++){
                    std::string subLine;
                    std::getline(file, subLine);
					std::stringstream sublineStream;
                    std::array<std::array<double,3>,8> xyz;
					sublineStream << subLine;
                    TiGER::CubicSourceParameters args;
					double tmp_hmin,tmp_beta;
					sublineStream >> tmp_hmin >> tmp_beta >> xyz[0][0] >> xyz[0][1] >> xyz[0][2]
                                                          >> xyz[1][0] >> xyz[1][1] >> xyz[1][2]
                                                          >> xyz[2][0] >> xyz[2][1] >> xyz[2][2]
                                                          >> xyz[3][0] >> xyz[3][1] >> xyz[3][2]
                                                          >> xyz[4][0] >> xyz[4][1] >> xyz[4][2]
                                                          >> xyz[5][0] >> xyz[5][1] >> xyz[5][2]
                                                          >> xyz[6][0] >> xyz[6][1] >> xyz[6][2]
                                                          >> xyz[7][0] >> xyz[7][1] >> xyz[7][2];
					args.setsizing_value(tmp_hmin);
                    args.setbeta(tmp_beta);
					TiGER::size_field::SizingFunction_addCubicSource(xyz,args,manager);
                }
            }

        }
        file.close();
		std::cout << "Finish create source!" << std::endl;
    }
	else{
		std::cout << "Failed create source!" << std::endl;
	}
	
}
