#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <memory>
#include <Eigen/Dense>
#include <unordered_map>
#include "sizing_field_interface.h"
#include "mesh_data_struct.h"
// #include <igl/AABB.h>
#include "../../sizingfunction/src/tiger_sizingfunction.h"
#include "../tools/binary_tree.hpp"
#include "alias.h"

#define NEWSOURCE
#ifdef NEWSOURCE
namespace TiGER {
    double backgroundsizefield::getsize(const double& x, const double& y, const double& z){
        // Eigen::MatrixXd sqr;
        // Eigen::MatrixXi I;
        // Eigen::MatrixXd C;
        // tree.squared_distance(V, F, P, sqr, I, C);
        // double dis = (P-C).norm();
        // Eigen::MatrixXd weights;
        // Eigen::MatrixXi f = F.row(I(0,0));
        // igl::barycentric_coordinates(C, V.row(f(0,0)), V.row(f(0,1)), V.row(f(0,2)), weights);
        // double s = 0;
        // for (int i=0; i<3; i++){
        //     s += backgroundmesh.point_attribute_double[f(0, i)] * weights(0, i);
        // }
        
        Eigen::RowVector3d P;
        P(0, 0) = x;
        P(0, 1) = y;
        P(0, 2) = z;
        size_t simplex_id;
        double distance;
        Eigen::Vector3d barycentric_coordinates;
        Eigen::RowVector3d closest_point;
        octree_.queryNearestTriangle(P, simplex_id, distance, closest_point,
                               barycentric_coordinates, -1);

        double s = 0;
        size_t p0 = backgroundmesh.tris[simplex_id][0];
        size_t p1 = backgroundmesh.tris[simplex_id][1];
        size_t p2 = backgroundmesh.tris[simplex_id][2];

        double s0 = backgroundmesh.point_attribute_double[p0];
        double s1 = backgroundmesh.point_attribute_double[p1];
        double s2 = backgroundmesh.point_attribute_double[p2];
        s = barycentric_coordinates.x() * s0 +
               barycentric_coordinates.y() * s1 +
               barycentric_coordinates.z() * s2;
        if (s<=0){
            std::cout << "Warning background mesh" << std::endl;
            s = 0;
        }
        s += distance*0.2;
        return s;
    }

    // void backgroundsizefield::exportSizeFieldToVTK(const std::string &filePath){
    //     FILE *file = fopen(filePath.c_str(), "w");

    //     if (file != nullptr) {
    //         fprintf(file, "# vtk DataFile Version 3.0\n");
    //         fprintf(file, "vtk output\n");
    //         fprintf(file, "ASCII\n");
    //         fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

    //         int numPnt = backgroundmesh.coords.size();
    //         fprintf(file, "POINTS %d double\n", numPnt);
    //         for (int i = 0; i < numPnt; i++) {
    //             fprintf(file, "%.15f %.15f %.15f\n", backgroundmesh.coords[i][0], backgroundmesh.coords[i][1], backgroundmesh.coords[i][2]);
    //         }

    //         int numTri = backgroundmesh.tris.size();
    //         fprintf(file, "CELLS %d %d\n", numTri, 4 * numTri);
    //         for (int i = 0; i < numTri; i++) {
    //             fprintf(file, "3 %d %d %d\n", backgroundmesh.tris[i][0], backgroundmesh.tris[i][2], backgroundmesh.tris[i][2]);
    //         }

    //         fprintf(file, "CELL_TYPES %d\n", numTri);
    //         for (int i = 0; i < numTri; i++) {
    //             fprintf(file, "5\n");
    //         }

    //         if (backgroundmesh.point_attribute_double.size() == numPnt) {
    //             fprintf(file, "POINT_DATA %d\n", numPnt);
    //             fprintf(file, "SCALARS sizing_value double 1\n");
    //             fprintf(file, "LOOKUP_TABLE default\n");
    //             for (int i = 0; i < numPnt; i++) {
    //                 fprintf(file, "%.15f\n", backgroundmesh.point_attribute_double[i]);
    //             }
    //         }
    //         fclose(file);
    //     } 
    //     else {
    //         perror("Unable to open file");
    //     }

        
    // }

    double constsizefield::getsize(const double& x, const double& y, const double& z) {
        return size;
    }

    double pointSource::getsize(const double& x, const double& y, const double& z) {
        double size;
        double dis = sqrt((xyz[0]-x)*(xyz[0]-x)+(xyz[1]-y)*(xyz[1]-y)+(xyz[2]-z)*(xyz[2]-z));
        if(args.getsourcetype()==0){
            double hmin=args.getsizing_value();
            double beta=args.getbeta();
            size = hmin+(beta-1)*dis;
        }
        else if(args.getsourcetype()==1){
            if(args.getspacingtype()==0){
                double spacing=args.getbeginSpacing();
                double decay=args.getbeginDecay();
                size = spacing+(-0.2*decay*decay*decay-0.1*decay*decay-0.2*decay+0.5)*dis;
            }
        }
        return size;
    }

    double lineSource::getsize(const double& x, const double& y, const double& z) {
        double dis;
        Eigen::Vector3d p0(xyz[0][0],xyz[0][1],xyz[0][2]);
        Eigen::Vector3d p1(xyz[1][0],xyz[1][1],xyz[1][2]);
        Eigen::Vector3d p(x,y,z);
        Eigen::Vector3d s=p1-p0,s0=p-p0,s1=p-p1;
        if((s.dot(s1)*s.dot(s0))<=0){
            dis = (s.cross(s1)).norm()/s.norm();
        }
        else{
            dis=std::min(s1.norm(),s0.norm());
        }
        double size;
        if(args.getsourcetype()==0){
            double hmin=args.getsizing_value();
            double beta=args.getbeta();
            size = hmin+(beta-1)*dis;
        }
        else if(args.getsourcetype()==1){
            if(args.getspacingtype()==0){
                double spacing=args.getbeginSpacing();
                double decay=args.getbeginDecay();
                size = spacing+(-0.2*decay*decay*decay-0.1*decay*decay-0.2*decay+0.5)*dis;
            }
            // else if(){
                
            // }
        }
        return size;
    }

    double triangleSource::getsize(const double& x, const double& y, const double& z) {
        double dis;
        double x0 = xyz[0][0], y0 = xyz[0][1], z0 = xyz[0][2];
        double x1 = xyz[1][0], y1 = xyz[1][1], z1 = xyz[1][2];
        double x2 = xyz[2][0], y2 = xyz[2][1], z2 = xyz[2][2];
        Eigen::Vector3d v0p(x-x0,y-y0,z-z0),v1p(x-x1,y-y1,z-z1),v2p(x-x2,y-y2,z-z2);
        Eigen::Vector3d v0v1(x1-x0,y1-y0,z1-z0),v1v0(x0-x1,y0-y1,z0-z1),
                        v1v2(x2-x1,y2-y1,z2-z1),v2v1(x1-x2,y1-y2,z1-z2),
                        v2v0(x0-x2,y0-y2,z0-z2),v0v2(x2-x0,y2-y0,z2-z0);
        Eigen::Vector3d n = (v0v1.cross(v0v2)).normalized();
        Eigen::Vector3d v0p_ = v0p - n * (v0p.dot(n));
        Eigen::Vector3d v1p_ = v1p - n * (v1p.dot(n));
        Eigen::Vector3d v2p_ = v2p - n * (v2p.dot(n));
        double tmp0 = (v1v2.cross(v1p_)).dot(n),
                tmp1 = (v2v0.cross(v2p_)).dot(n),
                tmp2 = (v0v1.cross(v0p_)).dot(n);

        if(tmp0>=0 && tmp1>=0 && tmp2>=0){
            dis = fabs(v0p.dot(n));          
        }
        else if(tmp0>=0 && tmp1<=0 && tmp2<=0){
            dis = v0p.norm();
        }
        else if(tmp0<=0 && tmp1>=0 && tmp2<=0){
            dis = v1p.norm();
        }
        else if(tmp0<=0 && tmp1<=0 && tmp2>=0){
            dis = v2p.norm();
        }
        else if(tmp0<=0 && tmp1>=0 && tmp2>=0){
            dis = (v1v2.cross(v1p)).norm()/v1v2.norm();
        }
        else if(tmp0>=0 && tmp1<=0 && tmp2>=0){
            dis = (v2v0.cross(v2p)).norm()/v2v0.norm();
        }
        else if(tmp0>=0 && tmp1>=0 && tmp2<=0){
            dis = (v0v1.cross(v0p)).norm()/v0v1.norm();
        }
        double size;
        if(args.getsourcetype()==0){
            double hmin=args.getsizing_value();
            double beta=args.getbeta();
            size = hmin+(beta-1)*dis;
        }
        else if(args.getsourcetype()==1){
            if(args.getspacingtype()==0){
                double spacing=args.getbeginSpacing();
                double decay=args.getbeginDecay();
                size = spacing+(-0.2*decay*decay*decay-0.1*decay*decay-0.2*decay+0.5)*dis;
            }
            // else if(){
                
            // }
        }
        return size;
    }

    double cubicSource::getsize(const double& x, const double& y, const double& z) {
        double dis;
        double xc=0,yc=0,zc=0;
        double x0=0,y0=0,z0=0;
        double x1=0,y1=0,z1=0;
        double x2=0,y2=0,z2=0;
        for(int i=0;i<8;i++){
            xc+=xyz[i][0];
            yc+=xyz[i][1];
            zc+=xyz[i][2];
            if(i==0||i==1||i==2||i==3){
                x0+=xyz[i][0];
                y0+=xyz[i][1];
                z0+=xyz[i][2];
            }
            if(i==0||i==1||i==4||i==5){
                x1+=xyz[i][0];
                y1+=xyz[i][1];
                z1+=xyz[i][2];
            }
            if(i==0||i==3||i==4||i==7){
                x2+=xyz[i][0];
                y2+=xyz[i][1];
                z2+=xyz[i][2];
            }
        }
        xc/=8;yc/=8;zc/=8;
        x0/=4;y0/=4;z0/=4;
        x1/=4;y1/=4;z1/=4;
        x2/=4;y2/=4;z2/=4;
        Eigen::Vector3d c0(x0-xc,y0-yc,z0-zc), c1(x1-xc,y1-yc,z1-zc), c2(x2-xc,y2-yc,z2-zc);
        Eigen::Vector3d p(x-xc,y-yc,z-zc);
        double d0=c0.norm(), d1=c1.norm(), d2=c2.norm();
        double p0=fabs(p.dot(c0)/d0), p1=fabs(p.dot(c1)/d1), p2=fabs(p.dot(c2)/d2);
        if( p0<=d0 && p1<=d1 && p2<=d2 ){
            dis = 0;
        }
        else if( p0>=d0 && p1<=d1 && p2<=d2 ){
            dis = p0-d0;
        }
        else if( p0<=d0 && p1>=d1 && p2<=d2 ){
            dis = p1-d1;
        }
        else if( p0<=d0 && p1<=d1 && p2>=d2 ){
            dis = p2-d2;
        }
        else if( p0<=d0 && p1>=d1 && p2>=d2 ){
            dis = sqrt((p1-d1)*(p1-d1)+(p2-d2)*(p2-d2));
        }
        else if( p0>=d0 && p1<=d1 && p2>=d2 ){
            dis = sqrt((p2-d2)*(p2-d2)+(p0-d0)*(p0-d0));
        }
        else if( p0>=d0 && p1>=d1 && p2<=d2 ){
            dis = sqrt((p0-d0)*(p0-d0)+(p1-d1)*(p1-d1));
        }
        else{
            dis = sqrt((p0-d0)*(p0-d0)+(p1-d1)*(p1-d1)+(p2-d2)*(p2-d2));
        }
        double size;
        if(args.getsourcetype()==0){
            double hmin=args.getsizing_value();
            double beta=args.getbeta();
            size = hmin+(beta-1)*dis;
        }
        else if(args.getsourcetype()==1){
            if(args.getspacingtype()==0){
                double spacing=args.getbeginSpacing();
                double decay=args.getbeginDecay();
                size = spacing+(-0.2*decay*decay*decay-0.1*decay*decay-0.2*decay+0.5)*dis;
            }
            // else if(){
                
            // }
        }
        return size;
    }


    double cylinderSource::getsize(const double& x, const double& y, const double& z){
        double dis;
        double x0 = xyz[0][0], y0 = xyz[0][1], z0 = xyz[0][2], r0 = radius[0];
        double x1 = xyz[1][0], y1 = xyz[1][1], z1 = xyz[1][2], r1 = radius[1];
        Eigen::Vector3d N01(x1-x0,y1-y0,z1-z0), v0p(x-x0,y-y0,z-z0);
        Eigen::Vector3d n01 = N01.normalized();
        double H = N01.norm();
        double r = (v0p.cross(n01)).norm(), h = v0p.dot(n01);
        Eigen::Vector2d r0r1(r1-r0,H), r1r0(r0-r1,-H), r0p(r-r0,h), r1p(r-r1,h-H);
        if( h<=0 && r<=r0 ){
            dis = -h;
        }
        else if( h>=H && r<=r1 ){
            dis = h-H;
        }
        else if( h>=0 && h<=H && r<=(r1-r0)*h/H+r0){
            dis = 0;
        }
        else if( r>=r0 && r1r0.dot(r0p)>=0){
            dis = r0p.norm();
        }
        else if( r>=r1 && r0r1.dot(r1p)>=0){
            dis = r1p.norm();
        }
        else{
            dis = fabs((r1-r0)*h-(r-r0)*H)/sqrt((r1-r0)*(r1-r0)+H*H);
        }

        double size;
        if(args.getsourcetype()==0){
            double hmin=args.getsizing_value();
            double beta=args.getbeta();
            size = hmin+(beta-1)*dis;
        }
        else if(args.getsourcetype()==1){
            if(args.getspacingtype()==0){         //constant
                double spacing=args.getbeginSpacing();
                double decay=args.getbeginDecay();
                if (dis != 0){
                    size = spacing+(-0.2*decay*decay*decay-0.1*decay*decay-0.2*decay+0.5)*dis;
                }
                else{
                    size = spacing;
                }
            }
            else if(args.getspacingtype()==1){    //Parametric
                double b_spacing=args.getbeginSpacing();
                double b_decay=args.getbeginDecay();
                double e_spacing=args.getendSpacing();
                double e_decay=args.getendDecay();
                if (dis != 0){
                    // size = spacing+(-0.2*decay*decay*decay-0.1*decay*decay-0.2*decay+0.5)*dis;
                }
                else{
                    double cur_dis = h;
                    double max_dis = H;
                    size = b_spacing + (e_spacing - b_spacing) * cur_dis / max_dis;
                }
            }
            else if(args.getspacingtype()==2){    //axis to perimeter
                double b_spacing=args.getbeginSpacing();
                double b_decay=args.getbeginDecay();
                double e_spacing=args.getendSpacing();
                double e_decay=args.getendDecay();
                if (dis != 0){
                    // size = spacing+(-0.2*decay*decay*decay-0.1*decay*decay-0.2*decay+0.5)*dis;
                }
                else{
                    double cur_dis = r;
                    double max_dis = r0;
                    size = b_spacing + (e_spacing - b_spacing) * cur_dis / max_dis;
                }
            }
            else if(args.getspacingtype()==3){    //center to perimeter
                double b_spacing=args.getbeginSpacing();
                double b_decay=args.getbeginDecay();
                double e_spacing=args.getendSpacing();
                double e_decay=args.getendDecay();
                if (dis != 0){
                    size = e_spacing+(-0.2*e_decay*e_decay*e_decay-0.1*e_decay*e_decay-0.2*e_decay+0.5)*dis;
                }
                else{
                    // size = spacing;
                }
            }
        }
        return size;

    }

    double sphereSource::getsize(const double& x, const double& y, const double& z){
        double dis;
        double x0 = xyz[0], y0 = xyz[1], z0 = xyz[2], r0 = radius[0];
        dis = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))-r0;

        double size;
        if(args.getsourcetype()==0){
            double hmin=args.getsizing_value();
            double beta=args.getbeta();
            size = hmin+(beta-1)*dis;
        }
        else if(args.getsourcetype()==1){
            if(args.getspacingtype()==0){
                double spacing=args.getbeginSpacing();
                double decay=args.getbeginDecay();
                size = spacing+(-0.2*decay*decay*decay-0.1*decay*decay-0.2*decay+0.5)*dis;
            }
            // else if(){
                
            // }
        }
        return size;

    }

    namespace size_field{
        double sizefunccovert(std::array<double,3> xyz){
            double size = API_Sizing_Query(xyz[0],xyz[1],xyz[2]);
            return size;
        }

        double SizingFunction_getSizeAtPoint(const double& x, const double& y, const double& z, SizingManager& manager){
            double size = std::numeric_limits<double>::max();
            if(manager.combinetype==0){
                for(int i=0;i<manager.sf.size();i++){                    //min value
                    double tmpsize = (*(manager.sf[i])).getsize(x,y,z);
                    if(tmpsize<size) size=tmpsize;
                }
            }
            // else if(manager.combinetype==1){    //Minimum distance   
            //     for(int i=0;i<manager.sf.size();i++){
            //         int obj;
            //         int max = 99999;
            //         double dis = (*(manager.sf[i])).getdistance(x,y,z);
            //         if(dis < max) {obj = i;}
            //     }
            //     size = (*(manager.sf[obj])).getsize(x,y,z);
            // }
            // else if(manager.combinetype==2){    //Inverse distance   
            //     double size = 0;
            //     for(int i=0;i<manager.sf.size();i++){
            //         double dis = (*(manager.sf[i])).getdistance(x,y,z);
            //         size += (1 / dis)*(*(manager.sf[i])).getsize(x,y,z);
            //     }
            // }
            // else if(manager.combinetype==3){    //Blend distance   
                    // weight = f(blend(mindis,inversedis))
            // }

            return size;
        }

        void setSizeField(const SurfaceMesh& SurfaceMesh,const BackGroundParameters& args,SizingManager& manager){
            int cvNums=0; int* cvPtNums=nullptr; double* cvPts=nullptr;
            int fcNums=0; int* fcPtNums=nullptr; double* fcPts=nullptr;
            int* lpCvNum=nullptr; int* lpCvs=nullptr;
            int bndPtNum=SurfaceMesh.coords.size();
            double* bndPts;
            bndPts = new double[bndPtNum*3];
            for(int i=0;i<bndPtNum;i++){
                bndPts[3*i] = SurfaceMesh.coords[i][0];
                bndPts[3*i+1] = SurfaceMesh.coords[i][1];
                bndPts[3*i+2] = SurfaceMesh.coords[i][2];
            }
            int bndFctNum=SurfaceMesh.tris.size();
            int* bndFcts; int* fctFace;
            bndFcts = new int[bndFctNum*3];
            fctFace = new int[bndFctNum];
            // std::cout << "region's size:" << SurfaceMesh.regions.size() <<std::endl;
            for(int i=0;i<bndFctNum;i++){
                bndFcts[3*i] = SurfaceMesh.tris[i][0];
                bndFcts[3*i+1] = SurfaceMesh.tris[i][1];
                bndFcts[3*i+2] = SurfaceMesh.tris[i][2];
                fctFace[i] = SurfaceMesh.regions[i];
            }
            int bdNum=0; int* bdFcNums=nullptr; int* bdFcs=nullptr;
            double spGlSettings[8];
            spGlSettings[0]=args.gethmax();
            spGlSettings[1]=args.gethmin();
            spGlSettings[2]=args.getbeta();
            spGlSettings[3]=args.gettheta();
            spGlSettings[4]=args.getproximity();
            spGlSettings[5]=args.getproximity();
            spGlSettings[6]=args.getproximity();
            spGlSettings[7]=-1;

            int spSettingNum[4] = {8, 0, 0, 0};
            double *spPtSettings = nullptr;
            double *spCvSettings = nullptr;
            double *spFcSettings = nullptr;
            double volmGrowthRatio = 1.2;
            int sfObjID;
            if(spGlSettings[0]<0||spGlSettings[1]<0){
                double tmpx = bndPts[0], tmpy = bndPts[1], tmpz = bndPts[2];
                double xmax = tmpx, ymax=tmpy, zmax = tmpz;
                double xmin = tmpx, ymin = tmpy, zmin = tmpz;
                for(int i=1; i<bndPtNum; i++){
                        int base =3*i;
                        double xx = bndPts[base];
                        double yy = bndPts[base+1];
                        double zz = bndPts[base+2];
                        xmax = fmax(xmax, xx);
                        xmin = fmin(xmin, xx);
                        ymax = fmax(ymax, yy);
                        ymin = fmin(ymin, yy);
                        zmax = fmax(zmax, zz);
                        zmin = fmin(zmin, zz);
                }
                double L = fmax(xmax-xmin, fmax((ymax-ymin), zmax-zmin));
                if(spGlSettings[0]<0) spGlSettings[0] = L/10;
                if(spGlSettings[1]<0) spGlSettings[1] = L*4/100000;
            }

            API_Create_SurfBKG_SF(
                cvNums,cvPtNums, cvPts,
                fcNums, fcPtNums, fcPts,
                lpCvNum, lpCvs,
                bndPtNum, bndPts,
                bndFctNum,  bndFcts, fctFace,
                bdNum, bdFcNums, bdFcs,
                spSettingNum,
                spGlSettings, spPtSettings,
                spCvSettings, spFcSettings,
                volmGrowthRatio,
                &sfObjID
            );

            // SizingFunction func = [&(mapBGMesh[sfObjID])](const double& x, const double& y, const double& z) { return mapBGMesh[sfObjID].getsize(x, y, z); };
            TiGER::backgroundsizefield bkgsf;
            std::cout << sfObjID << std::endl;
            std::cout << SIZING_FUNCTION::mapBGMesh[sfObjID].V.rows() << std::endl;
            std::cout << SIZING_FUNCTION::mapBGMesh[sfObjID].F.rows() << std::endl;
            std::cout << SIZING_FUNCTION::mapBGMesh[sfObjID].S.rows() << std::endl;

            bkgsf.backgroundmesh.coords.resize(SIZING_FUNCTION::mapBGMesh[sfObjID].V.rows());
            Eigen::MatrixXd mat_v(SIZING_FUNCTION::mapBGMesh[sfObjID].V.rows(), 3);
            for(int i=0;i<SIZING_FUNCTION::mapBGMesh[sfObjID].V.rows();i++){
                mat_v.row(i) << SIZING_FUNCTION::mapBGMesh[sfObjID].V(i,0), 
                    SIZING_FUNCTION::mapBGMesh[sfObjID].V(i,1), 
                    SIZING_FUNCTION::mapBGMesh[sfObjID].V(i,2);
                for(int j=0;j<3;j++){
                    bkgsf.backgroundmesh.coords[i][j]=SIZING_FUNCTION::mapBGMesh[sfObjID].V(i,j);
                    
                }
            }

            bkgsf.backgroundmesh.tris.resize(SIZING_FUNCTION::mapBGMesh[sfObjID].F.rows());
            Eigen::MatrixXi mat_f(SIZING_FUNCTION::mapBGMesh[sfObjID].F.rows(), 3);
            for(int i=0;i<SIZING_FUNCTION::mapBGMesh[sfObjID].F.rows();i++){
                mat_f.row(i) << SIZING_FUNCTION::mapBGMesh[sfObjID].F(i,0), 
                    SIZING_FUNCTION::mapBGMesh[sfObjID].F(i,1), 
                    SIZING_FUNCTION::mapBGMesh[sfObjID].F(i,2);
                for(int j=0;j<3;j++){
                    bkgsf.backgroundmesh.tris[i][j]=SIZING_FUNCTION::mapBGMesh[sfObjID].F(i,j);
                }
            }

            bkgsf.backgroundmesh.point_attribute_double.resize(SIZING_FUNCTION::mapBGMesh[sfObjID].S.rows());
            for(int i=0;i<SIZING_FUNCTION::mapBGMesh[sfObjID].S.rows();i++){
                bkgsf.backgroundmesh.point_attribute_double[i]=SIZING_FUNCTION::mapBGMesh[sfObjID].S(i,0);
            }

            bkgsf.backgroundmesh.regions.resize(SIZING_FUNCTION::mapBGMesh[sfObjID].tris_geoface.rows());
            for(int i=0;i<SIZING_FUNCTION::mapBGMesh[sfObjID].tris_geoface.rows();i++){
                bkgsf.backgroundmesh.regions[i]=SIZING_FUNCTION::mapBGMesh[sfObjID].tris_geoface(i,0);
            }

            // bkgsf.tree.init(SIZING_FUNCTION::mapBGMesh[sfObjID].V,SIZING_FUNCTION::mapBGMesh[sfObjID].F);
            // bkgsf.backgroundmesh.V=SIZING_FUNCTION::mapBGMesh[sfObjID].V;
            // bkgsf.backgroundmesh.F=SIZING_FUNCTION::mapBGMesh[sfObjID].F;
            bkgsf.octree_ = TiGER::common::BinaryTree<Eigen::Matrix3d>(mat_v, mat_f);
               //tree init
            

            std::shared_ptr<backgroundsizefield> ptr_to_bkgsf = std::make_shared<backgroundsizefield>(bkgsf);

            std::shared_ptr<sizeobj> ptr_to_sizeobj = std::static_pointer_cast<sizeobj>(ptr_to_bkgsf);

            manager.sf.emplace_back(ptr_to_sizeobj);
        }

        void SizingFunction_setUniformSize(double size,SizingManager& manager){
            TiGER::constsizefield constsf;
            constsf.size=size;
            std::shared_ptr<constsizefield> ptr_to_sf = std::make_shared<constsizefield>(constsf);
            std::shared_ptr<sizeobj> ptr_to_sizeobj = std::static_pointer_cast<sizeobj>(ptr_to_sf);
            manager.sf.emplace_back(ptr_to_sizeobj);
        }

        void SizingFunction_addPointSource(const std::array<double,3>& xyz,const PointSourceParameters& args,SizingManager& manager){
            TiGER::pointSource point_s;
            point_s.args=args;
            point_s.xyz=xyz;
            std::shared_ptr<pointSource> ptr_to_sf = std::make_shared<pointSource>(point_s);
            std::shared_ptr<sizeobj> ptr_to_sizeobj = std::static_pointer_cast<sizeobj>(ptr_to_sf);
            manager.sf.emplace_back(ptr_to_sizeobj);
        }

        void SizingFunction_addLineSource(const std::array<std::array<double,3>,2>& xyz,const LineSourceParameters& args,SizingManager& manager){
            TiGER::lineSource line_s;
            line_s.args=args;
            line_s.xyz=xyz;
            std::shared_ptr<lineSource> ptr_to_sf = std::make_shared<lineSource>(line_s);
            std::shared_ptr<sizeobj> ptr_to_sizeobj = std::static_pointer_cast<sizeobj>(ptr_to_sf);
            manager.sf.emplace_back(ptr_to_sizeobj);
        }

        void SizingFunction_addTriangleSource(const std::array<std::array<double,3>,3>& xyz,const TriangleSourceParameters& args,SizingManager& manager){
            TiGER::triangleSource triangle_s;
            triangle_s.args=args;
            triangle_s.xyz=xyz;
            std::shared_ptr<triangleSource> ptr_to_sf = std::make_shared<triangleSource>(triangle_s);
            std::shared_ptr<sizeobj> ptr_to_sizeobj = std::static_pointer_cast<sizeobj>(ptr_to_sf);
            manager.sf.emplace_back(ptr_to_sizeobj);
        }

        void SizingFunction_addCubicSource(const std::array<std::array<double,3>,8>& xyz,const CubicSourceParameters& args,SizingManager& manager){
            TiGER::cubicSource cubic_s;
            cubic_s.args=args;
            cubic_s.xyz=xyz;
            std::shared_ptr<cubicSource> ptr_to_sf = std::make_shared<cubicSource>(cubic_s);
            std::shared_ptr<sizeobj> ptr_to_sizeobj = std::static_pointer_cast<sizeobj>(ptr_to_sf);
            manager.sf.emplace_back(ptr_to_sizeobj);
        }

        void addCylinderSource(const std::array<std::array<double,3>,2>& xyz, const std::array<double,2>& radius, CylinderSourceParameters& args, SizingManager& manager){
            TiGER::cylinderSource cylinder_s;
            cylinder_s.args=args;
            cylinder_s.xyz=xyz;
            cylinder_s.radius=radius;
            std::shared_ptr<cylinderSource> ptr_to_sf = std::make_shared<cylinderSource>(cylinder_s);
            std::shared_ptr<sizeobj> ptr_to_sizeobj = std::static_pointer_cast<sizeobj>(ptr_to_sf);
            manager.sf.emplace_back(ptr_to_sizeobj);
        }

        void addSphereSource(const std::array<double,3>& xyz, const std::array<double,1>& radius, SphereSourceParameters& args, SizingManager& manager){
            TiGER::sphereSource sphere_s;
            sphere_s.args=args;
            sphere_s.xyz=xyz;
            sphere_s.radius=radius;
            std::shared_ptr<sphereSource> ptr_to_sf = std::make_shared<sphereSource>(sphere_s);
            std::shared_ptr<sizeobj> ptr_to_sizeobj = std::static_pointer_cast<sizeobj>(ptr_to_sf);
            manager.sf.emplace_back(ptr_to_sizeobj);
        }

        void setAttributeSize(const std::string& attribute_name,const std::vector<int>& attribute_index,
                              const double& sizing_value,backgroundsizefield& sf){
            if(attribute_name=="point"){
                sf.backgroundmesh.point_attribute_double[attribute_index[0]]=sizing_value;
            }
            else if(attribute_name=="line"){
                std::unordered_map<int,int> countmap;
                for(int i=0;i<sf.backgroundmesh.tris.size();i++){
                    if(sf.backgroundmesh.regions[i]==attribute_index[0]||sf.backgroundmesh.regions[i]==attribute_index[1]){
                        int v0=sf.backgroundmesh.tris[i][0];
                        int v1=sf.backgroundmesh.tris[i][1];
                        int v2=sf.backgroundmesh.tris[i][2];
                        countmap[v0]++;
                        countmap[v1]++;
                        countmap[v2]++;
                    }
                }
                for(auto &it:countmap){
                    if(it.second==2){
                        sf.backgroundmesh.point_attribute_double[it.first]=sizing_value;
                        // std::cout << "Get it" << it.first <<std::endl;
                    }
                }  
            }
            else if(attribute_name=="face"){
                for(int i=0;i<sf.backgroundmesh.tris.size();i++){
                    if(sf.backgroundmesh.regions[i]==attribute_index[0]){
                        int v0=sf.backgroundmesh.tris[i][0];
                        int v1=sf.backgroundmesh.tris[i][1];
                        int v2=sf.backgroundmesh.tris[i][2];
                        sf.backgroundmesh.point_attribute_double[v0]=sizing_value;
                        sf.backgroundmesh.point_attribute_double[v1]=sizing_value;
                        sf.backgroundmesh.point_attribute_double[v2]=sizing_value;
                    }
                }
            }
        }


    //     void SizingFunction_setSizeFromVolume(const VolumeMesh& VolumeMesh,SizingManager& manager,double scale){
    //         Eigen::MatrixXd V;
    //         Eigen::MatrixXi E;
    //         int sfObjID;
    //         V.resize(VolumeMesh.coords.size(),3);
    //         for(int i=0;i<VolumeMesh.coords.size();i++){
    //             for(int j=0;j<3;j++){
    //                 V(i,j)=VolumeMesh.coords[i][j];
    //             }
    //         }
    //         E.resize(VolumeMesh.tetras.size(),4);
    //         for(int i=0;i<VolumeMesh.tetras.size();i++){
    //             for(int j=0;j<4;j++){
    //                 E(i,j)=VolumeMesh.tetras[i][j];
    //             }
    //         }
    //         API_Create_SF_with_Volume_Mesh(V,E,&sfObjID,scale);
    //         SizingFunction func = API_Sizing_Query;
    //         manager.sf.emplace_back(func);
    //     }
    // }


    // void SizingFunction_addPointSource(const std::array<double, 3>& xyz,PointSourceParameters args,SizingManager& manager){
    //     double coord[3]={xyz[0],xyz[1],xyz[2]};
    //     double hmin=args.setendSpacing();
    //     double size=args.getbeta();
    //     API_Create_Point_Source(coord, hmin, beta);
    //     SizingFunction func = API_Sizing_Query; 
    //     manager.sf.emplace_back(func);
    // }


    
    }

}


#else
namespace TiGER {

    class SizingFieldImpl : public SizingField, public Source {
    public:
        // �޸Ĺ��캯�������� SizingField �� Source ������Ϊ����
        SizingFieldImpl(const SizingField& sizingField, const Source& source) : SizingField(sizingField), Source(source) {}

        const std::array<double, 9> sizeQuery(const std::array<double, 3>& xyz) const override {
            // ���� getSize �� Source ���еĺ�����ȷ����Ϊ const ��Ա����
            double size = getsize(xyz[0], xyz[1], xyz[2]);
            return {size, size, size, size, size, size, size, size, size};
        }
    };

    void SizingFunction_addPointSource(const std::array<double,3>& xyz,
                        const double& inner_radius,
                        const double& outter_radius,
                        const double& sizing_value,
                        const double& growth_value,
                        SizingManager& manager){
        Source pointSource;
        SizingField sizingField;
        pointSource.init(xyz[0], xyz[1], xyz[2], sizing_value, growth_value, inner_radius, outter_radius);
        SizingFieldImpl sizingFieldImpl(sizingField, pointSource);
        manager.sf.emplace_back(sizingField);
    }

    void SizingFunction_addLineSource(const std::array<std::array<double,3>,2>& xyz,
                        const double& inner_radius,
                        const double& outter_radius,
                        const double& sizing_value,
                        const double& growth_value,
                        SizingManager& manager){
        Source lineSource;
        SizingField sizingField;
        lineSource.init(xyz[0][0], xyz[0][1], xyz[0][2],
                        xyz[1][0], xyz[1][1], xyz[1][2],
                        sizing_value, growth_value, inner_radius, outter_radius);
        SizingFieldImpl sizingFieldImpl(sizingField, lineSource);
        manager.sf.emplace_back(sizingField);
    }

    void SizingFunction_addTriangleSource( const std::array<std::array<double,3>,3>& xyz,
                            const double& inner_radius,
                            const double& outter_radius,
                            const double& sizing_value,
                            const double& growth_value,
                            SizingManager& manager){
        Source triangleSource;
        SizingField sizingField;
        triangleSource.init(xyz[0][0], xyz[0][1], xyz[0][2],
                            xyz[1][0], xyz[1][1], xyz[1][2],
                            xyz[2][0], xyz[2][1], xyz[2][2],
                            sizing_value, growth_value, inner_radius, outter_radius);
        SizingFieldImpl sizingFieldImpl(sizingField, triangleSource);
        manager.sf.emplace_back(sizingField);
    }

    void SizingFunction_addFrustumSource(  const std::array<double,3>& upper_xyz,
                            const std::array<double,3>& lower_xyz,
                            const double& upper_radius,
                            const double& lower_radius,
                            const double& outter_radius,
                            const double& sizing_value,
                            const double& growth_value,
                            SizingManager& manager){
        Source frustumSource;
        SizingField sizingField;
        frustumSource.init_frustum_of_a_cone(   upper_xyz[0], upper_xyz[1], upper_xyz[2], upper_radius, 
                                                lower_xyz[0], lower_xyz[1], lower_xyz[2], lower_radius, 
                                                sizing_value, growth_value, outter_radius);
        SizingFieldImpl sizingFieldImpl(sizingField, frustumSource);
        manager.sf.emplace_back(sizingField);
    }

    void SizingFunction_addCubicSource(const std::array<double,3>& centerpoints,
                        const std::array<std::array<double,3>,3>& facecenterpoints,
                        const double& outter_radius,
                        const double& sizing_value,
                        const double& growth_value,
                        SizingManager& manager){
        Source cubicSource;
        SizingField sizingField;
        cubicSource.init_cuboid(centerpoints[0], centerpoints[1], centerpoints[2],
                                facecenterpoints[0][0], facecenterpoints[0][1], facecenterpoints[0][2],
                                facecenterpoints[1][0], facecenterpoints[1][1], facecenterpoints[1][2],
                                facecenterpoints[2][0], facecenterpoints[2][1], facecenterpoints[2][2],
                                sizing_value, growth_value, outter_radius);
        SizingFieldImpl sizingFieldImpl(sizingField, cubicSource);
        manager.sf.emplace_back(sizingField);
    }
    
    void SizingFunction_addPolygonSource(  const std::vector<std::array<double,3>>& xyz,
                            const double& inner_radius,
                            const double& outter_radius,
                            const double& sizing_value,
                            const double& growth_value,
                            SizingManager& manager){
        Source polygonSource;
        SizingField sizingField;
        polygonSource.init_polygon(xyz, sizing_value, growth_value, inner_radius, outter_radius);
        SizingFieldImpl sizingFieldImpl(sizingField, polygonSource);
        manager.sf.emplace_back(sizingField);
    }

    void SizingFunction_addEllipsoidSource(const std::array<double,3>& xyz1,
                            const std::array<double,3>& xyz2,
                            const double& a,
                            const double& outter_radius,
                            const double& sizing_value,
                            const double& growth_value,
                            SizingManager& manager){
        Source ellipsoidSource;
        SizingField sizingField;
        ellipsoidSource.init_ellipsoid( xyz1[0], xyz1[1], xyz1[2],
                                        xyz2[0], xyz2[1], xyz2[2],
                                        a, sizing_value, growth_value, outter_radius);
        SizingFieldImpl sizingFieldImpl(sizingField, ellipsoidSource);
        manager.sf.emplace_back(sizingField);
    }
}
#endif
