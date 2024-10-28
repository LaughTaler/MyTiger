#include"../include/extrude.h"
#include"../src/extrude_API.h"
#include "../include/extrude_parameters.h"
#include "alias.h"

namespace TiGER {
    namespace extrude {

        std::vector<std::vector<int>> find_adjacent(Eigen::MatrixXd V, Eigen::MatrixXi F) {
            std::vector<std::vector<int>> m;
            m.resize(V.rows());
            for (int j = 0; j < F.rows(); j++) {
                for (int k = 0; k < 3; k++) {
                    m[F(j, k)].push_back(F(j, (k + 1) % 3));
                }
            }
            return m;
        }

        std::vector<int> find_blpoint(Eigen::MatrixXd V, Eigen::MatrixXi F) {
            std::vector<int> x;
            x.resize(V.rows());
            std::map<pair<int, int>, int> y;
            for (int i = 0; i < F.rows(); i++) {
                int a = F(i, 0);
                int b = F(i, 1);
                int c = F(i, 2);
                if (a > b) {
                    int temp = b;
                    b = a;
                    a = temp;
                }
                if (a > c) {
                    int temp = c;
                    c = a;
                    a = temp;
                }
                if (b > c) {
                    int temp = c;
                    c = b;
                    b = temp;
                }
                y[std::pair<int, int>(a, b)]++;
                y[std::pair<int, int>(a, c)]++;
                y[std::pair<int, int>(b, c)]++;
            }
            for (auto it = y.begin(); it != y.end(); it++) {
                if (it->second == 1) {
                    x[it->first.first]++;
                    x[it->first.second]++;
                }
            }
            return x;
        }

        Eigen::MatrixXd normalsmooth(Eigen::MatrixXd N, std::vector<std::vector<int>> m) {
            Eigen::MatrixXd newN = N;
            for (int i = 0; i < N.rows(); i++) {
                int temp = m[i].size();
                Eigen::Vector3d sum(0, 0, 0);
                for (int j = 0; j < temp; j++) {
                    for (int k = 0; k < 3; k++) {
                        sum(k) += N(m[i][j], k);
                    }
                }
                for (int k = 0; k < 3; k++) {
                    sum(k) /= temp;
                    newN(i, k) = sum(k);
                }
            }
            return newN;
        }

        void eigen_sort_3d(const Eigen::MatrixXd& in_vec, Eigen::MatrixXd& out_vec, Eigen::VectorXi& ind, double eps)
        {
            ind = Eigen::VectorXi::LinSpaced(in_vec.rows(), 0, in_vec.rows() - 1);
            auto rule = [&in_vec, &eps](int i, int j)->bool {
                if (fabs(in_vec(i, 0) - in_vec(j, 0)) > eps) return in_vec(i, 0) < in_vec(j, 0);
                if (fabs(in_vec(i, 1) - in_vec(j, 1)) > eps) return in_vec(i, 1) < in_vec(j, 1);
                if (fabs(in_vec(i, 2) - in_vec(j, 2)) > eps) return in_vec(i, 2) < in_vec(j, 2);
                return i < j;
                };
            std::sort(ind.data(), ind.data() + ind.size(), rule);
            out_vec.resize(in_vec.rows(), in_vec.cols());
            for (int i = 0; i < in_vec.rows(); i++)
            {
                out_vec.row(i) << in_vec.row(ind(i));
            }
            return;
        }

        bool removeDulplicatePoint(Eigen::MatrixXd& V, Eigen::MatrixXi& F, double eps) {
            Eigen::MatrixXd V_in = V;
            Eigen::MatrixXi F_in = F;
            Eigen::MatrixXd V_sort;
            Eigen::VectorXi ind;
            eigen_sort_3d(V_in, V_sort, ind, eps);
            int cnt = 0;
            std::unordered_map<int, int> mpid;
            std::vector<int> new_point_lst;
            for (int i = 0; i < V_sort.rows(); i++)
            {
                while (i + 1 < V_sort.rows())
                {
                    if ((V_sort.row(i + 1) - V_sort.row(i)).norm() < eps)
                    {
                        mpid[ind[i]] = cnt;
                        i++;
                    }
                    else break;
                }

                mpid[ind[i]] = cnt;
                new_point_lst.push_back(ind[i]);
                cnt++;
            }

            V.resize(new_point_lst.size(), 3);
            F.resize(F_in.rows(), 3);

            for (int i = 0; i < new_point_lst.size(); i++)
            {
                V.row(i) = V_in.row(new_point_lst[i]);
            }

            for (int i = 0; i < F_in.rows(); i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    F(i, j) = mpid[F_in(i, j)];
                }
            }
            Eigen::VectorXi I;
            auto V_old = V;
            auto F_old = F;

            igl::remove_unreferenced(V_old, F_old, V, F, I);

            std::cout << "Remove " << V_in.rows() - V.rows() << " duplicate points\n";
            std::cout << "Remove " << F_in.rows() - F.rows() << " facets\n";
            return true;
        }

        bool volunme_judge(Eigen::MatrixXd NV, Eigen::MatrixXi NF) {
            bool b = 1;
            for (int i = 0; i < NF.rows(); i++) {
                double p1[3] = { NV(NF(i,0),0),NV(NF(i,0),1),NV(NF(i,0),2) };
                double p2[3] = { NV(NF(i,1),0),NV(NF(i,1),1),NV(NF(i,1),2) };
                double p3[3] = { NV(NF(i,2),0),NV(NF(i,2),1),NV(NF(i,2),2) };
                double p4[3] = { NV(NF(i,3),0),NV(NF(i,3),1),NV(NF(i,3),2) };
                double p5[3] = { NV(NF(i,4),0),NV(NF(i,4),1),NV(NF(i,4),2) };
                double p6[3] = { NV(NF(i,5),0),NV(NF(i,5),1),NV(NF(i,5),2) };
                double ori[12] = { 0 };
                ori[0] = TIGER_GEOM_FUNC::orient3d(p1, p2, p4, p3) / 6;
                ori[1] = TIGER_GEOM_FUNC::orient3d(p1, p2, p5, p3) / 6;
                ori[2] = TIGER_GEOM_FUNC::orient3d(p1, p2, p6, p3) / 6;
                ori[3] = TIGER_GEOM_FUNC::orient3d(p4, p5, p6, p1) / 6;
                ori[4] = TIGER_GEOM_FUNC::orient3d(p4, p5, p6, p2) / 6;
                ori[5] = TIGER_GEOM_FUNC::orient3d(p4, p5, p6, p3) / 6;

                ori[6] = TIGER_GEOM_FUNC::orient3d(p1, p4, p6, p2) / 6;
                ori[7] = TIGER_GEOM_FUNC::orient3d(p4, p6, p3, p2) / 6;
                ori[8] = TIGER_GEOM_FUNC::orient3d(p2, p5, p4, p3) / 6;
                ori[9] = TIGER_GEOM_FUNC::orient3d(p5, p4, p1, p3) / 6;
                ori[10] = TIGER_GEOM_FUNC::orient3d(p3, p6, p5, p1) / 6;
                ori[11] = TIGER_GEOM_FUNC::orient3d(p6, p5, p2, p1) / 6;
                for (int j = 0; j < 12; j++) {
                    if (ori[j] < 0)return b = 0;
                }
            }
            return b;
        }

        void initial_input(const TiGER::SurfaceMesh& surface_mesh, const TiGER::ExtrudeParameters& args, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& N,
            double& Layerthick, double& initial_length, double &height,short &LayerNum) {
            //VΪ������� FΪ�����α�� NΪ����
            
            auto gpp=args.getgeo_progre_para();
            auto gpp1 = args.getnum_layer();
            Layerthick = gpp.growth_rate;
            initial_length = gpp.initial_delta_s;
            height = initial_length;
            LayerNum = gpp1;
            F.resize(surface_mesh.tris.size(), 3);
            V.resize(surface_mesh.coord.size(), 3);
            N.resize(surface_mesh.coord.size(), 3);
            for (int i = 0; i < surface_mesh.tris.size(); i++) {
                for (int j = 0; j < 3; j++) {
                    F(i, j) = surface_mesh.tris[i][j];
                  
                }
            }
            for (int i = 0; i < surface_mesh.coord.size(); i++) {
                for (int j = 0; j < 3; j++) {
                    V(i, j) = surface_mesh.coord[i][j];
                    N(i, j) = surface_mesh.point_normal[i][j];
                }
            }
            return;
        }

        void extrude_realize(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd N, double Layerthick, double initial_length, double height, short LayerNum, Mesh &mesh) {
            TIGER_GEOM_FUNC::exactinit();
           // removeDulplicatePoint(V, F, 1e-9);
            //int LayerNum = 40;
            std::vector<std::vector<int>> m;
            Eigen::MatrixXi NF;
            Eigen::MatrixXd NV;
            NF.resize(F.rows() * LayerNum, 6);
            NV.resize(V.rows() * (LayerNum + 1), 3);
            NV.topRows(V.rows()) = V;
            std::vector<std::vector<int>> adpoint;
            adpoint = find_adjacent(V, F);
            Eigen::MatrixXd N1;
            N1 = normalsmooth(N, adpoint);
            Eigen::MatrixXd LastV;
            Eigen::MatrixXi OF = F;
            Eigen::MatrixXd tempN = N;
            for (int i = 0; i < LayerNum; i++) {
                //printf("%d %d\n", V.rows(), N1.rows());
                Eigen::MatrixXd OV = V + N1 * height;

                // igl::cat(1, V, OV, NV);
                for (int k = 0; k < V.rows(); k++) {
                    NV.row(V.rows() * (i + 1) + k) = OV.row(k);
                }
                for (int j = 0; j < F.rows(); j++) {
                    for (int k = 0; k < 3; k++) {
                        NF(j + F.rows() * i, k) = F(j, k) + i * V.rows();
                        NF(j + F.rows() * i, k + 3) = F(j, k) + (i + 1) * V.rows();
                        OF(j, k) = F(j, k) + (i + 1) * V.rows();
                    }
                }

                V = OV;
                height *= Layerthick;
                igl::per_vertex_normals(V, F, N);
                adpoint = find_adjacent(V, F);
                N1 = normalsmooth(N, adpoint);
                tempN = N1;
                /*for (int h = 0; h < 10; h++) {
                    N1 = normalsmooth(tempN, adpoint);
                    tempN = N1;
                }*/
                //�⻬������
            }
           
            mesh.Vertex = NV;
            mesh.Topo = NF;
            MESHIO::writeVTK("333.vtk", mesh);
            bool b;
            b = volunme_judge(NV, NF);
            return;
        }

        void change_output(TiGER::VolumeMesh& vol_mesh, Mesh mesh) {
            vol_mesh.coords.resize(mesh.Vertex.rows());
            vol_mesh.prisms.resize(mesh.Topo.rows());
            for (int i = 0; i < mesh.Vertex.rows(); i++) {
                for (int j = 0; j < 3; j++) {
                    vol_mesh.coords[i][j] = mesh.Vertex(i, j);
                }
            }
            for (int i = 0; i < mesh.Topo.rows(); i++) {
                for (int j = 0; j < 6; j++) {
                    vol_mesh.prisms[i][j] = mesh.Topo(i, j);
                }
            }
            return;
        }

        void readstltoSurfaceMesh(SurfaceMesh& surface_mesh, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd N) {
            surface_mesh.tris.resize(F.rows());
            surface_mesh.coord.resize(V.rows());
            surface_mesh.point_normal.resize(V.rows());
            for (int i = 0; i < F.rows(); i++) {
                for (int j = 0; j < 3; j++) {
                    surface_mesh.tris[i][j]= F(i, j);
                   
                }
            }
            for (int i = 0; i < V.rows(); i++) {
                for (int j = 0; j < 3; j++) {
                    surface_mesh.coord[i][j]= V(i, j);
                    surface_mesh.point_normal[i][j] = N(i, j);
                }
            }
            return;
        }

        /**
         * @brief ʹ��ǰ�ز�������ɶ��������
         *
         * @param line_mesh ������
         * @param args ����
         * @param surface_mesh �ı�������
         */
        void HexGrid_Extrude2D(const TiGER::LineMesh& line_mesh, const ExtrudeParameters& args, SurfaceMesh& surface_mesh) {
            return;
        }

        /**
         * @brief ʹ��ǰ�ز����������������
         *
         * @param surface_mesh �����Ρ��ı�������
         * @param args
         * @param vol_mesh ��������
         */
        void HexGrid_Extrude3D(const SurfaceMesh& surface_mesh, const ExtrudeParameters& args,VolumeMesh& vol_mesh) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            Eigen::MatrixXd N;
            double Layerthick ;
            double initial_length;
            double height;
            short LayerNum;
            Mesh mesh;
            initial_input(surface_mesh, args, V, F, N, Layerthick, initial_length, height,LayerNum);
            extrude_realize(V, F, N, Layerthick, initial_length, height,LayerNum,mesh);
            change_output(vol_mesh, mesh);
            return;
        }
        
        ;
    }
}


