#ifndef TiGER_COMMON_GEOMETRY_INTERSECTION_H_
#define TiGER_COMMON_GEOMETRY_INTERSECTION_H_

#include "Eigen/Dense"

namespace TiGER {
namespace common {
double direction(const Eigen::Vector2d& pi, const Eigen::Vector2d& pj,
                 const Eigen::Vector2d& pk);
/**
 * @brief:
 *        pj
 *      /
 *    /
 * pi----- pk
 * calculate the (pj-pi)x(pk-pi)
 */
double angleSin(const Eigen::Vector2d& pi, const Eigen::Vector2d& pj,
                const Eigen::Vector2d& pk);

/**
 * @brief:
 *        a
 *      /
 *    /
 * b ----- c
 * calculate the angle of (ba and bc)
 */
double angleDegree(const Eigen::RowVector3d& a, const Eigen::RowVector3d& b,
                   const Eigen::RowVector3d& c);

bool onSegment(const Eigen::Vector2d& pi, const Eigen::Vector2d& pj,
               const Eigen::Vector2d& pk);

bool checkIntersection(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                       const Eigen::Vector2d& c, const Eigen::Vector2d& d);

bool lineTriIntersection(const Eigen::RowVector3d& e1,
                         const Eigen::RowVector3d& e2,
                         const Eigen::RowVector3d& a,
                         const Eigen::RowVector3d& b,
                         const Eigen::RowVector3d& c);

bool TriTriIntersection(const std::array<Eigen::RowVector3d, 6>& T);

bool TriTriIntersection(const std::array<Eigen::RowVector3d, 3>& t1,
                        const std::array<Eigen::RowVector3d, 3>& t2);


bool TriTriIntersection(const Eigen::MatrixX3d& t1_tmp,
                        const Eigen::MatrixX3d& t2_tmp);

bool TriTriIntersection(const Eigen::RowVector3d& a1,
                        const Eigen::RowVector3d& b1,
                        const Eigen::RowVector3d& c1,
                        const Eigen::RowVector3d& a2,
                        const Eigen::RowVector3d& b2,
                        const Eigen::RowVector3d& c2);

bool lineTriIntersectionGeo(const Eigen::RowVector3d& e1,
                            const Eigen::RowVector3d& e2,
                            const Eigen::RowVector3d& a,
                            const Eigen::RowVector3d& b,
                            const Eigen::RowVector3d& c, const double d1,
                            const double d2);

bool TriTriIntersectionGeo(const Eigen::RowVector3d& a1,
                           const Eigen::RowVector3d& b1,
                           const Eigen::RowVector3d& c1,
                           const Eigen::RowVector3d& a2,
                           const Eigen::RowVector3d& b2,
                           const Eigen::RowVector3d& c2);

bool lineTri2DIntersection(const Eigen::Vector2d& e1, const Eigen::Vector2d& e2,
                           const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                           const Eigen::Vector2d& c);

bool LineLine2DIntersection(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                            const Eigen::Vector2d& c, const Eigen::Vector2d& d);

bool TriTriIntersection_tetgen(const Eigen::RowVector3d& a1,
                               const Eigen::RowVector3d& b1,
                               const Eigen::RowVector3d& c1,
                               const Eigen::RowVector3d& a2,
                               const Eigen::RowVector3d& b2,
                               const Eigen::RowVector3d& c2);

bool lineTriIntersection_tetgen(const Eigen::RowVector3d& e1,
                                const Eigen::RowVector3d& e2,
                                const Eigen::RowVector3d& a,
                                const Eigen::RowVector3d& b,
                                const Eigen::RowVector3d& c);

bool PointInTriangle(const Eigen::RowVector3d& p, const Eigen::RowVector3d& a,
                     const Eigen::RowVector3d& b, const Eigen::RowVector3d& c);

Eigen::RowVector3d closet_point_to_segment(const Eigen::RowVector3d& p,
                                           const Eigen::RowVector3d& s,
                                           const Eigen::RowVector3d& e);
}  // namespace common
}  // namespace TiGER
#include "geometry_intersection.cpp"
#endif  // !TiGER_COMMON_GEOMETRY_INTERSECTION_H_
