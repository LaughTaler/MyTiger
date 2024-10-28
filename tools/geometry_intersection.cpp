#define _USE_MATH_DEFINES
#include "geometry_intersection.h"

#include <cmath>
#include <iostream>
#include <vector>

// #include "tetgen.h"

#include "meshIO.h"
// #include "igl/predicates/predicates.h"
// #include "predicates.cpp"
// #include "geom_func.h"
#include <geom_func.h>
// #include "igl/tri_tri_intersect.h"

typedef double REAL;

namespace TiGER {
namespace common {
// using namespace predicates::adaptive;
bool lineTriIntersection(const Eigen::RowVector3d& e1,
                         const Eigen::RowVector3d& e2,
                         const Eigen::RowVector3d& a,
                         const Eigen::RowVector3d& b,
                         const Eigen::RowVector3d& c, double d1, double d2) {
  // if (std::fabs(d1) < 1e-10) d1 = 0;
  // if (std::fabs(d2) < 1e-10) d2 = 0;

  // std::cout<<d1<<" "<<d2<<std::endl;

  if (d1 * d2 > 0) return false;

  if (d1 == 0 && d2 != 0) {
    return PointInTriangle(e1, a, b, c);
  }

  if (d1 != 0 && d2 == 0) {
    return PointInTriangle(e2, a, b, c);
  }

  if ((d1 == 0 && d2 == 0)) {
    // return false;
    Eigen::RowVector3d abc_normal = (b - a).cross(c - a).normalized();
    std::vector<Eigen::RowVector3d> axis_normal = {Eigen::RowVector3d(1, 0, 0),
                                                   Eigen::RowVector3d(0, 1, 0),
                                                   Eigen::RowVector3d(0, 0, 1)};
    int axis = 0;
    double xx = 0;
    for (int i = 0; i < 3; ++i) {
      xx = abc_normal.dot(axis_normal[i]);
      if (fabs(abc_normal.dot(axis_normal[i])) > 0.5) {
        axis = i;
        break;
      }
    }

    int x = (axis + 1) % 3;
    int y = (axis + 2) % 3;
    Eigen::Vector2d e1_(e1[x], e1[y]);
    Eigen::Vector2d e2_(e2[x], e2[y]);
    Eigen::Vector2d a_(a[x], a[y]);
    Eigen::Vector2d b_(b[x], b[y]);
    Eigen::Vector2d c_(c[x], c[y]);

    return lineTri2DIntersection(e1_, e2_, a_, b_, c_);
  } else {
    REAL* real_e1 = const_cast<REAL*>(e1.data());
    REAL* real_e2 = const_cast<REAL*>(e2.data());
    REAL* real_a = const_cast<REAL*>(a.data());
    REAL* real_b = const_cast<REAL*>(b.data());
    REAL* real_c = const_cast<REAL*>(c.data());

    double a1 = TiGER_GEOM_FUNC::orient3d(real_e1, real_e2, real_a, real_b);
    double a2 = TiGER_GEOM_FUNC::orient3d(real_e1, real_e2, real_b, real_c);
    double a3 = TiGER_GEOM_FUNC::orient3d(real_e1, real_e2, real_c, real_a);

    // std::cout<<a1<<" "<<a2<<" "<<a3<<std::endl;

    // if (std::fabs(a1) < 1e-10) a1 = 0;
    // if (std::fabs(a2) < 1e-10) a2 = 0;
    // if (std::fabs(a3) < 1e-10) a3 = 0;

    // std::cout<<a1<<" "<<a2<<" "<<a3<<std::endl;

    if (a1 >= 0 && a2 >= 0 && a3 >= 0) return true;
    if (a1 <= 0 && a2 <= 0 && a3 <= 0) return true;

    return false;
  }
  return false;
}
// #pragma optimize("",off)
// bool lineTriIntersection_tetgen(const Eigen::RowVector3d& e1,
//                                 const Eigen::RowVector3d& e2,
//                                 const Eigen::RowVector3d& a,
//                                 const Eigen::RowVector3d& b,
//                                 const Eigen::RowVector3d& c) {
//   tetgenmesh t;
//   tetgenmesh::point pe1 = const_cast<double*>(e1.data());
//   tetgenmesh::point pe2 = const_cast<double*>(e2.data());
//   tetgenmesh::point pa = const_cast<double*>(a.data());
//   tetgenmesh::point pb = const_cast<double*>(b.data());
//   tetgenmesh::point pc = const_cast<double*>(c.data());
//   double d1 = TiGER_GEOM_FUNC::orient3d(e1.data(), a.data(), b.data(), c.data());
//   double d2 = TiGER_GEOM_FUNC::orient3d(e2.data(), a.data(), b.data(), c.data());
//   return t.tri_edge_inter_tail(pa, pb, pc, pe1, pe2, d1, d2);
// }

bool lineTriIntersection(const Eigen::RowVector3d& e1,
                         const Eigen::RowVector3d& e2,
                         const Eigen::RowVector3d& a,
                         const Eigen::RowVector3d& b,
                         const Eigen::RowVector3d& c) {
  if (e1 == a || e1 == b || e1 == c) return false;
  if (e2 == a || e2 == b || e2 == c) return false;
  REAL* real_e1 = const_cast<REAL*>(e1.data());
  REAL* real_e2 = const_cast<REAL*>(e2.data());
  REAL* real_a = const_cast<REAL*>(a.data());
  REAL* real_b = const_cast<REAL*>(b.data());
  REAL* real_c = const_cast<REAL*>(c.data());
  double d1 = TiGER_GEOM_FUNC::orient3d(real_e1, real_a, real_b, real_c);
  double d2 = TiGER_GEOM_FUNC::orient3d(real_e2, real_a, real_b, real_c);
  double t = (b - a).cross(c - a).dot(e2 - a);
  return lineTriIntersection(e1, e2, a, b, c, d1, d2);
}

bool lineTri2DIntersection(const Eigen::Vector2d& e1, const Eigen::Vector2d& e2,
                           const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                           const Eigen::Vector2d& c) {
  if (e1 == a || e1 == b || e1 == c) return false;
  if (e2 == a || e2 == b || e2 == c) return false;
  bool is_intersection = 0;
  is_intersection = LineLine2DIntersection(e1, e2, a, b);
  if (is_intersection) return true;
  is_intersection = LineLine2DIntersection(e1, e2, b, c);
  if (is_intersection) return true;
  is_intersection = LineLine2DIntersection(e1, e2, c, a);
  if (is_intersection) return true;

  REAL* real_e1 = const_cast<REAL*>(e1.data());
  REAL* real_e2 = const_cast<REAL*>(e2.data());
  REAL* real_a = const_cast<REAL*>(a.data());
  REAL* real_b = const_cast<REAL*>(b.data());
  REAL* real_c = const_cast<REAL*>(c.data());

  double d1 = TiGER_GEOM_FUNC::orient2d(real_e1, real_a, real_b);
  double d2 = TiGER_GEOM_FUNC::orient2d(real_e1, real_b, real_c);
  double d3 = TiGER_GEOM_FUNC::orient2d(real_e1, real_c, real_a);
  double d4 = TiGER_GEOM_FUNC::orient2d(real_e2, real_a, real_b);
  double d5 = TiGER_GEOM_FUNC::orient2d(real_e2, real_b, real_c);
  double d6 = TiGER_GEOM_FUNC::orient2d(real_e2, real_c, real_a);

  if ((d1 > 0 && d2 > 0 && d3 > 0 && d4 > 0 && d5 > 0 && d6 > 0) ||
      (d1 < 0 && d2 < 0 && d3 < 0 && d4 < 0 && d5 < 0 && d6 < 0)) {
    return true;
  }

  return false;
}

bool LineLine2DIntersection(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                            const Eigen::Vector2d& c,
                            const Eigen::Vector2d& d) {
  double eps = (a - b).norm();
  eps = std::min(eps, (c - d).norm());
  eps *= 1e-10;

  auto onSegment = [](const Eigen::Vector2d& pi, const Eigen::Vector2d& pj,
                      const Eigen::Vector2d& pk) {
    if ((pk - pi).norm() < 1e-10 || (pk - pj).norm() < 1e-10) return false;
    return ((std::min(pi.x(), pj.x()) <= pk.x() &&
             pk.x() <= std::max(pi.x(), pj.x()) &&
             std::min(pi.y(), pj.y()) <= pk.y() &&
             pk.y() <= std::max(pi.y(), pj.y())));
  };

  REAL* real_a = const_cast<REAL*>(a.data());
  REAL* real_b = const_cast<REAL*>(b.data());
  REAL* real_c = const_cast<REAL*>(c.data());
  REAL* real_d = const_cast<REAL*>(d.data());

  double d1 = TiGER_GEOM_FUNC::orient2d(real_c, real_d, real_a);
  double d2 = TiGER_GEOM_FUNC::orient2d(real_c, real_d, real_b);
  double d3 = TiGER_GEOM_FUNC::orient2d(real_a, real_b, real_c);
  double d4 = TiGER_GEOM_FUNC::orient2d(real_a, real_b, real_d);
  
  if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
      ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
    return true;
  }
  if (fabs(d1) < eps && onSegment(c, d, a)) {
    return true;
  }
  if (fabs(d2) < eps && onSegment(c, d, b)) {
    return true;
  }
  if (fabs(d3) < eps && onSegment(a, b, c)) {
    return true;
  }
  if (fabs(d4) < eps && onSegment(a, b, d)) {
    return true;
  }
  return false;
};

bool lineTriIntersectionGeo(const Eigen::RowVector3d& e1,
                            const Eigen::RowVector3d& e2,
                            const Eigen::RowVector3d& a,
                            const Eigen::RowVector3d& b,
                            const Eigen::RowVector3d& c, const double d1,
                            const double d2) {
  double eps = (e1 - e2).norm();
  eps = std::min(eps, (a - b).norm());
  eps = std::min(eps, (a - c).norm());
  eps = std::min(eps, (c - b).norm());
  eps *= 1e-3;
  if ((e1 - a).norm() < eps || (e1 - b).norm() < eps || (e1 - c).norm() < eps)
    return false;
  if ((e2 - a).norm() < eps || (e2 - b).norm() < eps || (e2 - c).norm() < eps)
    return false;
  // if (e1 == a || e1 == b || e1 == c) return false;
  // if (e2 == a || e2 == b || e2 == c) return false;
  return lineTriIntersection(e1, e2, a, b, c, d1, d2);
}
// #pragma optimize("", on)
bool TriTriIntersectionGeo(const Eigen::RowVector3d& a1,
                           const Eigen::RowVector3d& b1,
                           const Eigen::RowVector3d& c1,
                           const Eigen::RowVector3d& a2,
                           const Eigen::RowVector3d& b2,
                           const Eigen::RowVector3d& c2) {
  REAL* real_a1 = const_cast<REAL*>(a1.data());
  REAL* real_b1 = const_cast<REAL*>(b1.data());
  REAL* real_c1 = const_cast<REAL*>(c1.data());
  REAL* real_a2 = const_cast<REAL*>(a2.data());
  REAL* real_b2 = const_cast<REAL*>(b2.data());
  REAL* real_c2 = const_cast<REAL*>(c2.data());

  double da1 = TiGER_GEOM_FUNC::orient3d(real_a1, real_a2, real_b2, real_c2);
  double db1 = TiGER_GEOM_FUNC::orient3d(real_b1, real_a2, real_b2, real_c2);
  double dc1 = TiGER_GEOM_FUNC::orient3d(real_c1, real_a2, real_b2, real_c2);

  double da2 = TiGER_GEOM_FUNC::orient3d(real_a2, real_a1, real_b1, real_c1);
  double db2 = TiGER_GEOM_FUNC::orient3d(real_b2, real_a1, real_b1, real_c1);
  double dc2 = TiGER_GEOM_FUNC::orient3d(real_c2, real_a1, real_b1, real_c1);

  if (da1 * db1 <= 0) {
    return lineTriIntersectionGeo(a1, b1, a2, b2, c2, da1, db1);
  }

  if (dc1 * db1 <= 0) {
    return lineTriIntersectionGeo(c1, b1, a2, b2, c2, dc1, db1);
  }

  if (da1 * dc1 <= 0) {
    return lineTriIntersectionGeo(a1, c1, a2, b2, c2, da1, dc1);
  }

  if (da2 * db2 <= 0) {
    return lineTriIntersectionGeo(a2, b2, a1, b1, c1, da2, db2);
  }

  if (dc2 * db2 <= 0) {
    return lineTriIntersectionGeo(c2, b2, a1, b1, c1, dc2, db2);
  }

  if (da2 * dc2 <= 0) {
    return lineTriIntersectionGeo(a2, c2, a1, b1, c1, da2, dc2);
  }
  return false;
}

bool TriTriIntersection(const Eigen::MatrixX3d& t1_tmp,
                        const Eigen::MatrixX3d& t2_tmp) {
  std::array<Eigen::RowVector3d, 3> t1;
  std::array<Eigen::RowVector3d, 3> t2;
  for (int i = 0; i < 3; ++i) {
    t1[i] = t1_tmp.row(i);
    t2[i] = t2_tmp.row(i);
  }
  return TriTriIntersection(t1, t2);
}

bool TriTriIntersection(const std::array<Eigen::RowVector3d, 3>& t1,
                        const std::array<Eigen::RowVector3d, 3>& t2) {
  int coindex_count = 0;
  int coindex[2];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (t1[i] == t2[j]) {
        coindex[0] = i;
        coindex[1] = j;
        coindex_count++;
      }
    }
  }

  if (coindex_count == 0) {
    return TriTriIntersection(t1[0], t1[1], t1[2], t2[0], t2[1], t2[2]);
  } else if (coindex_count == 1) {
    return lineTriIntersection(t1[(coindex[0] + 1) % 3],
                               t1[(coindex[0] + 2) % 3], t2[0], t2[1], t2[2]) ||
           lineTriIntersection(t2[(coindex[1] + 1) % 3],
                               t2[(coindex[1] + 2) % 3], t1[0], t1[1], t1[2]);
  }
  return false;
}

bool TriTriIntersection(const std::array<Eigen::RowVector3d, 6>& T) {
  int coindex_count = 0;
  int coindex[2];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (T[i] == T[3 + j]) {
        coindex[0] = i;
        coindex[1] = j;
        coindex_count++;
      }
    }
  }

  if (coindex_count == 0) {
    return TriTriIntersection(T[0], T[1], T[2], T[3], T[4], T[5]);
    // bool coplane;
    // Eigen::Vector3d source, target;
    // return igl::tri_tri_intersection_test_3d(T[0], T[1], T[2], T[3], T[4], T[5],
    //                                          coplane, source, target);
  } else if (coindex_count == 1) {
    return lineTriIntersection(T[(coindex[0] + 1) % 3], T[(coindex[0] + 2) % 3],
                               T[3], T[4], T[5]) ||
           lineTriIntersection(T[(coindex[1] + 1) % 3 + 3],
                               T[(coindex[1] + 2) % 3 + 3], T[0], T[1], T[2]);
  }
  return false;
}
bool PointInTriangle(const Eigen::RowVector3d& p, const Eigen::RowVector3d& a,
                     const Eigen::RowVector3d& b, const Eigen::RowVector3d& c) {
  Eigen::RowVector3d abc_normal = (b - a).cross(c - a).normalized();
  std::vector<Eigen::RowVector3d> axis_normal = {Eigen::RowVector3d(1, 0, 0),
                                                 Eigen::RowVector3d(0, 1, 0),
                                                 Eigen::RowVector3d(0, 0, 1)};

  int axis = 0;

  double xx = 0;
  for (int i = 0; i < 3; ++i) {
    xx = abc_normal.dot(axis_normal[i]);
    if (fabs(abc_normal.dot(axis_normal[i])) > 0.5) {
      axis = i;
      break;
    }
  }

  int x = (axis + 1) % 3;
  int y = (axis + 2) % 3;

  Eigen::Vector2d p1, a1, b1, c1;
  p1(0) = p(x);
  p1(1) = p(y);
  a1(0) = a(x);
  a1(1) = a(y);
  b1(0) = b(x);
  b1(1) = b(y);
  c1(0) = c(x);
  c1(1) = c(y);

  double s1 = TiGER_GEOM_FUNC::orient2d(p1.data(), a1.data(), b1.data());
  double s2 = TiGER_GEOM_FUNC::orient2d(p1.data(), b1.data(), c1.data());
  double s3 = TiGER_GEOM_FUNC::orient2d(p1.data(), c1.data(), a1.data());
  if (s1 >= 0 && s2 >= 0 && s3 >= 0) return true;
  if (s1 <= 0 && s2 <= 0 && s3 <= 0) return true;
  // if (s1 * s2 *s3 == 0) return true; // intersected in a line // point in
  // edge of the extension
  // if (s1*s2>0&&s1*s3>0) {
  // return true;
  //}
  return false;
}
Eigen::RowVector3d closet_point_to_segment(const Eigen::RowVector3d& p,
                                           const Eigen::RowVector3d& s,
                                           const Eigen::RowVector3d& e) {
  /**    P(p)
  *     /|
  *    / |
  *   /  |
  *  /   |
  * /    |
  A(s)---C------------- B(e)
  */
  /* https://blog.csdn.net/qq_28087491/article/details/119239974 */
  Eigen::RowVector3d AB = e - s;
  Eigen::RowVector3d AP = p - s;
  double ratio = AP.dot(AB) / AB.squaredNorm();
  Eigen::RowVector3d AC = AB * ratio;
  if (ratio < 0) {
    return s;
  } else if (ratio > 1) {
    return e;
  } else {
    return AC + s;
  }

  return Eigen::RowVector3d();
}
// bool TriTriIntersection_tetgen(const Eigen::RowVector3d& a1,
//                                const Eigen::RowVector3d& b1,
//                                const Eigen::RowVector3d& c1,
//                                const Eigen::RowVector3d& a2,
//                                const Eigen::RowVector3d& b2,
//                                const Eigen::RowVector3d& c2) {
//   tetgenmesh::point pa1 = const_cast<double*>(a1.data());
//   tetgenmesh::point pb1 = const_cast<double*>(b1.data());
//   tetgenmesh::point pc1 = const_cast<double*>(c1.data());
//   tetgenmesh::point pa2 = const_cast<double*>(a2.data());
//   tetgenmesh::point pb2 = const_cast<double*>(b2.data());
//   tetgenmesh::point pc2 = const_cast<double*>(c2.data());
//   tetgenmesh t;
//   return t.tri_tri_inter(pa1, pb1, pc1, pa2, pb2, pc2);
// }
bool TriTriIntersection(const Eigen::RowVector3d& a1,
                        const Eigen::RowVector3d& b1,
                        const Eigen::RowVector3d& c1,
                        const Eigen::RowVector3d& a2,
                        const Eigen::RowVector3d& b2,
                        const Eigen::RowVector3d& c2) {
  REAL* real_a1 = const_cast<REAL*>(a1.data());
  REAL* real_b1 = const_cast<REAL*>(b1.data());
  REAL* real_c1 = const_cast<REAL*>(c1.data());
  REAL* real_a2 = const_cast<REAL*>(a2.data());
  REAL* real_b2 = const_cast<REAL*>(b2.data());
  REAL* real_c2 = const_cast<REAL*>(c2.data());

  double da1 = TiGER_GEOM_FUNC::orient3d(real_a1, real_a2, real_b2, real_c2);
  double db1 = TiGER_GEOM_FUNC::orient3d(real_b1, real_a2, real_b2, real_c2);
  double dc1 = TiGER_GEOM_FUNC::orient3d(real_c1, real_a2, real_b2, real_c2);

  double da2 = TiGER_GEOM_FUNC::orient3d(real_a2, real_a1, real_b1, real_c1);
  double db2 = TiGER_GEOM_FUNC::orient3d(real_b2, real_a1, real_b1, real_c1);
  double dc2 = TiGER_GEOM_FUNC::orient3d(real_c2, real_a1, real_b1, real_c1);

  if (da1 * db1 > 0 && da1 * dc1 > 0) return false;
  if (da2 * db2 > 0 && da2 * dc2 > 0) return false;

  if (lineTriIntersection(a1, b1, a2, b2, c2, da1, db1)) return true;
  if (lineTriIntersection(c1, b1, a2, b2, c2, dc1, db1)) return true;
  if (lineTriIntersection(a1, c1, a2, b2, c2, da1, dc1)) return true;
  if (lineTriIntersection(a2, b2, a1, b1, c1, da2, db2)) return true;
  if (lineTriIntersection(c2, b2, a1, b1, c1, dc2, db2)) return true;
  if (lineTriIntersection(a2, c2, a1, b1, c1, da2, dc2)) return true;
  return false;
}
double direction(const Eigen::Vector2d& pi, const Eigen::Vector2d& pj,
                 const Eigen::Vector2d& pk) {
  return (pi.x() - pk.x()) * (pi.y() - pj.y()) -
         (pi.x() - pj.x()) * (pi.y() - pk.y());
};

double angleDegree(const Eigen::RowVector3d& a, const Eigen::RowVector3d& b,
                   const Eigen::RowVector3d& c) {
  Eigen::RowVector3d ba = (a - b).normalized();
  Eigen::RowVector3d bc = (c - b).normalized();
  double cos_theta = ba.dot(bc);
  if (cos_theta > 1) cos_theta = 1;
  double ang = std::acos(cos_theta) * 180.0 / M_PI;
  return ang;
}

double angleSin(const Eigen::Vector2d& pi, const Eigen::Vector2d& pj,
                const Eigen::Vector2d& pk) {
  Eigen::Vector2d pkpi = (pk - pi).normalized();
  Eigen::Vector2d pjpi = (pj - pi).normalized();
  return (pkpi.x()) * (pjpi.y()) - (pjpi.x()) * (pkpi.y());
}

bool onSegment(const Eigen::Vector2d& pi, const Eigen::Vector2d& pj,
               const Eigen::Vector2d& pk) {
  if ((pk - pi).norm() < 1e-10 || (pk - pj).norm() < 1e-10) return false;
  return ((std::min(pi.x(), pj.x()) <= pk.x() &&
           pk.x() <= std::max(pi.x(), pj.x()) &&
           std::min(pi.y(), pj.y()) <= pk.y() &&
           pk.y() <= std::max(pi.y(), pj.y())));
};

bool checkIntersection(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                       const Eigen::Vector2d& c, const Eigen::Vector2d& d) {
  double eps = 1e-13;

  REAL* real_a = const_cast<REAL*>(a.data());
  REAL* real_b = const_cast<REAL*>(b.data());
  REAL* real_c = const_cast<REAL*>(c.data());
  REAL* real_d = const_cast<REAL*>(d.data());

  double d1 = TiGER_GEOM_FUNC::orient2d(real_c, real_d, real_a);
  double d2 = TiGER_GEOM_FUNC::orient2d(real_c, real_d, real_b);
  double d3 = TiGER_GEOM_FUNC::orient2d(real_a, real_b, real_c);
  double d4 = TiGER_GEOM_FUNC::orient2d(real_a, real_b, real_d);

  if (fabs(d1) < eps) d1 = 0;
  if (fabs(d2) < eps) d2 = 0;
  if (fabs(d3) < eps) d3 = 0;
  if (fabs(d4) < eps) d4 = 0;

  if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
      ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
    return true;
  }
  if (fabs(d1) == 0 && onSegment(c, d, a)) {
    return true;
  }
  if (fabs(d2) == 0 && onSegment(c, d, b)) {
    return true;
  }
  if (fabs(d3) == 0 && onSegment(a, b, c)) {
    return true;
  }
  if (fabs(d4) == 0 && onSegment(a, b, d)) {
    return true;
  }
  return false;
};
}  // namespace common
}  // namespace TiGER