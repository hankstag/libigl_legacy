// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "map_vertices_to_circle.h"
template <typename DerivedV, typename DerivedF>
  IGL_INLINE void igl::map_vertices_to_circle(
  	const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedF>& bnd,
  	Eigen::PlainObjectBase<DerivedV>& UV)
{
  // Get sorted list of boundary vertices
  typedef typename DerivedF::Scalar ScalarF;
  typedef typename DerivedV::Scalar ScalarV;
  
  std::vector<ScalarF> interior,map_ij;
  map_ij.resize(V.rows());
  std::vector<bool> isOnBnd(V.rows(),false);
  
  for (int i = 0; i < bnd.size(); i++)
  {
    isOnBnd[bnd[i]] = true;
    map_ij[bnd[i]] = static_cast<ScalarF>(i);
  }

  for (int i = 0; i < (int)isOnBnd.size(); i++)
  {
    if (!isOnBnd[i])
    {
      map_ij[i] = static_cast<ScalarF>(interior.size());
      interior.push_back(static_cast<ScalarF>(i));
    }
  }
  
  // Map boundary to unit circle
  std::vector<ScalarV> len(bnd.size());
  len[0] = static_cast<ScalarV>(0.);

  for (int i = 1; i < bnd.size(); i++)
  {
    len[i] = len[i-1] + (V.row(bnd[i-1]) - V.row(bnd[i])).norm();
  }

  ScalarV total_len = len[len.size()-1] + (V.row(bnd[0]) - V.row(bnd[bnd.size()-1])).norm();

  UV.resize(bnd.rows(),2);
  for (int i = 0; i < bnd.size(); i++)
  {
    ScalarV frac = len[i] * 2. * static_cast<ScalarV>(M_PI) / total_len;
    UV.row(map_ij[bnd[i]]) << cos(frac), sin(frac);
  }

}

#ifdef IGL_STATIC_LIBRARY
template void igl::map_vertices_to_circle<Eigen::Matrix<mpfr::mpreal, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<mpfr::mpreal, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<mpfr::mpreal, -1, -1, 0, -1, -1> >&);
template void igl::map_vertices_to_circle<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
//template void igl::map_vertices_to_circle<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
//template void igl::map_vertices_to_circle<Eigen::Matrix<mpfr::mpreal, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<mpfr::mpreal, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<mpfr::mpreal, -1, -1, 0, -1, -1> >&);
#endif