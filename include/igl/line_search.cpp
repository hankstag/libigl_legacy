// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Michael Rabinovich
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "line_search.h"
#include <iostream>

IGL_INLINE double igl::line_search(
        Eigen::MatrixXd& x,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXd& d,
        double step_size,
        std::function<double(Eigen::MatrixXd&)> energy,
        double cur_energy) {
  double old_energy;
//  if (cur_energy > 0)
//  {
//    old_energy = cur_energy;
//  }
//  else
//  {
  old_energy = energy(x); // no energy was given -> need to compute the current energy
//  }
  assert(!std::isinf(old_energy) && !std::isnan(old_energy));
  double new_energy = old_energy;
  int cur_iter = 0;
  int MAX_STEP_SIZE_ITER = 80;

  while (new_energy >= old_energy) {
    if (cur_iter > MAX_STEP_SIZE_ITER) {
            std::cout << "line_search.cpp runs out of iterations! " << std::endl;
      break;
    }

    Eigen::MatrixXd new_x = x + step_size * d;

    // explicitly check orientation
    bool flipped = false;
    // for 2D case
    if(x.cols()==2) {
      Eigen::Vector3d a, b;
      assert(F.rows() > 0);
      a << (new_x.row(F(0, 0)) - new_x.row(F(0, 1))).transpose(), 0;
      b << (new_x.row(F(0, 0)) - new_x.row(F(0, 2))).transpose(), 0;
      double k = a.cross(b)(2);
      for (int i = 0; i < F.rows(); i++) {
        a << (new_x.row(F(i, 0)) - new_x.row(F(i, 1))).transpose(), 0;
        b << (new_x.row(F(i, 0)) - new_x.row(F(i, 2))).transpose(), 0;
        double o = a.cross(b)(2);
        if ((o > 0 && k < 0) || (o < 0 && k > 0)) {
          flipped = true;
          break;
        }
      }
    }

    double cur_e = energy(new_x);
//        std::cerr<<"it "<<cur_iter<<": "<<cur_e<<" "<<old_energy<<std::endl;
    if (std::isnan(cur_e) || std::isinf(cur_e) || cur_e > old_energy || flipped) {
      step_size /= 2;
    } else {
      x = new_x;
      new_energy = cur_e;
      break;
    }
    cur_iter++;
  }

  return new_energy;
}


#ifdef IGL_STATIC_LIBRARY
#endif