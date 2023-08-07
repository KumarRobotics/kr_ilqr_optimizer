//
// Created by Zixin Zhang on 03/17/23.
// Copyright (c) 2022 Robotic Exploration Lab. All rights reserved.
//

#include "quaternion_utils.hpp"
#include <iostream>

namespace altro {
Eigen::Matrix3d skew(Eigen::Vector3d vec) {
  Eigen::Matrix3d skew;
  skew << 0, -vec[2], vec[1], vec[2], 0, -vec[0], -vec[1], vec[0], 0;
  return skew;
}

Eigen::Vector4d cayley_map(Eigen::Vector3d phi) {
  Eigen::Vector4d phi_quat;
  phi_quat << 1, phi[0], phi[1], phi[2];
  return 1 / (sqrt(1 + phi.norm() * phi.norm())) * phi_quat;
}

Eigen::Vector3d inv_cayley_map(Eigen::Vector4d q) { return q.tail(3) / q[0]; }

Eigen::Vector4d quat_conj(Eigen::Vector4d q) {
  Eigen::Matrix4d T;
  T << 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1;
  return T * q;
}

Eigen::Matrix4d T(){
  Eigen::Matrix4d T = -Eigen::Matrix4d::Identity();
  T(0, 0) = 1;
  return T;
}

Eigen::MatrixXd H(){
  Eigen::Matrix<double, 4,3> H_con;
  H_con << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;
  return H_con;
}

Eigen::MatrixXd dQdq(Eigen::Vector4d q){
  double qw = q[0];
  double qx = q[1];
  double qy = q[2];
  double qz = q[3];
  Eigen::Matrix<double, 9, 4> dR_dq;
//row major
  dR_dq(0, 0) = 0;
  dR_dq(0, 1) = 0;
  dR_dq(0, 2) = -4 * qy;
  dR_dq(0, 3) = -4 * qz;

  dR_dq(1, 0) = -2 * qz;
  dR_dq(1, 1) = 2 * qy;
  dR_dq(1, 2) = 2 * qx;
  dR_dq(1, 3) = -2 * qw;

  dR_dq(2, 0) = 2 * qy;
  dR_dq(2, 1) = 2 * qz;
  dR_dq(2, 2) = 2 * qw;
  dR_dq(2, 3) = 2 * qx;

  dR_dq(3, 0) = 2 * qz;
  dR_dq(3, 1) = 2 * qy;
  dR_dq(3, 2) = 2 * qx;
  dR_dq(3, 3) = 2 * qw;

  dR_dq(4, 0) = 0;
  dR_dq(4, 1) = -4 * qx;
  dR_dq(4, 2) = 0;
  dR_dq(4, 3) = -4 * qz;

  dR_dq(5, 0) = -2 * qx;
  dR_dq(5, 1) = -2 * qw;
  dR_dq(5, 2) = 2 * qz;
  dR_dq(5, 3) = 2 * qy;

  dR_dq(6, 0) = -2 * qy;
  dR_dq(6, 1) = 2 * qz;
  dR_dq(6, 2) = -2 * qw;
  dR_dq(6, 3) = 2 * qx;

  dR_dq(7, 0) = 2 * qx;
  dR_dq(7, 1) = 2 * qw;
  dR_dq(7, 2) = 2 * qz;
  dR_dq(7, 3) = 2 * qy;

  dR_dq(8, 0) = 0;
  dR_dq(8, 1) = -4 * qx;
  dR_dq(8, 2) = -4 * qy;
  dR_dq(8, 3) = 0;
  return dR_dq;
}

Eigen::Matrix4d L(Eigen::Vector4d q) {
  Eigen::Matrix4d L;
  L(0, 0) = q[0];
  L.block<1, 3>(0, 1) = -q.tail(3).transpose();
  L.block<3, 1>(1, 0) = q.tail(3);
  L.block<3, 3>(1, 1) = q[0] * Eigen::Matrix3d::Identity() + skew(q.tail(3));
  // std::cout << "L: success" << std::endl;
  return L;
}

Eigen::Matrix4d R(Eigen::Vector4d q) {
  Eigen::Matrix4d R;
  R(0, 0) = q[0];
  R.block<1, 3>(0, 1) = -q.tail(3).transpose();
  R.block<3, 1>(1, 0) = q.tail(3);
  R.block<3, 3>(1, 1) = q[0] * Eigen::Matrix3d::Identity() - skew(q.tail(3));
  return R;
}

// T = Diagonal([1; -ones(3)])
// H = [zeros(1,3); I]
// function qtoQ(q)
//     return H'*T*L(q)*T*L(q)*H
// end
Eigen::Matrix3d Q(Eigen::Vector4d q) {
  return H().transpose()* T() * L(q) *T() * L(q) * H();
}

Eigen::MatrixXd G(Eigen::Vector4d q) {
  return L(q) * H();
}

Eigen::Vector4d quat_mult(Eigen::Vector4d q1, Eigen::Vector4d q2) { return L(q1) * q2; }

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove) {
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();

  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
        matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

  matrix.conservativeResize(numRows, numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove) {
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols() - 1;

  if (colToRemove < numCols)
    matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
        matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);

  matrix.conservativeResize(numRows, numCols);
}

Eigen::Matrix<double, -1, 1> removeElement(const Eigen::Matrix<double, -1, 1>& vec,
                                           unsigned int idxToRemove) {
  unsigned int length = vec.rows() - 1;
  Eigen::Matrix<double, -1, 1> vec_;
  vec_ = vec;

  if (idxToRemove < length) {
    vec_.block(idxToRemove, 0, length - idxToRemove, 1) =
        vec_.block(idxToRemove + 1, 0, length - idxToRemove, 1);
  }
  vec_.conservativeResize(length, 1);

  return vec_;
}
}  // namespace altro
