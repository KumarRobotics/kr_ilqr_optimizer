//
// Created by Brian Jackson on 9/27/22.
// Copyright (c) 2022 Robotic Exploration Lab. All rights reserved.
//

#include "kr_ilqr_optimizer/test_utils.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "kr_ilqr_optimizer/finitediff.hpp"

namespace fs = std::filesystem;
// using json = nlohmann::json;

void discrete_double_integrator_dynamics(
    double* xnext, const double* x, const double* u, float h, int dim) {
  double b = h * h / 2;
  for (int i = 0; i < dim; ++i) {
    xnext[i] = x[i] + x[i + dim] * h + u[i] * b;
    xnext[i + dim] = x[i + dim] + u[i] * h;
  }
}

void discrete_double_integrator_jacobian(
    double* jac, const double* x, const double* u, float h, int dim) {
  (void)x;
  (void)u;
  Eigen::Map<Eigen::MatrixXd> J(jac, 2 * dim, 3 * dim);
  J.setZero();
  double b = h * h / 2;
  for (int i = 0; i < dim; ++i) {
    J(i, i) = 1.0;
    J(i + dim, i + dim) = 1.0;
    J(i, i + dim) = h;
    J(i, 2 * dim + i) = b;
    J(i + dim, 2 * dim + i) = h;
  }
}

const double kPendulumMass = 1.0;
const double kPendulumLength = 0.5;
const double kPendulumFrictionCoeff = 0.1;
// const double kPendulumInertia = 0.25;
const double kPendulumGravity = 9.81;

void pendulum_dynamics(double* xnext, const double* x, const double* u) {
  double l = kPendulumLength;
  double g = kPendulumGravity;
  double b = kPendulumFrictionCoeff;
  double m = kPendulumMass * l * l;

  double theta = x[0];
  double omega = x[1];

  double omega_dot = u[0] / m - g * std::sin(theta) / l - b * omega / m;
  xnext[0] = omega;
  xnext[1] = omega_dot;
}

void pendulum_jacobian(double* jac, const double* x, const double* u) {
  (void)u;
  double l = kPendulumLength;
  double g = kPendulumGravity;
  double b = kPendulumFrictionCoeff;
  double m = kPendulumMass * l * l;

  double domega_dtheta = 0.0;
  double domega_domega = 1.0;
  double domega_du = 0.0;
  double dalpha_dtheta = -g * std::cos(x[0]) / l;
  double dalpha_domega = -b / m;
  double dalpha_du = 1 / m;
  jac[0] = domega_dtheta;
  jac[1] = dalpha_dtheta;
  jac[2] = domega_domega;
  jac[3] = dalpha_domega;
  jac[4] = domega_du;
  jac[5] = dalpha_du;
}

altro::ExplicitDynamicsFunction MidpointDynamics(int n,
                                                 int m,
                                                 ContinuousDynamicsFunction f) {
  (void)m;
  auto fd = [n, f](double* xn, const double* x, const double* u, float h) {
    static Eigen::VectorXd xm(n);
    Eigen::Map<Eigen::VectorXd> xn_vec(xn, n);
    Eigen::Map<const Eigen::VectorXd> x_vec(x, n);
    Eigen::Map<const Eigen::VectorXd> u_vec(u, n);
    f(xm.data(), x, u);
    xm *= h / 2;
    xm.noalias() += x_vec;
    f(xn, xm.data(), u);
    xn_vec = x_vec + h * xn_vec;
  };
  return fd;
}

altro::ExplicitDynamicsJacobian MidpointJacobian(
    int n, int m, ContinuousDynamicsFunction f, ContinuousDynamicsJacobian df) {
  auto fd = [n, m, f, df](
                double* jac, const double* x, const double* u, float h) {
    static Eigen::MatrixXd A(n, n);
    static Eigen::MatrixXd B(n, m);
    static Eigen::MatrixXd Am(n, n);
    static Eigen::MatrixXd Bm(n, m);
    static Eigen::VectorXd xm(n);
    static Eigen::MatrixXd In = Eigen::MatrixXd::Identity(n, n);

    Eigen::Map<Eigen::MatrixXd> J(jac, n, n + m);
    Eigen::Map<const Eigen::VectorXd> x_vec(x, n);
    Eigen::Map<const Eigen::VectorXd> u_vec(u, n);

    // Evaluate the midpoint
    f(xm.data(), x, u);
    xm = x_vec + h / 2 * xm;

    // Evaluate the Jacobian
    df(J.data(), x, u);
    A = J.leftCols(n);
    B = J.rightCols(m);

    // Evaluate the Jacobian at the midpoint
    df(J.data(), xm.data(), u);
    Am = J.leftCols(n);
    Bm = J.rightCols(m);

    // Apply the chain rule
    J.leftCols(n) = In + h * Am * (In + h / 2 * A);
    J.rightCols(m) = h * (Am * h / 2 * B + Bm);
  };
  return fd;
}

altro::ExplicitDynamicsFunction ForwardEulerDynamics(
    int n, int m, const ContinuousDynamicsFunction f) {
  (void)m;
  auto fd = [n, f](double* xn, const double* x, const double* u, float h) {
    Eigen::Map<Eigen::VectorXd> xn_vec(xn, n);
    Eigen::Map<const Eigen::VectorXd> x_vec(x, n);
    Eigen::Map<const Eigen::VectorXd> u_vec(u, n);
    f(xn, x, u);  // xn is actually x_dot here
    xn_vec = x_vec + h * xn_vec;
  };
  return fd;
}

altro::ExplicitDynamicsJacobian ForwardEulerJacobian(
    int n,
    int m,
    const ContinuousDynamicsFunction f,
    const ContinuousDynamicsJacobian df) {
  auto fd = [n, m, f, df](
                double* jac, const double* x, const double* u, float h) {
    Eigen::Map<Eigen::MatrixXd> J(jac, n, n + m);
    Eigen::Map<const Eigen::VectorXd> x_vec(x, n);
    Eigen::Map<const Eigen::VectorXd> u_vec(u, n);

    static Eigen::MatrixXd In = Eigen::MatrixXd::Identity(n, n);

    df(J.data(), x, u);
    J.leftCols(n) = In + h * J.leftCols(n);
    J.rightCols(m) = h * J.rightCols(m);
  };
  return fd;
}

void BicycleModel::Dynamics(double* x_dot,
                            const double* x,
                            const double* u) const {
  double v = u[0];          // longitudinal velocity (m/s)
  double delta_dot = u[1];  // steering angle rage (rad/s)
  double theta = x[2];      // heading angle (rad) relative to x-axis
  double delta = x[3];      // steering angle (rad)

  double beta = 0;
  double omega = 0;
  double stheta = 0;
  double ctheta = 0;
  switch (reference_frame_) {
    case ReferenceFrame::CenterOfGravity:
      beta = std::atan2(distance_to_rear_wheels_ * delta, length_);
      omega = v * std::cos(beta) * std::tan(delta) / length_;
      stheta = std::sin(theta + beta);
      ctheta = std::cos(theta + beta);
      break;
    case ReferenceFrame::Rear:
      omega = v * tan(delta) / length_;
      stheta = std::sin(theta);
      ctheta = std::cos(theta);
      break;
    case ReferenceFrame::Front:
      omega = v * std::sin(delta) / length_;
      stheta = std::sin(theta + delta);
      ctheta = std::cos(theta + delta);
      break;
  };
  double px_dot = v * ctheta;
  double py_dot = v * stheta;
  x_dot[0] = px_dot;
  x_dot[1] = py_dot;
  x_dot[2] = omega;
  x_dot[3] = delta_dot;
}

void BicycleModel::Jacobian(double* jac,
                            const double* x,
                            const double* u) const {
  double v = u[0];      // longitudinal velocity (m/s)
  double theta = x[2];  // heading angle (rad) relative to x-axis
  double delta = x[3];  // steering angle (rad)

  Eigen::Map<Eigen::Matrix<double, 4, 6>> J(jac);
  double beta = 0;
  double dbeta_ddelta = 0;
  double by = 0;
  double bx = 0;
  double domega_ddelta = 0;
  double domega_dv = 0;

  double stheta = 0;
  double ctheta = 0;
  double ds_dtheta = 0;
  double dc_dtheta = 0;
  double ds_ddelta = 0;
  double dc_ddelta = 0;
  switch (reference_frame_) {
    case ReferenceFrame::CenterOfGravity:
      by = distance_to_rear_wheels_ * delta;
      bx = length_;
      beta = std::atan2(by, bx);
      dbeta_ddelta = bx / (bx * bx + by * by) * distance_to_rear_wheels_;
      domega_ddelta = v / length_ *
                      (-std::sin(beta) * std::tan(delta) * dbeta_ddelta +
                       std::cos(beta) / (std::cos(delta) * std::cos(delta)));
      domega_dv = std::cos(beta) * std::tan(delta) / length_;

      stheta = std::sin(theta + beta);
      ctheta = std::cos(theta + beta);
      ds_dtheta = +std::cos(theta + beta);
      dc_dtheta = -std::sin(theta + beta);
      ds_ddelta = +std::cos(theta + beta) * dbeta_ddelta;
      dc_ddelta = -std::sin(theta + beta) * dbeta_ddelta;
      break;
    case ReferenceFrame::Rear:
      domega_ddelta = v / length_ / (std::cos(delta) * std::cos(delta));
      domega_dv = std::tan(delta) / length_;

      stheta = std::sin(theta);
      ctheta = std::cos(theta);
      ds_dtheta = +std::cos(theta);
      dc_dtheta = -std::sin(theta);
      break;
    case ReferenceFrame::Front:
      domega_ddelta = v / length_ * std::cos(delta);
      domega_dv = std::sin(delta) / length_;

      stheta = std::sin(theta + delta);
      ctheta = std::cos(theta + delta);
      ds_dtheta = +std::cos(theta + delta);
      dc_dtheta = -std::sin(theta + delta);
      ds_ddelta = ds_dtheta;
      dc_ddelta = dc_dtheta;
      break;
  };
  J.setZero();
  J(0, 2) = v * dc_dtheta;  // dxdot_dtheta
  J(0, 3) = v * dc_ddelta;  // dxdot_ddelta
  J(0, 4) = ctheta;         // dxdot_dv
  J(1, 2) = v * ds_dtheta;  // dydot_dtheta
  J(1, 3) = v * ds_ddelta;  // dydot_ddelta
  J(1, 4) = stheta;         // dydot_dv
  J(2, 3) = domega_ddelta;
  J(2, 4) = domega_dv;
  J(3, 5) = 1.0;
}

void SimpleQuaternionModel::Dynamics(double* x_dot,
                                     const double* x,
                                     const double* u) const {
  Eigen::Map<Eigen::VectorXd> x_dot_vec(x_dot, 4);
  Eigen::Map<const Eigen::VectorXd> x_vec(x, 4);
  Eigen::Map<const Eigen::VectorXd> u_vec(u, 3);

  x_dot_vec = 0.5 * altro::G(x_vec) * u_vec;
}

void SimpleQuaternionModel::Jacobian(double* jac,
                                     const double* x,
                                     const double* u) const {
  Eigen::Map<const Eigen::VectorXd> x_vec(x, 4);
  Eigen::Map<const Eigen::VectorXd> u_vec(u, 3);

  double qs = x_vec[0];
  double qa = x_vec[1];
  double qb = x_vec[2];
  double qc = x_vec[3];

  double wx = u_vec[0];
  double wy = u_vec[1];
  double wz = u_vec[2];

  Eigen::Map<Eigen::Matrix<double, 4, 7>> J(jac);
  J.setZero();
  J << 0, -wx / 2, -wy / 2, -wz / 2, -qa / 2, -qb / 2, -qc / 2, wx / 2, 0,
      wz / 2, -wy / 2, qs / 2, -qc / 2, qb / 2, wy / 2, -wz / 2, 0, wx / 2,
      qc / 2, qs / 2, -qa / 2, wz / 2, wy / 2, -wx / 2, 0, -qb / 2, qa / 2,
      qs / 2;
}

Eigen::Vector3d quadModel::moments(const Eigen::VectorXd& u) const {
  double L = motor_dist_;

  double w1 = u[0];
  double w2 = u[1];
  double w3 = u[2];
  double w4 = u[3];

  double F1 = std::max(0.0, kf_ * w1);
  double F2 = std::max(0.0, kf_ * w2);
  double F3 = std::max(0.0, kf_ * w3);
  double F4 = std::max(0.0, kf_ * w4);

  double M1 = km_ * w1;
  double M2 = km_ * w2;
  double M3 = km_ * w3;
  double M4 = km_ * w4;

  Eigen::Vector3d tau;
  tau << L * (F2 - F4), L * (F3 - F1), (M1 - M2 + M3 - M4);

  return tau;
}

Eigen::Vector3d quadModel::forces(const Eigen::VectorXd& u) const {
  double w1 = u[0];
  double w2 = u[1];
  double w3 = u[2];
  double w4 = u[3];

  double F1 = std::max(0.0, kf_ * w1);
  double F2 = std::max(0.0, kf_ * w2);
  double F3 = std::max(0.0, kf_ * w3);
  double F4 = std::max(0.0, kf_ * w4);

  Eigen::Vector3d F;
  F << 0.0, 0.0, F1 + F2 + F3 + F4;  // total rotor force in body frame
  return F;
}

void quadModel::Dynamics(double* x_dot,
                         const double* x,
                         const double* u) const {
  // std::cout << "section 0 success "  << std::endl;

  Eigen::Map<Eigen::VectorXd> x_dot_vec(x_dot, 13);
  x_dot_vec.segment<13>(0) = Eigen::VectorXd::Zero(13);
  Eigen::Map<const Eigen::VectorXd> x_vec(x, 13);
  Eigen::Map<const Eigen::VectorXd> u_vec(u, 4);

  const double robot_mass = mass_;
  Eigen::Vector3d g_vec;
  g_vec << 0, 0, -grav_;

  // std::cout << "section 1 success "  << std::endl;

  Eigen::Vector3d moment_body = moments(u_vec);
  Eigen::Vector3d force_body = forces(u_vec);

  // std::cout << "section 2 success "  << std::endl;

  // state: pos(world) quat(world->body) vel(body) omega(body)
  //        0 1 2      3  4  5  6        7 8 9     10 11 12
  Eigen::Vector4d q = x_vec.segment<4>(3).normalized();
  Eigen::Vector3d vb = x_vec.segment<3>(7);
  Eigen::Vector3d w = x_vec.segment<3>(10);
  // std::cout << "section 3 success "  << std::endl;

  Eigen::Matrix3d Q = altro::Q(q);
  // change rate of position
  x_dot_vec.segment<3>(0) = Q * vb;

  // change rate of quaternion
  // std::cout << "section 4 success "  << std::endl;

  x_dot_vec.segment<4>(3) = 0.5 * altro::G(q) * w;
  // std::cout << "section 5 success "  << std::endl;

  // change rate of linear velocity

  x_dot_vec.segment<3>(7) =
      Q.transpose() * g_vec + force_body / robot_mass - w.cross(vb);
  // std::cout << "section 6 success "  << std::endl;

  // change rate of angular velocity
  x_dot_vec.segment<3>(10) = moment_of_inertia_.inverse() *
                             (moment_body - w.cross(moment_of_inertia_ * w));
  // std::cout << "section 7 success " << std::endl;
}

void quadModel::finite_jacobian_quad_xu(double* jac,
                                        const double* x,
                                        const double* u) const {
  // Eigen::Map<Eigen::VectorXd> x_dot_vec(x_dot, 13);
  Eigen::Map<const Eigen::VectorXd> x_vec(x, 13);
  Eigen::Map<const Eigen::VectorXd> u_vec(u, 4);
  Eigen::Map<Eigen::Matrix<double, 13, 17>> J(jac);  // jac = [dfc_dx, dfc_du]

  Eigen::MatrixXd dfc_dx(13, 13);
  Eigen::MatrixXd dfc_du(13, 4);

  const double eps = 1e-10;
  fd::AccuracyOrder accuracy = fd::EIGHTH;
  const std::vector<double> external_coeffs = fd::get_external_coeffs(accuracy);
  const std::vector<double> internal_coeffs = fd::get_interior_coeffs(accuracy);
  assert(external_coeffs.size() == internal_coeffs.size());
  const size_t inner_steps = internal_coeffs.size();
  const double denom = get_denominator(accuracy) * eps;
  J.setZero();
  dfc_dx.setZero();
  dfc_du.setZero();
  // std::cout<< "number of rows" << x_vec.rows()<<std::endl;
  // std::cout<< "Inner steps" << inner_steps<<std::endl;
  // std::cout<< "dfc_dx.col(7)"<<dfc_dx.col(7)<<std::endl;

  Eigen::VectorXd x_mutable = x_vec;
  for (size_t i = 0; i < x_vec.rows(); i++) {
    for (size_t ci = 0; ci < inner_steps; ci++) {
      x_mutable[i] += internal_coeffs[ci] * eps;
      // normalize quaternion
      // x_mutable.segment<4>(3).normalize();
      Eigen::VectorXd fx(13);

      Dynamics(fx.data(), x_mutable.data(), u);
      // fx = Eigen::VectorXd::Zero(13);
      // std::cout << "fx" << fx << std::endl;

      dfc_dx.col(i) += external_coeffs[ci] * fx;
      for (size_t j = 0; j < x_vec.rows(); j++) {
        x_mutable[j] = x_vec[j];
        // reset x_mutable, it is unclear if quat will make things weird, so
        // just reset everything
      }
    }
    dfc_dx.col(i) /= denom;
  }

  Eigen::VectorXd u_mutable = u_vec;
  for (size_t i = 0; i < u_vec.rows(); i++) {
    for (size_t ci = 0; ci < inner_steps; ci++) {
      u_mutable[i] += internal_coeffs[ci] * eps;
      Eigen::VectorXd fx(13);

      Dynamics(fx.data(), x, u_mutable.data());
      dfc_du.col(i) += external_coeffs[ci] * fx;
      u_mutable[i] = u_vec[i];
    }
    dfc_du.col(i) /= denom;
  }

  // make all values < 1e-10 equal to zero
  for (size_t i = 0; i < dfc_dx.rows(); i++) {
    for (size_t j = 0; j < dfc_dx.cols(); j++) {
      if (std::abs(dfc_dx(i, j)) < 1e-8) {
        dfc_dx(i, j) = 0;
      }
    }
  }
  J.block<13, 13>(0, 0) = dfc_dx;
  J.block<13, 4>(0, 13) = dfc_du;
}

void quadModel::Jacobian_fd(double* jac,
                            const double* x,
                            const double* u) const {
  finite_jacobian_quad_xu(jac, x, u);
}
Eigen::Matrix<double, 3, 4> quat_a_jacobian(const Eigen::Vector4d& q,
                                            const Eigen::Vector3d& g) {
  Eigen::Matrix<double, 3, 4> mat;

  mat(0, 0) = -2 * g(2) * q(2);
  mat(0, 1) = 2 * g(2) * q(3);
  mat(0, 2) = -2 * g(2) * q(0);
  mat(0, 3) = 2 * g(2) * q(1);

  mat(1, 0) = 2 * g(2) * q(1);
  mat(1, 1) = 2 * g(2) * q(0);
  mat(1, 2) = 2 * g(2) * q(3);
  mat(1, 3) = 2 * g(2) * q(2);

  mat(2, 0) = 2 * g(2) * q(0);
  mat(2, 1) = -2 * g(2) * q(1);
  mat(2, 2) = -2 * g(2) * q(2);
  mat(2, 3) = 2 * g(2) * q(3);
  return mat;
}

Eigen::Matrix<double, 3, 4> quat_v_jacobian(const Eigen::Vector4d& q,
                                            const Eigen::Vector3d& v) {
  Eigen::Matrix<double, 3, 4> mat;  // MATLAB symbolic toolbox

  mat(0, 0) = 2 * q(0) * v(0) + 2 * q(2) * v(2) - 2 * q(3) * v(1);
  mat(0, 1) = 2 * q(1) * v(0) + 2 * q(2) * v(1) + 2 * q(3) * v(2);
  mat(0, 2) = 2 * q(0) * v(2) + 2 * q(1) * v(1) - 2 * q(2) * v(0);
  mat(0, 3) = 2 * q(1) * v(2) - 2 * q(0) * v(1) - 2 * q(3) * v(0);

  mat(1, 0) = 2 * q(0) * v(1) - 2 * q(1) * v(2) + 2 * q(3) * v(0);
  mat(1, 1) = 2 * q(2) * v(0) - 2 * q(1) * v(1) - 2 * q(0) * v(2);
  mat(1, 2) = 2 * q(1) * v(0) + 2 * q(2) * v(1) + 2 * q(3) * v(2);
  mat(1, 3) = 2 * q(0) * v(0) + 2 * q(2) * v(2) - 2 * q(3) * v(1);

  mat(2, 0) = 2 * q(0) * v(2) + 2 * q(1) * v(1) - 2 * q(2) * v(0);
  mat(2, 1) = 2 * q(0) * v(1) - 2 * q(1) * v(2) + 2 * q(3) * v(0);
  mat(2, 2) = 2 * q(3) * v(1) - 2 * q(2) * v(2) - 2 * q(0) * v(0);
  mat(2, 3) = 2 * q(1) * v(0) + 2 * q(2) * v(1) + 2 * q(3) * v(2);
  return mat;
}

void quadModel::Jacobian(double* jac, const double* x, const double* u) const {
  Eigen::Map<Eigen::Matrix<double, 13, 17>> J(jac);  // jac = [dfc_dx, dfc_du]
  J.setZero();

  Eigen::Map<const Eigen::VectorXd> x_vec(x, 13);
  Eigen::Map<const Eigen::VectorXd> u_vec(u, 4);

  Eigen::Vector4d q = x_vec.segment<4>(3);
  Eigen::Vector3d vb = x_vec.segment<3>(7);
  Eigen::Vector3d w = x_vec.segment<3>(10);
  Eigen::Matrix3d Q = altro::Q(q);
  Eigen::Vector4d inv_q = altro::quat_conj(q);

  const double robot_mass = mass_;
  Eigen::Vector3d g_vec;
  g_vec << 0, 0, -grav_;

  // Calculate dfc_dx
  Eigen::MatrixXd dfc_dx(13, 13);

  dfc_dx.setZero();
  // dvw/dq ????
  // Eigen::MatrixXd dQdq = altro::dQdq(q);
  //                           //1 x 3             3x4
  // dfc_dx.block<1, 4>(0, 3) = (vb.transpose() * dQdq.block<3, 4>(0, 0));
  // dfc_dx.block<1, 4>(1, 3) = (vb.transpose() * dQdq.block<3, 4>(3, 0));
  // dfc_dx.block<1, 4>(2, 3) = (vb.transpose() * dQdq.block<3, 4>(6, 0));
  double q_w = q[0];
  Eigen::Vector3d q_vec = q.segment<3>(1);
  // Eigen::Matrix<double, 3, 4> JJ;
  // JJ.setZero();

  // JJ.block<3,1>(0,3) = 2.0*(q_w*vb+q_vec.cross(vb));
  // JJ.block<3,3>(0,4) = 2.0*(q_vec.transpose()*vb*Eigen::MatrixXd::Identity(3,
  // 3)+q_vec*vb.transpose()- vb*q_vec.transpose()-q_w*altro::skew(vb));
  dfc_dx.block<3, 4>(0, 3) = quat_v_jacobian(x_vec.segment<4>(3), vb);

  // dvw/dvb
  dfc_dx.block<3, 3>(0, 7) = Q;
  // dqdot/dq
  dfc_dx.block<1, 3>(3, 4) = -0.5 * w.transpose();
  dfc_dx.block<3, 1>(4, 3) = 0.5 * w;
  dfc_dx.block<3, 3>(4, 4) = -0.5 * altro::skew(w);
  // dqdot/domega // is this just 0.5*G(q)?
  dfc_dx(3, 10) = -0.5 * x_vec(4);  // -0.5qa
  dfc_dx(3, 11) = -0.5 * x_vec(5);  // -0.5qb
  dfc_dx(3, 12) = -0.5 * x_vec(6);  // -0.5qc
  dfc_dx(4, 10) = 0.5 * x_vec(3);   // 0.5qs
  dfc_dx(4, 11) = -0.5 * x_vec(6);  // -0.5qc
  dfc_dx(4, 12) = 0.5 * x_vec(5);   // 0.5qb
  dfc_dx(5, 10) = 0.5 * x_vec(6);   // 0.5qc
  dfc_dx(5, 11) = 0.5 * x_vec(3);   // 0.5qs
  dfc_dx(5, 12) = -0.5 * x_vec(4);  // -0.5qa
  dfc_dx(6, 10) = -0.5 * x_vec(5);  // -0.5qb
  dfc_dx(6, 11) = 0.5 * x_vec(4);   // 0.5qa
  dfc_dx(6, 12) = 0.5 * x_vec(3);   // 0.5qs

  // dvbdot/dq // ADDED - sign here, not sure if correct
  //  Eigen::MatrixXd dQinv_dq = altro::dQdq(inv_q);
  //  dfc_dx.block<1, 4>(7, 3) = -(g_vec.transpose() * dQinv_dq.block<3, 4>(0,
  //  0)); dfc_dx.block<1, 4>(8, 3) = -(g_vec.transpose() * dQinv_dq.block<3,
  //  4>(3, 0)); dfc_dx.block<1, 4>(9, 3) = -(g_vec.transpose() *
  //  dQinv_dq.block<3, 4>(6, 0));

  // Eigen::Matrix<double, 3, 4> JJ2;
  // JJ2.setZero();
  // JJ2.block<3, 1>(0, 3) = 2.0 * (q_w * g_vec + q_vec.cross(g_vec));
  // JJ2.block<3, 3>(0, 4) =
  //     2.0 * (-q_vec.transpose() * g_vec * Eigen::MatrixXd::Identity(3, 3) +
  //            -q_vec * g_vec.transpose() - g_vec * -q_vec.transpose() - q_w *
  //            altro::skew(g_vec));
  dfc_dx.block<3, 4>(7, 3) = quat_a_jacobian(x_vec.segment<4>(3), g_vec);

  // dvbdot/dvb
  dfc_dx.block<3, 3>(7, 7) = -altro::skew(w);
  // dvbdot/domega
  dfc_dx.block<3, 3>(7, 10) = altro::skew(vb);

  // domegadot/domega
  dfc_dx.block<3, 3>(10, 10) =
      -moment_of_inertia_.inverse() * (altro::skew(w) * moment_of_inertia_ -
                                       altro::skew(moment_of_inertia_ * w));

  /////////////////////////////////////////////////////////////////////////////////////
  // Calculate dfc_du
  Eigen::MatrixXd dfc_du(13, 4);
  dfc_du.setZero();

  // dvw/du
  Eigen::Vector4d u_vec_sign;
  for (int i = 0; i < 4; ++i) {
    u_vec_sign(i) = u_vec(i) > 0 ? 1 : 0;
  }
  dfc_du.block<1, 4>(9, 0) = (1 / robot_mass * kf_) * u_vec_sign;

  // dw_dot/du
  Eigen::Matrix<double, 3, 4> Fmat = forceMatrix().block<3, 4>(1, 0);
  Fmat.block<2, 4>(0, 0) =
      Fmat.block<2, 4>(0, 0).array().rowwise() * u_vec_sign.transpose().array();

  dfc_du.block<3, 4>(10, 0) = moment_of_inertia_.inverse() * Fmat;
  //

  // Get Jacobian
  J.block<13, 13>(0, 0) = dfc_dx;
  J.block<13, 4>(0, 13) = dfc_du;
}

void QuadrupedQuaternionModel::Dynamics(
    double* x_dot,
    const double* x,
    const double* u,
    Eigen::Matrix<double, 3, 4> foot_pos_body,
    Eigen::Matrix3d inertia_body) const {
  Eigen::Map<Eigen::VectorXd> x_dot_vec(x_dot, 13);
  Eigen::Map<const Eigen::VectorXd> x_vec(x, 13);
  Eigen::Map<const Eigen::VectorXd> u_vec(u, 12);

  double robot_mass = 13;
  Eigen::Vector3d g_vec;
  g_vec << 0, 0, -9.81;

  Eigen::Vector3d moment_body;
  moment_body =
      altro::skew(foot_pos_body.block<3, 1>(0, 0)) * u_vec.segment<3>(0) +
      altro::skew(foot_pos_body.block<3, 1>(0, 1)) * u_vec.segment<3>(3) +
      altro::skew(foot_pos_body.block<3, 1>(0, 2)) * u_vec.segment<3>(6) +
      altro::skew(foot_pos_body.block<3, 1>(0, 3)) * u_vec.segment<3>(9);

  // change rate of position
  x_dot_vec.segment<3>(0) = x_vec.segment<3>(7);
  // change rate of quaternion
  x_dot_vec.segment<4>(3) =
      0.5 * altro::G(x_vec.segment<4>(3)) * x_vec.segment<3>(10);
  // change rate of linear velocity
  x_dot_vec.segment<3>(7) = (u_vec.segment<3>(0) + u_vec.segment<3>(3) +
                             u_vec.segment<3>(6) + u_vec.segment<3>(9)) /
                                robot_mass +
                            g_vec;
  // change rate of angular velocity
  x_dot_vec.segment<3>(10) =
      inertia_body.inverse() *
      (moment_body -
       x_vec.segment<3>(10).cross(inertia_body * x_vec.segment<3>(10)));
}

void QuadrupedQuaternionModel::Jacobian(
    double* jac,
    const double* x,
    const double* u,
    Eigen::Matrix<double, 3, 4> foot_pos_body,
    Eigen::Matrix3d inertia_body) const {
  (void)u;
  Eigen::Map<Eigen::Matrix<double, 13, 25>> J(jac);  // jac = [dfc_dx, dfc_du]
  J.setZero();

  Eigen::Map<const Eigen::VectorXd> x_vec(x, 13);

  double robot_mass = 13;

  // Calculate dfc_dx
  Eigen::MatrixXd dfc_dx(13, 13);
  dfc_dx.setZero();
  // dv/dv
  dfc_dx.block<3, 3>(0, 7) = Eigen::Matrix3d::Identity();
  // dqdot/dq
  dfc_dx.block<1, 3>(3, 4) = -0.5 * x_vec.segment<3>(10).transpose();
  dfc_dx.block<3, 1>(4, 3) = 0.5 * x_vec.segment<3>(10);
  dfc_dx.block<3, 3>(4, 4) = -0.5 * altro::skew(x_vec.segment<3>(10));
  // dqdot/domega
  dfc_dx(3, 10) = -0.5 * x_vec(4);  // -0.5qa
  dfc_dx(3, 11) = -0.5 * x_vec(5);  // -0.5qb
  dfc_dx(3, 12) = -0.5 * x_vec(6);  // -0.5qc
  dfc_dx(4, 10) = 0.5 * x_vec(3);   // 0.5qs
  dfc_dx(4, 11) = -0.5 * x_vec(6);  // -0.5qc
  dfc_dx(4, 12) = 0.5 * x_vec(5);   // 0.5qb
  dfc_dx(5, 10) = 0.5 * x_vec(6);   // 0.5qc
  dfc_dx(5, 11) = 0.5 * x_vec(3);   // 0.5qs
  dfc_dx(5, 12) = -0.5 * x_vec(4);  // -0.5qa
  dfc_dx(6, 10) = -0.5 * x_vec(5);  // -0.5qb
  dfc_dx(6, 11) = 0.5 * x_vec(4);   // 0.5qa
  dfc_dx(6, 12) = 0.5 * x_vec(3);   // 0.5qs
  // domegadot/domega
  dfc_dx.block<3, 3>(10, 10) =
      -inertia_body.inverse() *
      (altro::skew(x_vec.segment<3>(10)) * inertia_body -
       altro::skew(inertia_body * x_vec.segment<3>(10)));

  // Calculate dfc_du
  Eigen::MatrixXd dfc_du(13, 12);
  dfc_du.setZero();

  for (int i = 0; i < 4; ++i) {
    dfc_du.block<3, 3>(7, 3 * i) =
        (1 / robot_mass) * Eigen::Matrix3d::Identity();
    dfc_du.block<3, 3>(10, 3 * i) =
        inertia_body.inverse() * altro::skew(foot_pos_body.block<3, 1>(0, i));
  }

  // Get Jacobian
  J.block<13, 13>(0, 0) = dfc_dx;
  J.block<13, 12>(0, 13) = dfc_du;
}

Eigen::Vector4d Slerp(Eigen::Vector4d q1, Eigen::Vector4d q2, double t) {
  if (q1 == q2) {
    return q1;
  } else {
    double dot = q1.dot(q2);
    double theta = acos(dot);
    double sinTheta = sin(theta);
    double a = sin((1 - t) * theta) / sinTheta;
    double b = sin(t * theta) / sinTheta;

    return a * q1 + b * q2;
  }
}