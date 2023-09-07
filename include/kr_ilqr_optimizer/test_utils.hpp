//
// Created by Brian Jackson on 9/27/22.
// Copyright (c) 2022 Robotic Exploration Lab. All rights reserved.
//

#pragma once

#include <vector>

#include "Eigen/Dense"
#include "altro/typedefs.hpp"
#include "kr_ilqr_optimizer/quaternion_utils.hpp"

void discrete_double_integrator_dynamics(
    double* xnext, const double* x, const double* u, float h, int dim);

void discrete_double_integrator_jacobian(
    double* jac, const double* x, const double* u, float h, int dim);

void cartpole_dynamics_midpoint(double* xnext,
                                const double* x,
                                const double* u,
                                float h);
void cartpole_jacobian_midpoint(double* xnext,
                                const double* x,
                                const double* u,
                                float h);

void pendulum_dynamics(double* xnext, const double* x, const double* u);
void pendulum_jacobian(double* jac, const double* x, const double* u);

using ContinuousDynamicsFunction =
    std::function<void(double*, const double*, const double*)>;
using ContinuousDynamicsJacobian =
    std::function<void(double*, const double*, const double*)>;

altro::ExplicitDynamicsFunction MidpointDynamics(int n,
                                                 int m,
                                                 ContinuousDynamicsFunction f);
altro::ExplicitDynamicsJacobian MidpointJacobian(
    int n, int m, ContinuousDynamicsFunction f, ContinuousDynamicsJacobian jac);

altro::ExplicitDynamicsFunction ForwardEulerDynamics(
    int n, int m, const ContinuousDynamicsFunction f);
altro::ExplicitDynamicsJacobian ForwardEulerJacobian(
    int n,
    int m,
    const ContinuousDynamicsFunction f,
    const ContinuousDynamicsJacobian df);

class BicycleModel {
 public:
  enum class ReferenceFrame { CenterOfGravity, Rear, Front };

  explicit BicycleModel(ReferenceFrame frame = ReferenceFrame::CenterOfGravity)
      : reference_frame_(frame) {}

  void Dynamics(double* x_dot, const double* x, const double* u) const;
  void Jacobian(double* jac, const double* x, const double* u) const;

  void SetLengths(double length, double length_to_rear_wheel_from_cg) {
    length_ = length;
    distance_to_rear_wheels_ = length_to_rear_wheel_from_cg;
  }

  static constexpr int NumStates = 4;
  static constexpr int NumInputs = 2;

 private:
  ReferenceFrame reference_frame_;
  double length_ = 2.7;
  double distance_to_rear_wheels_ = 1.5;
};

class quadModel {
 public:
  enum class ReferenceFrame { CenterOfGravity, Rear, Front };

  explicit quadModel(double mass,
                     double grav_const,
                     Eigen::Matrix3d inertia,
                     double min_w_sq,
                     double max_w_sq,
                     double kf,
                     double km,
                     double L,
                     ReferenceFrame frame = ReferenceFrame::CenterOfGravity)
      : reference_frame_(frame),
        mass_(mass),
        grav_(grav_const),
        moment_of_inertia_(inertia),
        min_w_sq(min_w_sq),
        max_w_sq(max_w_sq),
        kf_(kf),
        km_(km),
        motor_dist_(L) {}

  void Dynamics(double* x_dot, const double* x, const double* u) const;
  void Jacobian(double* jac, const double* x, const double* u) const;
  void Jacobian_fd(double* jac, const double* x, const double* u) const;
  Eigen::Vector3d moments(const Eigen::VectorXd& u) const;
  Eigen::Vector3d forces(const Eigen::VectorXd& u) const;
  // Eigen::VectorXd f_quad(const Eigen::VectorXd& x_vec, const Eigen::VectorXd&
  // u_vec) const; void finite_jacobian(const Eigen::Ref<const Eigen::VectorXd>&
  // x,
  //   const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& f,
  //   Eigen::MatrixXd& jac,
  //   const double eps);
  void finite_jacobian_quad_xu(double* jac,
                               const double* x,
                               const double* u) const;

  static constexpr int NumStates = 13;
  static constexpr int NumInputs = 4;
  double get_hover_input() const { return mass_ * grav_ / 4.0 / kf_; }
  float max_w_sq;
  float min_w_sq;
  double kf_;
  double mass_;
  double grav_;
  Eigen::Matrix3d moment_of_inertia_ = Eigen::Matrix3d::Identity();
  Eigen::Matrix4d forceMatrix_inv() const { return forceMatrix().inverse(); }
  Eigen::Matrix4d forceMatrix() const {
    double L = motor_dist_;

    Eigen::Matrix4d SMatrix;
    SMatrix << kf_, kf_, kf_, kf_, 0, L * kf_, 0, -L * kf_, -L * kf_, 0,
        L * kf_, 0, km_, -km_, km_, -km_;

    return SMatrix;
  }

 private:
  ReferenceFrame reference_frame_;
  double motor_dist_ = 0.4;
  double km_ = 0.0245;
};

class SimpleQuaternionModel {
 public:
  void Dynamics(double* x_dot, const double* x, const double* u) const;
  void Jacobian(double* jac, const double* x, const double* u) const;

  static constexpr int NumStates = 4;  // Quaternion: [qs, qa, qb, qc]
  static constexpr int NumInputs =
      3;  // Angular Velocity: [omega_x, omega_y, omega_z]

  static constexpr int NumErrorStates = 3;
  static constexpr int NumErrorInputs = 3;
};

class QuadrupedQuaternionModel {
 public:
  void Dynamics(double* x_dot,
                const double* x,
                const double* u,
                Eigen::Matrix<double, 3, 4> foot_pos_body,
                Eigen::Matrix3d inertia_body) const;
  void Jacobian(double* jac,
                const double* x,
                const double* u,
                Eigen::Matrix<double, 3, 4> foot_pos_body,
                Eigen::Matrix3d inertia_body) const;

  static constexpr int NumStates = 13;  // r, q, v, w
  static constexpr int NumInputs = 12;  // f1, f2, f3, f4

  static constexpr int NumErrorStates = 12;
  static constexpr int NumErrorInputs = 12;
};

// void ReadScottyTrajectory(int *Nref, float *tref,
// std::vector<Eigen::Vector4d> *xref,
//                           std::vector<Eigen::Vector2d> *uref);

Eigen::Vector4d Slerp(Eigen::Vector4d quat1, Eigen::Vector4d quat2, double t);