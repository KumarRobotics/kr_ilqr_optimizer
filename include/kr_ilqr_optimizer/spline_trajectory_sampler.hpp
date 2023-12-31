#pragma once
#include <eigen_conversions/eigen_msg.h>
#include <kr_planning_msgs/SplineTrajectory.h>
#include <kr_planning_msgs/TrajectoryDiscretized.h>
#include <visualization_msgs/Marker.h>

#include "ros/ros.h"
// #include <kr_planning_rviz_plugins/spline_trajectory_visual.h>
#include <Eigen/Eigen>
#include <boost/bind.hpp>

//
#include <std_msgs/ColorRGBA.h>
#include <std_msgs/Float64.h>
#include <std_msgs/Header.h>
#include <std_msgs/Int32.h>

#include <fstream>
#include <iostream>
#include <kr_ilqr_optimizer/quadMPC.hpp>
#include <sstream>

class SplineTrajSampler {
 protected:
  const double dt_ = 0.1;

  bool compute_altro_ = true;
  bool write_summary_to_file_ = false;
  bool subscribe_to_traj_ = true;
  bool publish_optimized_traj_ = true;
  bool publish_viz_ = true;

  uint N_iter_ = 0;
  std::ofstream poly_array_file;

  ros::Subscriber sub_;
  ros::Publisher sampled_traj_pub_;
  ros::Publisher opt_traj_pub_;
  ros::Publisher opt_viz_pub_;

  double g_;
  double mass_;
  Eigen::Matrix3d inertia_;
  double min_thrust_;  // in Newtons
  double max_thrust_;  // in Newtons

  ros::NodeHandle nhp_;

 public:
  std::unique_ptr<quadMPC> mpc_solver;

  using Ptr = std::unique_ptr<SplineTrajSampler>;

  std::vector<Eigen::VectorXd> opt_traj_;  // empty if no solution yet

  SplineTrajSampler(ros::NodeHandle& nh,
                    bool subscribe_to_traj_,
                    bool publish_optimized_traj_,
                    bool publish_viz_,
                    double dt)
      : subscribe_to_traj_(subscribe_to_traj_),
        publish_optimized_traj_(publish_optimized_traj_),
        publish_viz_(publish_viz_),
        dt_(dt) {
    nhp_ = nh;
    bool use_quat = false;
    mpc_solver = std::make_unique<quadMPC>(dt, use_quat);

    nhp_.param("VehicleMass", mass_, 1.5);
    nhp_.param("GravAcc", g_, 9.81);
    double Ixx, Iyy, Izz;
    nhp_.param("Ixx", Ixx, 1.0);
    nhp_.param("Ixx", Iyy, 1.0);
    nhp_.param("Ixx", Izz, 1.0);
    inertia_ << Ixx, 0, 0, 0, Iyy, 0, 0, 0, Izz;
    nhp_.param("MinThrust", min_thrust_, 0.0);
    nhp_.param("MaxThrust", max_thrust_, 8.0 * 4);  // still in force N here
    double arm_length, kf, km;
    nhp_.param("arm_length", arm_length, 0.4);
    nhp_.param("kf", kf, 1.0);
    nhp_.param("km", km, 0.0245);

    // mass_ = 1.5;
    // kf = 1.0;
    // km = 0.0245;
    // arm_length = 0.4;
    // inertia_ = Eigen::Matrix3d::Identity();
    // min_thrust_ = 0.0;
    // max_thrust_ = 8.0*4;

    if (publish_optimized_traj_) {
      opt_traj_pub_ = nhp_.advertise<kr_planning_msgs::TrajectoryDiscretized>(
          "optimized_traj_samples", 1);
    }
    if (publish_viz_) {
      opt_viz_pub_ = nhp_.advertise<visualization_msgs::Marker>(
          "optimized_traj_samples_viz", 1);
    }
    if (subscribe_to_traj_) {
      sub_ = nhp_.subscribe(
          "local_plan_server/trajectory_planner/search_trajectory",
          1,
          &SplineTrajSampler::callbackWrapper,
          this);
    }
    sampled_traj_pub_ = nhp_.advertise<kr_planning_msgs::TrajectoryDiscretized>(
        "spline_traj_samples", 1);
    // opt_traj_pub_ =
    // n.advertise<kr_planning_msgs::TrajectoryDiscretized>("optimized_traj_samples",
    // 1); opt_viz_pub_ =
    // n.advertise<visualization_msgs::Marker>("optimized_traj_samples_viz", 1);
    // sub_ =
    // n.subscribe("/local_plan_server/trajectory_planner/search_trajectory", 1,
    // &SplineTrajSampler::callbackWrapper, this);
    if (write_summary_to_file_) {
      poly_array_file.open("/home/yifei/planning_summary.csv",
                           std::ios_base::app);
      poly_array_file << "Planning Iteration, Starting x, Jerk Norm Sum, Traj "
                         "Time, Total_Ref_Fdt, Total_Ref_Mdt \n";
    }
    mpc_solver->SetUp(
        mass_, g_, inertia_, min_thrust_, max_thrust_, kf, km, arm_length);
  }
  void callbackWrapper(const kr_planning_msgs::SplineTrajectory::ConstPtr& msg);

  kr_planning_msgs::TrajectoryDiscretized sample_and_refine_trajectory(
      const Eigen::VectorXd& start_state,
      const kr_planning_msgs::SplineTrajectory& traj,
      const std::vector<Eigen::MatrixXd>& hPolys,
      const Eigen::VectorXd& allo_ts);

 private:
  kr_planning_msgs::TrajectoryDiscretized publish_altro(
      const Eigen::VectorXd& start_state,
      std::vector<Eigen::Vector3d> pos,
      std::vector<Eigen::Vector3d> vel,
      std::vector<Eigen::Vector3d> acc,
      std::vector<double> yaw_ref,
      std::vector<double> thrust,
      std::vector<Eigen::Vector3d> moment,
      double t,
      std::vector<Eigen::Quaterniond> q_ref,
      std::vector<Eigen::Vector3d> w_ref,
      const std::vector<std::pair<Eigen::MatrixXd, Eigen::VectorXd>>& h_polys,
      const Eigen::VectorXd& allo_ts,
      std_msgs::Header header);
  Eigen::Vector3d compute_ref_inputs(const Eigen::Vector3d& pos,
                                     const Eigen::Vector3d& vel,
                                     const Eigen::Vector3d& acc,
                                     const Eigen::Vector3d& jerk,
                                     const Eigen::Vector3d& snap,
                                     const Eigen::Vector3d& yaw_dyaw_ddyaw,
                                     double& thrust,
                                     geometry_msgs::Point& moment,
                                     Eigen::Quaterniond& q_return,
                                     Eigen::Vector3d& omega);

  void unpack(geometry_msgs::Point& p,
              const Eigen::VectorXd& v,
              int start_idx = 0) {
    p.x = v(start_idx);
    p.y = v(start_idx + 1);
    p.z = v(start_idx + 2);
  }

  std::vector<double> linspace(double start, double end, int numPoints) {
    std::vector<double> result;
    if (numPoints <= 1) {
      std::cerr << "linspace: numPoints must be > 1\n";
    }

    double step = (end - start) / (numPoints - 1);
    for (int i = 0; i < numPoints; ++i) {
      result.push_back(start + i * step);
    }

    return result;
  }
  Eigen::Matrix3d skew(const Eigen::Vector3d& x) {
    Eigen::Matrix3d s;
    s << 0, -x[2], x[1], x[2], 0, -x[0], -x[1], x[0], 0;
    return s;
  }

  /**
   * This tutorial demonstrates simple sending of messages over the ROS system.
   */
  std::vector<float> differentiate(const std::vector<float>& p,
                                   float segment_time) {
    if (p.size() < 2) return std::vector<float>();
    std::vector<float> v;
    for (int i = 1; i < p.size(); i++) {
      v.push_back(p[i] * static_cast<float>(i) / segment_time);
    }
    return v;
  }
  Eigen::VectorXd evaluate(const kr_planning_msgs::SplineTrajectory& msg,
                           double t,
                           uint deriv_num) {
    Eigen::VectorXd result(msg.dimensions);

    for (int dim = 0; dim < msg.dimensions; dim++) {
      auto spline = msg.data[dim];
      double time_elapsed = 0;
      for (auto poly : spline.segs) {
        auto poly_coeffs = poly.coeffs;
        for (int d = 0; d < deriv_num; d++) {
          poly_coeffs = differentiate(poly_coeffs, poly.dt);
        }
        if (poly_coeffs.size() == 0) {
          result(dim) = 0;
          break;
        }
        result(dim) = poly_coeffs[0];

        if (t < time_elapsed + poly.dt || poly == spline.segs.back()) {
          for (int j = 1; j < poly_coeffs.size(); j++) {
            result(dim) +=
                poly_coeffs[j] * std::pow((t - time_elapsed) / poly.dt, j);
          }
          break;
        }
        time_elapsed += poly.dt;
      }
    }
    return result;
  }

  std::vector<Eigen::Vector3d> sample(
      const kr_planning_msgs::SplineTrajectory& msg,
      double total_time,
      int N_state,
      int deriv_num) {
    std::vector<Eigen::Vector3d> ps(N_state);  // 101
    double dt = total_time / (N_state - 1);    // dt = 1 / 100 = 0.01
    for (int i = 0; i < N_state; i++)          // 0, 1, 2, ..., 100 * 0.01 = 1
      ps.at(i) = evaluate(msg, i * dt, deriv_num);
    return ps;
  }
};
