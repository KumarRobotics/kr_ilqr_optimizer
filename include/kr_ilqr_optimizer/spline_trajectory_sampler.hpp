#pragma once
#include "ros/ros.h"
#include <kr_planning_msgs/SplineTrajectory.h>
#include <kr_planning_msgs/TrajectoryDiscretized.h>
#include <visualization_msgs/Marker.h>
#include <eigen_conversions/eigen_msg.h>
// #include <kr_planning_rviz_plugins/spline_trajectory_visual.h>
#include <boost/bind.hpp>
#include <Eigen/Eigen>

//
#include <iostream>
#include <fstream>

#include "ros/ros.h"
#include <std_msgs/Float64.h>
#include <std_msgs/Int32.h>
#include <std_msgs/Header.h>
#include <std_msgs/ColorRGBA.h>

#include <sstream>
#include "kr_ilqr_optimizer/quadMPC.hpp"


class SplineTrajSampler
{
  protected:
  int N_sample_pts_ = 80;
  double time_limit_ = 8;

  std::unique_ptr<quadMPC> mpc_solver;
  bool compute_altro_ = true;
  bool write_summary_to_file_ = false;
  uint N_iter_ = 0;
  std::ofstream poly_array_file;
  
  ros::Subscriber sub_;
  ros::Publisher sampld_traj_pub_;
  ros::Publisher opt_traj_pub_;
  ros::Publisher opt_viz_pub_;
  double g_ = 9.81;
  double mass_ = 0.5;
  Eigen::DiagonalMatrix<double, 3> inertia_ = Eigen::DiagonalMatrix<double, 3> (0.0023, 0.0023, 0.004);

  public:
    using Ptr = std::unique_ptr<SplineTrajSampler>;

  SplineTrajSampler(){
    bool use_quat = false;
    mpc_solver = std::make_unique<quadMPC>(N_sample_pts_, time_limit_, use_quat);
    ros::NodeHandle n;
    sampld_traj_pub_ = n.advertise<kr_planning_msgs::TrajectoryDiscretized>("spline_traj_samples", 1);
    opt_traj_pub_ = n.advertise<kr_planning_msgs::TrajectoryDiscretized>("optimized_traj_samples", 1);
    opt_viz_pub_ = n.advertise<visualization_msgs::Marker>("optimized_traj_samples_viz", 1);
    // sub_ = n.subscribe("/local_plan_server/trajectory", 1, &SplineTrajSampler::callbackWrapper, this);
    sub_ = n.subscribe("/local_plan_server/trajectory_planner/search_trajectory", 1, &SplineTrajSampler::callbackWrapper, this);
    if (write_summary_to_file_){
      poly_array_file.open("/home/yifei/planning_summary.csv",std::ios_base::app);
      poly_array_file << "Planning Iteration, Starting x, Jerk Norm Sum, Traj Time, Total_Ref_Fdt, Total_Ref_Mdt \n";
    }
    mpc_solver->SetUp();
  
  }
  void callbackWrapper(const kr_planning_msgs::SplineTrajectory::ConstPtr& msg);

  

  private:

  void publish_altro(std::vector<Eigen::Vector3d> pos, std::vector<Eigen::Vector3d> vel, std::vector<Eigen::Vector3d> acc, std::vector<double> yaw_ref,
                                std::vector<double> thrust, std::vector<Eigen::Vector3d> moment, std::vector<double> t, int N, std::vector<Eigen::Quaterniond> q_ref, std::vector<Eigen::Vector3d> w_ref, std_msgs::Header header);
  Eigen::Vector3d compute_ref_inputs(Eigen::Vector3d pos, Eigen::Vector3d vel, Eigen::Vector3d acc, Eigen::Vector3d jerk, Eigen::Vector3d snap, 
        Eigen::Vector3d yaw_dyaw_ddyaw, double& thrust, geometry_msgs::Point& moment, Eigen::Quaterniond& q_return, Eigen::Vector3d& initial_w);

  void unpack(geometry_msgs::Point& p, const Eigen::VectorXd& v, int start_idx = 0) {
    p.x = v(start_idx);
    p.y = v(start_idx + 1);
    p.z = v(start_idx + 2);
  }

  std::vector<double> linspace(double start, double end, int numPoints) {
    std::vector<double> result;
    if (numPoints <= 1) {
        std::cerr<<"linspace: numPoints must be > 1\n";
    }

    double step = (end - start) / (numPoints - 1);
    for (int i = 0; i < numPoints; ++i) {
        result.push_back(start + i * step);
    }

    return result;
  }

    /**
     * This tutorial demonstrates simple sending of messages over the ROS system.
     */
    std::vector<float> differentiate(const std::vector<float>& p, float segment_time)  {
      if (p.size() < 2) return std::vector<float>();
      std::vector<float> v;
      for (int i = 1; i < p.size(); i++) {
        v.push_back(p[i] * static_cast<float>(i) / segment_time);
      }
      return v;
    }
    Eigen::VectorXd evaluate(
        const kr_planning_msgs::SplineTrajectory& msg,
        double t,
        uint deriv_num) {
      Eigen::VectorXd result(msg.dimensions);

      for (int dim = 0; dim < msg.dimensions; dim++) {
        auto spline = msg.data[dim];
        double dt = 0;
        for (auto poly : spline.segs) {
          auto poly_coeffs = poly.coeffs;
          for (int d = 0; d < deriv_num; d++) {
            poly_coeffs = differentiate(poly_coeffs, poly.dt);
          }
          result(dim) = poly_coeffs[0];

          if (t < dt + poly.dt || poly == spline.segs.back()) {
            for (int j = 1; j < poly_coeffs.size(); j++) {
              result(dim) += poly_coeffs[j] * std::pow((t - dt) / poly.dt, j);
            }
            break;
          }
          dt += poly.dt;
        }
      }
      return result;
    }


    std::vector<Eigen::Vector3d> sample(
        const kr_planning_msgs::SplineTrajectory& msg, int N, int deriv_num) {
      std::vector<Eigen::Vector3d> ps(N + 1);
      // std::cout<<"total_time: "<<total_time<<std::endl;
      // if (total_time > time_limit_) total_time = time_limit_; //for a max of 6 seconds, too long make ilqr sovler fail
      double dt = time_limit_ / N;
      for (int i = 0; i <= N; i++) ps.at(i) = evaluate(msg, i * dt, deriv_num);

      return ps;
    }

};

