#include "ros/ros.h"
#include <kr_planning_msgs/SplineTrajectory.h>
#include <kr_planning_msgs/TrajectoryDiscretized.h>
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

#include <sstream>
#include "altro/altro.hpp"
#include "test_utils.hpp"
#include "finitediff.hpp"
#include "quadMPC.hpp"
class SplineTrajSampler
{
  public:
  SplineTrajSampler(){
    ros::NodeHandle n;
    pub_ = n.advertise<kr_planning_msgs::TrajectoryDiscretized>("spline_traj_samples", 1);
    // sub_ = n.subscribe("/local_plan_server/trajectory", 1, &SplineTrajSampler::callbackWrapper, this);
    sub_ = n.subscribe("/local_plan_server/trajectory_planner/search_trajectory", 1, &SplineTrajSampler::callbackWrapper, this);
    if (write_summary_to_file_){
      poly_array_file.open("/home/yifei/planning_summary.csv",std::ios_base::app);
      poly_array_file << "Planning Iteration, Starting x, Jerk Norm Sum, Traj Time, Total_Ref_Fdt, Total_Ref_Mdt \n";
    }
  
  }
  void callbackWrapper(const kr_planning_msgs::SplineTrajectory::ConstPtr& msg)
{
    kr_planning_msgs::SplineTrajectory traj= *msg;
    int N_sample_pts = 100;
    double total_time1 = 0.0;
    double total_time2 = 0.0;
    double total_time3 = 0.0;
    double total_time4 = 0.0;
    double total_time5 = 0.0;
    auto pos = sample(traj, N_sample_pts, 0, total_time1);
    auto vel = sample(traj, N_sample_pts, 1, total_time2);
    auto acc = sample(traj, N_sample_pts, 2, total_time3);
    auto jerk = sample(traj, N_sample_pts, 3, total_time4);
    auto snap = sample(traj, N_sample_pts, 4, total_time5);
    double dt = total_time1/N_sample_pts;

    kr_planning_msgs::TrajectoryDiscretized traj_discretized; 
    traj_discretized.header = traj.header;
    double jerk_abs_sum = 0.0;
    double total_ref_F = 0.0;
    double total_ref_M = 0.0;
    Eigen::Quaterniond initial_q;
    Eigen::Vector3d inital_w;
    geometry_msgs::Quaternion ini_q_msg;
    geometry_msgs::Point ini_w_msg;

    for (int i = 0; i < N_sample_pts; i++) {

      geometry_msgs::Point pos_t;
      geometry_msgs::Point vel_t;
      geometry_msgs::Point acc_t;
      geometry_msgs::Point jerk_t;
      geometry_msgs::Point snap_t;

      double thrust;
      geometry_msgs::Point moment;

      Eigen::Vector3d yaw_three_der = Eigen::Vector3d(0, 0, 0);
      Eigen::Vector3d M = compute_ref_inputs(pos[i], vel[i], acc[i], jerk[i], snap[i], yaw_three_der, thrust, moment, initial_q, inital_w);
      
      if (i == 0){
        tf::quaternionEigenToMsg	(	initial_q, ini_q_msg);
        tf::pointEigenToMsg	(	inital_w, ini_w_msg);
      }

      jerk_abs_sum += acc[i].squaredNorm()*dt;
      total_ref_F += abs(thrust)*dt;
      total_ref_M += M.norm()*dt;


      unpack(pos_t, pos[i]);
      unpack(vel_t, vel[i]);
      unpack(acc_t, acc[i]);
      // unpack(jerk_t, jerk[i]);
      // unpack(snap_t, snap[i]);
      
     

      // Assuming the Eigen vector contains 3D points (x, y, z)

      traj_discretized.pos.push_back(pos_t);
      traj_discretized.vel.push_back(vel_t);
      traj_discretized.acc.push_back(acc_t);
      // traj_discretized.jerk.push_back(jerk_t);
      // traj_discretized.snap.push_back(snap_t);
      traj_discretized.thrust.push_back(thrust);
      traj_discretized.moment.push_back(moment);

      // std::cout << thrust << std::endl;
    }
    if (write_summary_to_file_){
      //get start location as key to which ones successed
      poly_array_file << std::to_string(N_iter_)+ ", " + std::to_string(pos[0](0)) + ", " + std::to_string(jerk_abs_sum) + ", " + std::to_string(total_time1) +"," + std::to_string(total_ref_F)+","+std::to_string(total_ref_M)+"\n";
      //compute the total snap at 100 points
    }
    traj_discretized.t =  linspace(0.0, total_time1, N_sample_pts);
    traj_discretized.N = N_sample_pts;
    traj_discretized.inital_attitude = ini_q_msg;
    traj_discretized.initial_omega = ini_w_msg;
    ros::Duration(0.5).sleep(); // this is for yuwei's planner to finish so we have polytope info
    this->pub_.publish(traj_discretized);
    N_iter_++;
  }

  protected:
  bool write_summary_to_file_ = false;
  uint N_iter_ = 0;
  std::ofstream poly_array_file;
  
  ros::Subscriber sub_;
  ros::Publisher pub_;
  double time_limit_ = 5.0;
  double g_ = 9.81;
  double mass_ = 0.5;
  Eigen::DiagonalMatrix<double, 3> inertia_ = Eigen::DiagonalMatrix<double, 3> (0.0023, 0.0023, 0.004);

  private:

  Eigen::Vector3d compute_ref_inputs(Eigen::Vector3d pos, Eigen::Vector3d vel, Eigen::Vector3d acc, Eigen::Vector3d jerk, Eigen::Vector3d snap, 
  Eigen::Vector3d yaw_dyaw_ddyaw, double& thrust, geometry_msgs::Point& moment, Eigen::Quaterniond& initial_q, Eigen::Vector3d& initial_w){

    // Desired force vector.
    Eigen::Vector3d t = acc + Eigen::Vector3d(0, 0, g_);
    Eigen::Vector3d b3 = t.normalized();
    Eigen::Vector3d F_des = mass_ * t;

    // Desired thrust is force projected onto b3 axis.
    double u1 = F_des.dot(b3);

    // Desired orientation to obtain force vector.
    Eigen::Vector3d b3_des = F_des.normalized();
    double yaw_des = yaw_dyaw_ddyaw[0];
    Eigen::Vector3d c1_des(std::cos(yaw_des), std::sin(yaw_des), 0);
    Eigen::Vector3d b2_des = (b3_des.cross(c1_des)).normalized();
    Eigen::Vector3d b1_des = b2_des.cross(b3_des);
    Eigen::Matrix3d R_des;
    R_des << b1_des, b2_des, b3_des;

    initial_q = R_des;

    Eigen::Matrix3d R = R_des; // assume we have perfect tracking on rotation

    // Following section follows Mellinger paper to compute reference angular velocity
    double dot_u1 = b3.dot(jerk);
    Eigen::Vector3d hw = mass_ / u1 * (jerk - dot_u1 * b3);
    double p = -hw.dot(b2_des);
    double q = hw.dot(b1_des);
    Eigen::Vector3d w_des(0, 0, yaw_dyaw_ddyaw[1]);
    double r = w_des.dot(b3_des);
    Eigen::Vector3d Omega(p, q, r);

    initial_w = Omega;

    Eigen::Vector3d wwu1b3 = Omega.cross(Omega.cross(u1 * b3) );
    double ddot_u1 = b3.dot(mass_ * snap) - b3.dot(wwu1b3);
    Eigen::Vector3d ha = 1.0 / u1 * (mass_ * snap - ddot_u1 * b3 - 2 * Omega.cross( dot_u1 * b3) - wwu1b3);
    double p_dot = -ha.dot(b2_des);
    double q_dot = ha.dot(b1_des);
    double r_dot = yaw_dyaw_ddyaw[2]* b3_des.dot(Eigen::Vector3d(0, 0, 1.0));
    Eigen::Vector3d Alpha(p_dot, q_dot, r_dot);

    Eigen::Vector3d u2 = inertia_ * Alpha + Omega.cross( inertia_ * Omega);
    Eigen::Vector3d TM;
    unpack(moment, u2);
    thrust = u1;
    return u2;
  }

  void unpack(geometry_msgs::Point& p, const Eigen::VectorXd& v) {
    p.x = v(0);
    p.y = v(1);
    p.z = v(2);
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


    std::vector<Eigen::VectorXd> sample(
        const kr_planning_msgs::SplineTrajectory& msg, int N, int deriv_num, double& total_time) {
      std::vector<Eigen::VectorXd> ps(N + 1);
      total_time = msg.data[0].t_total;
      // std::cout<<"total_time: "<<total_time<<std::endl;
      // if (total_time > time_limit_) total_time = time_limit_; //for a max of 6 seconds, too long make ilqr sovler fail
      double dt = total_time / N;
      for (int i = 0; i <= N; i++) ps.at(i) = evaluate(msg, i * dt, deriv_num);

      return ps;
    }

};


int main(int argc, char **argv)
{
  /**
   * The ros::init() function needs to see argc and argv so that it can perform
   * any ROS arguments and name remapping that were provided at the command line.
   * For programmatic remappings you can use a different version of init() which takes
   * remappings directly, but for most command-line programs, passing argc and argv is
   * the easiest way to do it.  The third argument to init() is the name of the node.
   *
   * You must call one of the versions of ros::init() before using any other
   * part of the ROS system.
   */
  ros::init(argc, argv, "spline_traj_sampler");

  /**
   * NodeHandle is the main access point to communications with the ROS system.
   * The first NodeHandle constructed will fully initialize this node, and the last
   * NodeHandle destructed will close down the node.
   */

  // std::shared_ptr<SplineTrajectoryVisual> visual_;
  SplineTrajSampler sampler;
  altro::ALTROSolver solver(10);
  solver.SetDimension(13, 4, 0, altro::LastIndex);
  fmt::print("Solver Initialized!\n");

  quadMPC mpc_solver;
  mpc_solver.SetUp();

  
    fmt::print("Setup success\n");

  mpc_solver.eg1();

  ros::spin();
  return 0;
}