#include <iostream>

#include "kr_ilqr_optimizer/spline_trajectory_sampler.hpp"
#include "ros/ros.h"
int main(int argc, char** argv) {
  /**
   * The ros::init() function needs to see argc and argv so that it can perform
   * any ROS arguments and name remapping that were provided at the command
   * line. For programmatic remappings you can use a different version of init()
   * which takes remappings directly, but for most command-line programs,
   * passing argc and argv is the easiest way to do it.  The third argument to
   * init() is the name of the node.
   *
   * You must call one of the versions of ros::init() before using any other
   * part of the ROS system.
   */
  ros::init(argc, argv, "spline_traj_sampler");

  /**
   * NodeHandle is the main access point to communications with the ROS system.
   * The first NodeHandle constructed will fully initialize this node, and the
   * last NodeHandle destructed will close down the node.
   */
  ros::NodeHandle nh;
  // std::shared_ptr<SplineTrajectoryVisual> visual_;
  SplineTrajSampler sampler(nh, true, true, true, 40);
  // sampler.mpc_solver->SetUp();
  std::cout << "sampler setup Complete" << std::endl;
  ros::spin();
  return 0;
}