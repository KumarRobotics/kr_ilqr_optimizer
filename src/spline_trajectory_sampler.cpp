#include "kr_ilqr_optimizer/spline_trajectory_sampler.hpp"

#include <algorithm>

#include "altro/altro.hpp"
#include "kr_ilqr_optimizer/finitediff.hpp"
#include "kr_ilqr_optimizer/test_utils.hpp"

// #include "spline_trajectory_sampler.hpp"
using Vector = Eigen::Matrix<a_float, Eigen::Dynamic, 1>;

// helpers
// Function Converted from Julia Code
std::pair<Eigen::MatrixXd, Eigen::VectorXd> pt_normal_to_halfspace(
    const Eigen::MatrixXd& hPoly) {
  int size_n = hPoly.cols();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(size_n, 3);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(size_n);

  for (int ii = 0; ii < size_n; ++ii) {
    Eigen::Vector3d p = hPoly.block<3, 1>(0, ii);
    Eigen::Vector3d n = hPoly.block<3, 1>(3, ii);
    double c = n.dot(p);

    A.row(ii) = n;
    b(ii) = c;
  }

  return std::make_pair(A, b);
}

Eigen::Vector3d SplineTrajSampler::compute_ref_inputs(
    const Eigen::Vector3d& pos,
    const Eigen::Vector3d& vel,
    const Eigen::Vector3d& acc,
    const Eigen::Vector3d& jerk,
    const Eigen::Vector3d& snap,
    const Eigen::Vector3d& yaw_dyaw_ddyaw,
    double& thrust,
    geometry_msgs::Point& moment,
    Eigen::Quaterniond& q_return,
    Eigen::Vector3d& omega) {
  // std::cout << "Computing ref inputs" << std::endl;
  // std::cout << mpc_solver->model_ptr->mass_ << std::endl;
  // Eigen::Matrix3d inertia_;
  // double mass_;
  // double g_;  // TODO: Bad place to define this
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

  q_return = R_des;

  Eigen::Matrix3d R = R_des;  // assume we have perfect tracking on rotation

  double dot_u1 = mass_ * b3.dot(jerk);
  Eigen::Vector3d hw =
      mass_ / u1 *
      (jerk);  // removing  - dot_u1 * b3 since thomas' does not have it
  double p = -hw.dot(b2_des);
  double q = hw.dot(b1_des);

  Eigen::Vector3d e3(0, 0, 1);

  double r = ((1 - std::pow(e3.dot(b1_des), 2)) * yaw_dyaw_ddyaw[1] -
              e3.dot(b2_des) * q) /
             e3.dot(b3_des);
  Eigen::Vector3d Omega(p, q, r);

  Eigen::Matrix<double, 2, 3> mat;
  mat << -b2_des, b1_des;
  Eigen::Vector2d pq_dot =
      (mass_ / u1 * mat * snap - 2 * dot_u1 / u1 * Eigen::Vector2d(p, q) +
       r * Eigen::Vector2d(q, -p));

  Eigen::Matrix3d b_dot = R * skew(Omega);
  Eigen::Vector3d b1_dot = b_dot.col(0);
  Eigen::Vector3d b2_dot = b_dot.col(1);
  Eigen::Vector3d b3_dot = b_dot.col(2);

  double r_dot =
      -(e3.dot(b3_dot) * r + e3.dot(b2_dot) * q + e3.dot(b2_des) * pq_dot[1] +
        2 * e3.dot(b1_des) * e3.dot(b1_dot) * yaw_dyaw_ddyaw[1] +
        (std::pow(e3.dot(b1_des), 2) - 1) * yaw_dyaw_ddyaw[2]) /
      e3.dot(b3_des);

  Eigen::Vector3d Alpha(pq_dot[0], pq_dot[1], r_dot);
  Eigen::Vector3d u2 = inertia_ * Alpha + Omega.cross(inertia_ * Omega);

  // OLD from mellinger
  //  Eigen::Vector3d w_des(0, 0, yaw_dyaw_ddyaw[1]);
  //  double r = w_des.dot(b3_des);
  //  Eigen::Vector3d Omega(p, q, r);

  omega = Omega;

  // Eigen::Vector3d wwu1b3 = Omega.cross(Omega.cross(u1 * b3));
  // double ddot_u1 = b3.dot(mass_ * snap) - b3.dot(wwu1b3);
  // Eigen::Vector3d ha =
  //     1.0 / u1 *
  //     (mass_ * snap - ddot_u1 * b3 - 2 * Omega.cross(dot_u1 * b3) - wwu1b3);
  // double p_dot = -ha.dot(b2_des);
  // double q_dot = ha.dot(b1_des);
  // double r_dot = yaw_dyaw_ddyaw[2] * b3_des.dot(Eigen::Vector3d(0, 0, 1.0));
  // Eigen::Vector3d Alpha(p_dot, q_dot, r_dot);

  Eigen::Vector3d TM;
  unpack(moment, u2);
  thrust = u1;
  return u2;
}

kr_planning_msgs::TrajectoryDiscretized SplineTrajSampler::publish_altro(
    const Eigen::VectorXd& start_state,
    std::vector<Eigen::Vector3d> pos,
    std::vector<Eigen::Vector3d> vel,
    std::vector<Eigen::Vector3d> acc,
    std::vector<double> yaw_ref,
    std::vector<double> thrust,
    std::vector<Eigen::Vector3d> moment,
    double dt,
    std::vector<Eigen::Quaterniond> q_ref,
    std::vector<Eigen::Vector3d> w_ref,
    const std::vector<std::pair<Eigen::MatrixXd, Eigen::VectorXd>>& h_polys,
    const Eigen::VectorXd& allo_ts,
    std_msgs::Header header) {
  std::vector<Vector> X_sim;  // for return
  std::vector<Vector> U_sim;  // for return
  std::vector<double> t_sim;  // for return
  ROS_WARN("[iLQR Optimizer]: Solving");
  int N_controls = pos.size() - 1;
  uint solver_status = mpc_solver->solve_problem(start_state,
                                                 pos,
                                                 vel,
                                                 acc,
                                                 yaw_ref,
                                                 thrust,
                                                 moment,
                                                 dt,
                                                 q_ref,
                                                 w_ref,
                                                 h_polys,
                                                 allo_ts,
                                                 X_sim,
                                                 U_sim,
                                                 t_sim);
  // DEAL WITH VIZ
  visualization_msgs::Marker viz_msg;
  viz_msg.type = 4;
  viz_msg.header = header;
  viz_msg.header.stamp = ros::Time::now();
  viz_msg.scale.x = 0.2;
  viz_msg.color.a = 1.0;
  // DEAL WITH Actual MSG
  kr_planning_msgs::TrajectoryDiscretized traj;
  traj.header = header;  // TODO: ASK LAURA
  traj.header.stamp = ros::Time::now();
  traj.N_ctrl = N_controls;
  geometry_msgs::Point pos_t;
  geometry_msgs::Point vel_t;
  geometry_msgs::Point acc_t;
  geometry_msgs::Point moment_t;
  Eigen::Vector3d moment_t_v =
      Eigen::Vector3d::Zero();  // TODO:(Yifei) return proper moment
  unpack(moment_t, moment_t_v);
  opt_traj_.clear();

  if (solver_status == 0) {
    // no state here, need to calc from force, stay at 0 for now
    for (int i = 0; i < N_controls + 1; i++) {
      // deal with last index outside the loop
      // unpack points . THIS CODE IS SIMILAR TO NEXT SECTION FOR UNPACKING,
      // write a function!
      Eigen::Quaterniond q_unpack(X_sim[i][3],
                                  X_sim[i][4],
                                  X_sim[i][5],
                                  X_sim[i][6]);  // not eigen quat
      Eigen::Vector3d a_body;
      Eigen::Vector3d v_body;

      a_body << 0.0, 0.0,
          mpc_solver->model_ptr->kf_ * U_sim[i].sum() /
              mpc_solver->model_ptr->mass_;

      auto a_world = (q_unpack * a_body) - Eigen::Vector3d(0, 0, 9.81);
      // std::cout << "a_world = " << a_world << std::endl;
      // std::cout << "a_world = " << a_world << std::endl;

      v_body << X_sim[i][7], X_sim[i][8], X_sim[i][9];
      auto v_world = q_unpack * v_body;

      auto RPY = q_unpack.toRotationMatrix().eulerAngles(0, 1, 2);
      unpack(pos_t, X_sim[i]);  // x y z q1 q2 q3 q4 v1 v2 v3 w1 w2 w3
      unpack(vel_t, v_world);
      unpack(acc_t, a_world);
      Eigen::VectorXd pva_temp(9);
      pva_temp << X_sim[i].segment<3>(0), v_world, a_world;
      opt_traj_.push_back(pva_temp);
      // DEAL WITH VIZ
      viz_msg.points.push_back(pos_t);
      // DEAL WITH Actual MSG
      traj.pos.push_back(pos_t);
      traj.vel.push_back(vel_t);
      traj.acc.push_back(acc_t);
      traj.yaw.push_back(RPY[2]);
      if (i != N_controls) {
        traj.thrust.push_back(U_sim[i].sum());
        traj.moment.push_back(moment_t);  // 0
      }
      traj.dt = dt;
      // initial attitude and initial omega not used
    }
  }
  // DEAL WITH VIZ
  if (publish_viz_) opt_viz_pub_.publish(viz_msg);
  // DEAL WITH Actual MSG
  if (publish_optimized_traj_) opt_traj_pub_.publish(traj);

  if (solver_status == 0)
    return traj;
  else
    return kr_planning_msgs::TrajectoryDiscretized();
}

// if compute altro is on, then we will compute altro, otherwise we will just
// publish the discretized traj from passed in
kr_planning_msgs::TrajectoryDiscretized
SplineTrajSampler::sample_and_refine_trajectory(
    const Eigen::VectorXd& start_state,
    const kr_planning_msgs::SplineTrajectory& traj,
    const std::vector<Eigen::MatrixXd>& hPolys,
    const Eigen::VectorXd& allo_ts) {
  ROS_WARN("[iLQR]: sample_and_refine_trajectory Started");
  std::vector<std::pair<Eigen::MatrixXd, Eigen::VectorXd>> Abpolys(
      hPolys.size());
  std::transform(
      hPolys.begin(), hPolys.end(), Abpolys.begin(), pt_normal_to_halfspace);
  ROS_WARN("[iLQR]: msg Converted");

  double total_time = traj.data[0].t_total;
  int N_controls = total_time / dt_ / 5;
  N_controls *= 5;  // make sure a multiple of 5 for traj to be correctly fitted
                    // to poly outside
  int N_state = N_controls + 1;

  std::vector<Eigen::Vector3d> pos = sample(traj, total_time, N_state, 0);
  std::vector<Eigen::Vector3d> vel = sample(traj, total_time, N_state, 1);
  std::vector<Eigen::Vector3d> acc = sample(traj, total_time, N_state, 2);
  std::vector<Eigen::Vector3d> jerk = sample(traj, total_time, N_state, 3);
  std::vector<Eigen::Vector3d> snap = sample(traj, total_time, N_state, 4);
  ROS_INFO_STREAM("[iLQR]: Sampled " << pos.size() << " points"
                                     << "total time = " << total_time);
  double dt = dt_;
  ROS_WARN("[iLQR]: Sample Done");

  kr_planning_msgs::TrajectoryDiscretized traj_discretized;
  traj_discretized.header = traj.header;
  double jerk_abs_sum = 0.0;
  double total_ref_F = 0.0;
  double total_ref_M = 0.0;
  Eigen::Quaterniond q_return;
  Eigen::Vector3d w_return;
  geometry_msgs::Quaternion ini_q_msg;
  geometry_msgs::Point ini_w_msg;
  std::vector<double> yaw_ref(N_state);
  std::vector<double> thrust_ref(N_state);
  std::vector<Eigen::Vector3d> moment_ref(N_state);
  std::vector<Eigen::Quaterniond> q_ref(N_state);
  std::vector<Eigen::Vector3d> w_ref(N_state);
  geometry_msgs::Point pos_t;
  geometry_msgs::Point vel_t;
  geometry_msgs::Point acc_t;
  geometry_msgs::Point jerk_t;
  geometry_msgs::Point snap_t;
  double thrust;                // for return
  geometry_msgs::Point moment;  // for return

  for (int i = 0; i < N_state; i++) {
    Eigen::Vector3d yaw_three_der =
        Eigen::Vector3d(0, 0, 0);  // TODO: Assign from traj
    yaw_ref[i] = yaw_three_der[0];
    Eigen::Vector3d M = compute_ref_inputs(pos[i],
                                           vel[i],
                                           acc[i],
                                           jerk[i],
                                           snap[i],
                                           yaw_three_der,
                                           thrust,
                                           moment,
                                           q_return,
                                           w_return);
    moment_ref[i] = M;
    thrust_ref[i] = thrust;

    if (i == 0) {
      tf::quaternionEigenToMsg(q_return, ini_q_msg);
      tf::pointEigenToMsg(w_return, ini_w_msg);
    }
    // std::cout<<"rot q = "<<q_return.w()<<q_return.vec() <<std::endl;
    q_ref[i] = q_return;
    w_ref[i] = w_return;

    if (i != N_state - 1) {
      jerk_abs_sum += acc[i].squaredNorm() * dt;
      total_ref_F += abs(thrust) * dt;
      total_ref_M += M.norm() * dt;

      traj_discretized.thrust.push_back(thrust);
      traj_discretized.moment.push_back(moment);
    }

    unpack(pos_t, pos[i]);
    unpack(vel_t, vel[i]);
    unpack(acc_t, acc[i]);

    traj_discretized.pos.push_back(pos_t);
    traj_discretized.vel.push_back(vel_t);
    traj_discretized.acc.push_back(acc_t);
    traj_discretized.yaw.push_back(yaw_ref[i]);
  }

  if (write_summary_to_file_) {
    // get start location as key to which ones successed
    poly_array_file << std::to_string(N_iter_) + ", " +
                           std::to_string(pos[0](0)) + ", " +
                           std::to_string(jerk_abs_sum) + ", " +
                           std::to_string(total_time) + "," +
                           std::to_string(total_ref_F) + "," +
                           std::to_string(total_ref_M) + "\n";
    // compute the total snap at 100 points
  }

  // this is input t, has one less element since you need ref everywhere
  // std::vector<double> t = linspace(0.0, total_time, N_controls_ + 1);
  traj_discretized.dt = dt;
  traj_discretized.N_ctrl = N_controls;
  traj_discretized.inital_attitude = ini_q_msg;
  traj_discretized.initial_omega = ini_w_msg;

  ROS_WARN("[iLQR]: Sampling Finished");

  if (compute_altro_) {
    kr_planning_msgs::TrajectoryDiscretized refined_traj =
        publish_altro(start_state,
                      pos,
                      vel,
                      acc,
                      yaw_ref,
                      thrust_ref,
                      moment_ref,
                      dt,
                      q_ref,
                      w_ref,
                      Abpolys,
                      allo_ts,
                      traj.header);

    ROS_WARN("[iLQR]: Altro Finished");
    return refined_traj;
  } else {
    return traj_discretized;
  }
}

void SplineTrajSampler::callbackWrapper(
    const kr_planning_msgs::SplineTrajectory::ConstPtr& msg) {
  kr_planning_msgs::SplineTrajectory traj = *msg;

  ROS_ERROR("[iLQR]: Callback not working, need to import as class");
  // kr_planning_msgs::TrajectoryDiscretized traj_discretized =
  //     sample_and_refine_trajectory(traj);

  // ros::Duration(0.5).sleep(); // this is for yuwei's planner to finish so we
  // have polytope info ROS_ERROR("SLEEPING 0.5 sec to wait for polytope info");
  // this->sampled_traj_pub_.publish(traj_discretized);
  // do optimization here if needed
  N_iter_++;
}
