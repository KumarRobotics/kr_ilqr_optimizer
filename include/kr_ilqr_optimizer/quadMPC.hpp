#include <ros/ros.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "Eigen/Dense"
#include "altro/altro.hpp"
#include "kr_ilqr_optimizer/test_utils.hpp"
namespace fs = std::filesystem;
// using json = nlohmann::json;

using Eigen::MatrixXd;
// using Eigen::Vector;

using namespace altro;
using Vector = Eigen::Matrix<a_float, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<a_float, Eigen::Dynamic, Eigen::Dynamic>;
using Eigen::VectorXi;

// const std::string ALTRO_TEST_DIR
// ="/home/yifei/Documents/optimal_ctrl/try_altro";

class quadMPC {
 public:
  quadMPC(int N_ctrl, bool use_quaternion) : N_ctrl_(N_ctrl), solver(N_ctrl) {}
  void SetUp(double mass,
             double grav_const,
             Eigen::Matrix3d inertia,
             double min_thrust,
             double max_thrust,
             double kf,
             double km,
             double L) {
    model_ptr = std::make_shared<const quadModel>(
        mass,
        grav_const,
        inertia,
        min_thrust / 4 / kf,  // converted to w^2
        max_thrust / 4 / kf,
        kf,
        km,
        L);

    // Objective
    std::cout << "Size of Problem N_ctrl = " << N_ctrl_ << std::endl;
  }
  std::shared_ptr<const quadModel> model_ptr;

 protected:
  int N_ctrl_;

  const int n = quadModel::NumStates;
  const int m = quadModel::NumInputs;
  float h;

  Vector Qd;
  Vector Rd;
  Vector Qdf;
  bool use_quaternion;

  ExplicitDynamicsFunction dyn;
  ExplicitDynamicsJacobian jac;

  ALTROSolver solver;

  ErrorCodes err;
  // Reference Trajectory (the "Scotty Dog")
  std::vector<Eigen::Matrix<double, 13, 1>> x_ref;
  std::vector<Eigen::Vector4d> u_ref;
  Vector u0;
  Eigen::Vector4d u_ref_single;

 public:
  uint solve_problem(
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
      std::vector<Vector>& X_sim,
      std::vector<Vector>& U_sim,
      std::vector<double>& t_sim) {
    solver = ALTROSolver(N_ctrl_);
    Qd = Vector::Constant(n, 0);
    Rd = Vector::Constant(m, 0.1);
    Qdf = Vector::Constant(n, 0);
    const int three = 3;
    Qd.segment<3>(0) = Vector::Constant(three, 0.1);
    Qdf.segment<3>(0) = Vector::Constant(three, 0.1);
    Qdf.segment<3>(10) = Vector::Constant(three, 0.2);
    Qdf.segment<3>(7) = Vector::Constant(three, 0.2);

    // Dynamics
    ContinuousDynamicsFunction dyn0 =
        [this](double* x_dot, const double* x, const double* u) {
          model_ptr->Dynamics(x_dot, x, u);
        };
    ContinuousDynamicsJacobian jac0 =
        [this](double* jac, const double* x, const double* u) {
          model_ptr->Jacobian(jac, x, u);
        };

    ContinuousDynamicsJacobian jac_dt =
        [this](double* jac, const double* x, const double* u) {
          model_ptr->Jacobian_fd(jac, x, u);
        };
    dyn = MidpointDynamics(n, m, dyn0);
    jac = MidpointJacobian(n, m, dyn0, jac0);

    // Dimension and Time step
    err = solver.SetDimension(n, m);
    const int en = 12;
    const int em = 4;
    err = solver.SetErrorDimension(en, em);

    // Dynamics
    err = solver.SetExplicitDynamics(dyn, jac);

    // Read Reference Trajectory
    int N_state = N_ctrl_ + 1;
    // float t_total = t_total_;

    u_ref_single = Vector::Constant(m, model_ptr->get_hover_input());
    std::cout << "u_ref_single = " << u_ref_single << std::endl;

    // Set time step equal to the reference trajectory
    printf("h = %f\n", dt);
    err = solver.SetTimeStep(dt);

    Eigen::Matrix<double, 13, 1> xf;
    xf << pos.back(), 1.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0;
    // Suprise ! Specifying COST FUNCTION TWICE is useful!!
    for (int k = 0; k < N_ctrl_; ++k) {
      if (use_quaternion) {
        err = solver.SetQuaternionCost(
            n, m, Qd.data(), Rd.data(), 1.0, xf.data(), u_ref_single.data(), k);
      } else {
        err = solver.SetLQRCost(
            n, m, Qd.data(), Rd.data(), xf.data(), u_ref_single.data(), k);
      }
    }
    // Set FINAL Cost
    if (use_quaternion) {
      err = solver.SetQuaternionCost(
          n,
          m,
          Qdf.data(),
          Rd.data(),
          1.0,
          xf.data(),
          u_ref_single.data(),
          N_state - 1);  // Set it for N_state, but index -1
    } else {
      err = solver.SetLQRCost(n,
                              m,
                              Qdf.data(),
                              Rd.data(),
                              xf.data(),
                              u_ref_single.data(),
                              N_state - 1);
    }
    // err = solver.SetQuaternionCost(n, m, Qdf.data(), Rd.data(), 1.0,
    // xf.data(), u_ref_single.data(),N);
    err = solver.SetInitialState(xf.data(), n);

    // Initialize Solver
    // err = solver.Initialize();
    // std::cout << "Solver Initialized!\n" << std::endl;

    // Solve
    AltroOptions opts;
    opts.verbose = Verbosity::Silent;
    opts.iterations_max = 40;
    opts.use_backtracking_linesearch = true;
    opts.quat_start_index = 3;  // THIS IS VERY IMPORTANT!
    opts.use_quaternion = use_quaternion;
    solver.SetOptions(opts);
    // Constraints
    const a_float max_thrust = model_ptr->max_w_sq;
    const a_float min_thrust = model_ptr->min_w_sq;
    std::cout<< "min max propeller limits = " << min_thrust<< "  "<<max_thrust<< std::endl;
    auto actuator_con = [max_thrust, min_thrust](
                            a_float* c, const a_float* x, const a_float* u) {
      (void)x;
      // < 0 format
      c[0] = u[0] - max_thrust;
      c[1] = u[1] - max_thrust;
      c[2] = u[2] - max_thrust;
      c[3] = u[3] - max_thrust;
      c[4] = min_thrust - u[0];
      c[5] = min_thrust - u[1];
      c[6] = min_thrust - u[2];
      c[7] = min_thrust - u[3];
    };
    auto actuator_jac = [](a_float* jac, const a_float* x, const a_float* u) {
      (void)x;
      (void)u;
      Eigen::Map<Eigen::Matrix<a_float, 8, 17>> J(jac);
      J.setZero();
      J(0, 13) = 1.0;
      J(1, 14) = 1.0;
      J(2, 15) = 1.0;
      J(3, 16) = 1.0;
      J(4, 13) = -1.0;
      J(5, 14) = -1.0;
      J(6, 15) = -1.0;
      J(7, 16) = -1.0;
    };
    auto state_con = [xf](a_float* c, const a_float* x, const a_float* u) {
      (void)u;
      // Eigen::Map<const Vector> X(x, 13);
      // Eigen::Map<Vector> C(c, 6);
      // constrain velocity and angular velocity to 0
      c[0] = x[7];  // v = 0
      c[1] = x[8];
      c[2] = x[9];
      c[3] = x[10];  // w= 0
      c[4] = x[11];
      c[5] = x[12];
      // c[6] = x[0] - xf(0);  // xf = xf
      // c[7] = x[1] - xf(1);
      // c[8] = x[2] - xf(2);
      // c[9] = x[3] - 1.0;  // qf = 0 #We do not want to constrain the quad to
      // face a direction c[10] = x[4]; c[11] = x[5]; c[12] = x[6];
    };
    auto state_jac = [](a_float* jac, const a_float* x, const a_float* u) {
      (void)u;
      Eigen::Map<Eigen::Matrix<a_float, 6, 17>> J(jac);
      J.setZero();
      J(0, 7) = 1.0;
      J(1, 8) = 1.0;
      J(2, 9) = 1.0;
      J(3, 10) = 1.0;
      J(4, 11) = 1.0;
      J(5, 12) = 1.0;
      // J(6, 0) = 1.0;  // to get to goal properly
      // J(7, 1) = 1.0;
      // J(8, 2) = 1.0;
    };
    err = solver.SetConstraint(state_con,
                               state_jac,
                               6,
                               ConstraintType::EQUALITY,
                               "State",
                               N_state - 1);

    auto state_con_z = [](a_float* c, const a_float* x, const a_float* u) {
      (void)u;
      // Eigen::Map<const Vector> X(x, 13);
      // Eigen::Map<Vector> C(c, 6);
      // constrain velocity and angular velocity to 0
      c[0] = -x[2];
    };
    auto state_jac_z = [](a_float* jac, const a_float* x, const a_float* u) {
      (void)u;
      Eigen::Map<Eigen::Matrix<a_float, 1, 17>> J(jac);
      J.setZero();
      J(0, 2) = -1.0;
    };

    err = solver.SetConstraint(state_con_z,
                               state_jac_z,
                               1,
                               ConstraintType::INEQUALITY,
                               "StateZ",
                               AllIndices);
    err = solver.SetConstraint(actuator_con,
                               actuator_jac,
                               8,
                               ConstraintType::INEQUALITY,
                               "Actuator Min Max",
                               AllIndices);

    // int cur_poly_idx = 0;
    std::cout << "allots= " << allo_ts << std::endl;
    std::cout << "dt" << dt << std::endl;
    uint index_so_far = 0;
    double time_so_far = 0.0;
    for (int k = 0; k < h_polys.size(); k++) {
      auto h_poly = h_polys[k];
      const int num_con = h_poly.second.rows();
      auto poly_con = [h_poly](a_float* c, const a_float* x, const a_float* u) {
        (void)u;
        const int n_con = h_poly.second.rows();
        Eigen::Map<const Vector> X(x, 13);
        Eigen::Map<Vector> C(c, n_con);
        C = h_poly.first * X.segment<3>(0) - h_poly.second;
      };
      auto poly_jac = [h_poly](
                          a_float* jac, const a_float* x, const a_float* u) {
        (void)u;
        const int n_con = h_poly.second.rows();
        Eigen::Map<Eigen::MatrixXd> J(jac, n_con, 17);
        J.setZero();
        J.block(0, 0, n_con, 3) = h_poly.first;
      };
      double time_end = time_so_far + allo_ts[k];
      uint idx_end = int(time_end / dt);
      if (idx_end > N_state) {
        ROS_ERROR_STREAM("idx_end = " << idx_end << " > N_state");
        idx_end = N_state;
      }
      if (idx_end == index_so_far) {
        ROS_ERROR("polytope contains no points");
      } else {
        err = solver.SetConstraint(poly_con,
                                   poly_jac,
                                   num_con,
                                   ConstraintType::INEQUALITY,
                                   "polytope",
                                   index_so_far,
                                   idx_end);
      }
      std::cout << "k= " << k << " startidx = " << index_so_far
                << " endidx = " << idx_end << " time_end = " << time_end
                << std::endl;
      index_so_far = idx_end;
      time_so_far = time_end;
    }
    // auto poly_con = [h_polys](a_float* c, const a_float* x, const a_float* u)
    // {
    //   (void)x;
    //   // < 0 format
    //   c =
    // };
    solver.Initialize();
    // ErrorCodes err;
    // solver.Initialize();
    err = solver.SetTimeStep(dt);
    PrintErrorCode(err);

    const int n_const = 13;

    Eigen::Matrix<double, n_const, 1> x0;
    // x0 << 2.0, 1.0, -1.0,   -0.752,  0.443,  0.443, -0.206,  0.5,-0.5,1.0,
    // 0.8,0.8,0.8;
    Eigen::AngleAxisd yawAngle(start_state[6], Eigen::Vector3d::UnitZ());
    Eigen::Quaterniond q0(yawAngle);
    x0 << start_state.segment<3>(0), q0.w(), q0.vec(),
        start_state.segment<3>(3), 0, 0, 0;
    // x0 << 2.0, 1.0, -1.0,  1.0,  0,   0, 0,  0.5,-0.5,0, 0.,0.,0.;

    // Eigen::Matrix<double, n_const, 1> xf;

    std::cout << "x0 = " << x0 << std::endl;
    std::cout << "xf = " << xf << std::endl;
    // xf << 0.0, 0.0, 0.0,   1.0, 0.0, 0.0, 0.0,  0,0,0, 0,0,0;

    // std::cout <<"size of pos" << pos.size() << std::endl;
    // std::cout << "x0 = " << x0 << std::endl;
    // std::cout << "xf = " << xf << std::endl;
    // Set Tracking Cost Function
    for (int k = 0; k < N_ctrl_; ++k) {
      Eigen::Matrix<double, 13, 1> x_ref_k;
      // std::cout<< "ref k = " << k << std::endl;
      // std::cout<< "pos[k] = " << pos[k] << std::endl;
      // std::cout<< "q_ref[k] = " << q_ref[k][1] << std::endl;

      x_ref_k << pos[k], q_ref[k].w(), q_ref[k].vec(), vel[k], w_ref[k];
      // std::cout << k << std::endl;
      //

      if (use_quaternion) {
        err = solver.SetQuaternionCost(n,
                                       m,
                                       Qd.data(),
                                       Rd.data(),
                                       1.0,
                                       x_ref_k.data(),
                                       u_ref_single.data(),
                                       k);
      } else {
        err = solver.SetLQRCost(
            n, m, Qd.data(), Rd.data(), x_ref_k.data(), u_ref_single.data(), k);
      }
      // PrintErrorCode(err);
      // fmt::print("Print error\n");
      err = solver.SetState(x_ref_k.data(), n, k);
      // Set Initial Trajectory
      Eigen::Vector4d fm(4);  // Initialize a 4D Eigen vector
      // std::cout << "k = " << k << std::endl;
      // std::cout << "thrust[k] = " << thrust[k] << std::endl;
      // std::cout << "moment[k] = " << moment[k] << std::endl;

      fm << thrust[k], moment[k];
      Eigen::Vector4d reference_input = model_ptr->forceMatrix_inv() * fm;
      // std::cout << "reference_input = " << reference_input << std::endl;
      // we can change to reference input if needed but it is not getting
      // anywhere close to the goal, when there are obstacles it's even worse
      // since some constraints are violated
      err = solver.SetInput(u_ref_single.data(), m, k);
    }
    // Set Final COST!!
    if (use_quaternion) {
      err = solver.SetQuaternionCost(n,
                                     m,
                                     Qdf.data(),
                                     Rd.data(),
                                     1.0,
                                     xf.data(),
                                     u_ref_single.data(),
                                     N_state - 1);
    } else {
      err = solver.SetLQRCost(n,
                              m,
                              Qdf.data(),
                              Rd.data(),
                              xf.data(),
                              u_ref_single.data(),
                              N_state - 1);
    }
    //
    err = solver.SetState(xf.data(), n, N_state - 1);

    // Initial State
    err = solver.SetInitialState(x0.data(), n);

    // fmt::print("Set Input Finished!\n");

    solver.OpenLoopRollout();
    double cost_initial = solver.CalcCost();
    // fmt::print("Initial cost = {}\n", cost_initial);
    SolveStatus status = SolveStatus::Success;
    try {
      status = solver.Solve();
    } catch (std::exception& e) {
      std::cout << "Exception caught: " << e.what() << std::endl;
      status = SolveStatus::Unsolved;
    }

    std::cout << "Solve status is: " << static_cast<uint>(status) << std::endl;
    cost_initial = solver.CalcCost();
    // fmt::print("Final cost = {}\n", cost_initial);

    // std::vector<Vector> X_sim;
    // std::vector<Vector> U_sim;
    // std::vector<float> t_sim;
    if (static_cast<uint>(status) != 0){
      ROS_ERROR_STREAM("Solve failed with status: " << static_cast<uint>(status));
      return static_cast<uint>(status);
    }

    float t_now = 0;
    for (int k = 0; k < N_state; k++) {  // should be of size N + 1
      Eigen::VectorXd x(n);
      solver.GetState(x.data(), k);
      X_sim.emplace_back(x);
      t_sim.emplace_back(t_now);
      t_now += solver.GetTimeStep(k);
      if (k != N_state - 1) {
        Eigen::VectorXd u(m);
        solver.GetInput(u.data(), k);  // this is in w_sq now!!! not force
        U_sim.emplace_back(u);
      } else {
        std::cout << "last state = " << x << std::endl;
      }
    }
    // ros::shutdown();
    return static_cast<uint>(status);
  }
};
