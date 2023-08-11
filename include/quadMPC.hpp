
#include "Eigen/Dense"
#include "altro/altro.hpp"
#include "fmt/core.h"
#include "test_utils.hpp"
#include <chrono>
#include <fstream>
#include <filesystem>
#include <string>
#include <iostream>
// #include "nlohmann/json.hpp"
// #include <unsupported/Eigen/AutoDiff>
// #include "gtest/gtest.h"
namespace fs = std::filesystem;
// using json = nlohmann::json;


using Eigen::MatrixXd;
using Eigen::MatrixXd;
// using Eigen::Vector;

using namespace altro;
using Vector = Eigen::Matrix<a_float, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<a_float, Eigen::Dynamic, Eigen::Dynamic>;
using Eigen::VectorXi;

// const std::string ALTRO_TEST_DIR ="/home/yifei/Documents/optimal_ctrl/try_altro";


class quadMPC {
 public:
  quadMPC(int N, double t_ref_, bool use_quaternion) : N(N), t_ref_(t_ref_), solver(N) {}
  void SetUp() {
    // Objective
    std::cout<< "Size of Problem N = "<< N << std::endl;
    Qd = Vector::Constant(n, 0);
    Rd = Vector::Constant(m, 0.1);
    Qdf = Vector::Constant(n, 0);
    const int three = 3;
    Qd.segment<3>(0) = Vector::Constant(three, 0.1);
    Qdf.segment<3>(0) = Vector::Constant(three, 0.1);
    Qdf.segment<3>(7) = Vector::Constant(three, 0.1);

    // Dynamics
    auto model_ptr = std::make_shared<quadModel>();
    ContinuousDynamicsFunction dyn0 = [model_ptr](double *x_dot, const double *x, const double *u) {
      model_ptr->Dynamics(x_dot, x, u);
    };
    ContinuousDynamicsJacobian jac0 = [model_ptr](double *jac, const double *x, const double *u) {
      model_ptr->Jacobian(jac, x, u);
    };

    ContinuousDynamicsJacobian jac_dt = [model_ptr](double *jac, const double *x, const double *u) {
      model_ptr->Jacobian_fd(jac, x, u);
    };
    dyn = MidpointDynamics(n, m, dyn0);
    jac = MidpointJacobian(n, m, dyn0, jac_dt);



    // Dimension and Time step
    err = solver.SetDimension(n, m);
    const int en = 12;
    const int em = 4;
    err = solver.SetErrorDimension(en, em);


    // Dynamics
    err = solver.SetExplicitDynamics(dyn, jac);

    // Read Reference Trajectory
    int N_ref = N;
    float t_ref = t_ref_;

    u_ref_single = Vector::Constant(m, model_ptr->get_hover_input());
    std::cout << "u_ref_single = " << u_ref_single << std::endl;
    // u_ref_single[0] = -1;
    // ReadScottyTrajectory(&N_ref, &t_ref, &x_ref, &u_ref); // this is where we should input the dispersion planner data

    // a_float * Qd_data = this->Qd.data();
    // std::cout << "Qd = " << Qd_data[0] << std::endl;
    // std::cout << "N_ref" << N_ref << std::endl;

    // Set time step equal to the reference trajectory
    h = t_ref / static_cast<double>(N_ref);
    printf("h = %f\n", h);
    err = solver.SetTimeStep(h);

    fmt::print("Set TimeStep Finished!\n");

    Eigen::Matrix<double,13,1> xf;
    xf << 0.,0.,0., 1.0, 0.0, 0.0, 0.0,  0,0,0, 0,0,0;
    for (int k = 0; k <= N-1; ++k) {
      // std::cout << k << std::endl;
      if (use_quaternion) {err = solver.SetQuaternionCost(n, m, Qd.data(), Rd.data(), 1.0, xf.data(), u_ref_single.data(), k);}
      else {err = solver.SetLQRCost(n, m, Qd.data(), Rd.data(), xf.data(), u_ref_single.data(), k);}
      // fmt::print("Print error\n");
    }
    // err = solver.SetLQRCost(n, m, Qdf.data(), Rd.data(), xf.data(), u_ref_single.data(), N);
    if (use_quaternion) {err = solver.SetQuaternionCost(n, m, Qdf.data(), Rd.data(), 1.0, xf.data(), u_ref_single.data(), N);}
    else {err = solver.SetLQRCost(n, m, Qdf.data(), Rd.data(), xf.data(), u_ref_single.data(), N);}
    // err = solver.SetQuaternionCost(n, m, Qdf.data(), Rd.data(), 1.0, xf.data(), u_ref_single.data(),N);  
    err = solver.SetInitialState(xf.data(), n);


    // Constraints
    const a_float max_thrust = model_ptr->max_thrust_per_prop;
    const a_float min_thrust = model_ptr->min_thrust_per_prop;
    auto actuator_con = [max_thrust, min_thrust](a_float *c, const a_float *x, const a_float *u) {
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
    auto actuator_jac = [](a_float *jac, const a_float *x, const a_float *u) {
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
    err = solver.SetConstraint(actuator_con, actuator_jac, 8, ConstraintType::INEQUALITY,
                               "Actuator Min Max", 0, N + 1);

    // Initialize Solver
    err = solver.Initialize();
    fmt::print("Solver Initialized!\n");

    // Solve
    AltroOptions opts;
    opts.verbose = Verbosity::Silent;
    opts.iterations_max = 40;
    opts.use_backtracking_linesearch = true;
    opts.quat_start_index = 3;  // THIS IS VERY IMPORTANT!
    solver.SetOptions(opts);

  }
  protected:
     int N; 
     double t_ref_;

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
    void solve_problem(std::vector<Eigen::Vector3d> pos, std::vector<Eigen::Vector3d> vel, std::vector<Eigen::Vector3d> acc, std::vector<double> yaw_ref,
                                std::vector<double> thrust, std::vector<Eigen::Vector3d> moment, std::vector<double> t, int N_input,  std::vector<Eigen::Quaterniond> q_ref, std::vector<Eigen::Vector3d> w_ref,
                                std::vector<Vector> & X_sim, std::vector<Vector> & U_sim, std::vector<double> & t_sim){
      (void) t;
      (void) N_input;

      // ErrorCodes err;
      // solver.Initialize();


      const int n_const = 13;

      Eigen::Matrix<double,n_const,1> x0;
      // x0 << 2.0, 1.0, -1.0,   -0.752,  0.443,  0.443, -0.206,  0.5,-0.5,1.0, 0.8,0.8,0.8;
      q_ref[0].normalize();
      x0 << pos[0], q_ref[0].w(), q_ref[0].vec(), vel[0], w_ref[0];
      // x0 << 2.0, 1.0, -1.0,  1.0,  0,   0, 0,  0.5,-0.5,0, 0.,0.,0.;

      Eigen::Matrix<double,n_const,1> xf;
      xf << pos.back(), 1.0, 0.0, 0.0, 0.0,  0,0,0, 0,0,0;
      // xf << 0.0, 0.0, 0.0,   1.0, 0.0, 0.0, 0.0,  0,0,0, 0,0,0;

      // std::cout <<"size of pos" << pos.size() << std::endl;
      // std::cout << "x0 = " << x0 << std::endl;
      // std::cout << "xf = " << xf << std::endl;
      // Set Tracking Cost Function
      for (int k = 0; k <= N-1; ++k) {
        Eigen::Matrix<double,13,1> x_ref_k;
        // std::cout<< "ref k = " << k << std::endl;
        // std::cout<< "pos[k] = " << pos[k] << std::endl;
        // std::cout<< "q_ref[k] = " << q_ref[k][1] << std::endl;

        x_ref_k << pos[k], q_ref[k].w(), q_ref[k].vec(), vel[k], w_ref[k];
        // std::cout << k << std::endl;
        // 
        
        if (use_quaternion) {err = solver.SetQuaternionCost(n, m, Qd.data(), Rd.data(), 1.0, x_ref_k.data(), u_ref_single.data(), k);}
        else{err = solver.SetLQRCost(n, m, Qd.data(), Rd.data(), x_ref_k.data(), u_ref_single.data(), k);}        
        // PrintErrorCode(err);
        // fmt::print("Print error\n");
        err = solver.SetState(x_ref_k.data(), n, k);

      }
      if (use_quaternion) {err = solver.SetQuaternionCost(n, m, Qdf.data(), Rd.data(), 1.0, xf.data(), u_ref_single.data(), N);}
      else {err = solver.SetLQRCost(n, m, Qdf.data(), Rd.data(), xf.data(), u_ref_single.data(), N);}
        // 
      err = solver.SetState(xf.data(), n, N);

        // Initial State
      err = solver.SetInitialState(x0.data(), n);

      // Set Initial Trajectory
      err = solver.SetInput(u_ref_single.data(), m);
      // fmt::print("Set Input Finished!\n");

      solver.OpenLoopRollout();
      double cost_initial = solver.CalcCost();
      // fmt::print("Initial cost = {}\n", cost_initial);

      
      SolveStatus status = solver.Solve();

      std::cout << "Solve status is: " << static_cast<unsigned int>(status) << std::endl; 
      cost_initial = solver.CalcCost();
      // fmt::print("Final cost = {}\n", cost_initial);
      
      // std::vector<Vector> X_sim;
      // std::vector<Vector> U_sim;
      // std::vector<float> t_sim;
      float t_now = 0;
      for (int k = 0; k <= N; k++) {
        Eigen::VectorXd x(n);
        solver.GetState(x.data(), k);
        X_sim.emplace_back(x);
        t_sim.emplace_back(t_now);
        t_now += solver.GetTimeStep(k);
        if (k != N) {
          Eigen::VectorXd u(m);
          solver.GetInput(u.data(), k);
          U_sim.emplace_back(u);
        }
      }
      // ROS_ERROR("SHUTTING DOWN SINGLE TEST");
      // ros::shutdown();
    }
   
    
};
