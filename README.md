# KR iLQR Optimizer
## Maintainer: 
Yifei Shao (yishao at seas.upenn.edu)
## Description
This package hooks up the Altro solver with the KR_autonomous flight stack to do constrained trajectory optimization.
It can work in 2 ways
1. Run the solver in a separate node
2. Functions in this package will are called in the action_planner package to run the solver easily
## Dependencies
1. [Altro](https://github.com/bjack205/altro): When building altro, remember to not include the -march=native flag in the CMakeLists.txt file. After building, install it as a staic library so find_packages() can find it. **Forked Altro is not yet publically avaible, and I don't know if the original version works.**
2. Eigen > 3.4.0 
## ToDos
1. Don't use finite difference to compute the jacobian of dynamics.
2. Move Altro Fork as a dependency inside this repo
3. Build planning_detail option to run this.
4. Build python binding for generating motion primitives
5. MPC style initialization, make it much faster
6. Add polyotpe constraints
