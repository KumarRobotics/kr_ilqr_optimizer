# KR iLQR Optimizer
## Maintainer: 
Yifei Simon Shao (yishao at seas.upenn.edu)
## Description
This package connects the Altro solver with the KR_autonomous flight stack to do constrained trajectory optimization.
It can work in 2 ways
1. Run the solver in a separate node (not exactly working)
2. Functions in this package will are called in the action_planner package to run the solver easily
## Dependencies
1. [Altro](https://github.com/shaoyifei96/altro.git) This should automatically pull all dependencies from the internet
## ToDos
1. Build python binding for generating motion primitives
2. MPC style initialization, make it much faster
3. Do proper control input initialization
