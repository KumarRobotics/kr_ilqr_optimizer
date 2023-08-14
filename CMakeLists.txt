cmake_minimum_required(VERSION 3.10)
project(kr_ilqr_optimizer)


## Compile as C++11, supported in ROS Kinetic and newer
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  kr_planning_msgs
  roscpp
  eigen_conversions
  visualization_msgs
  altro
)



## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)
find_package(Eigen3 3.4 REQUIRED) # for ilqr stuff

catkin_package(
 INCLUDE_DIRS include ${catkin_INCLUDE_DIRS}
 LIBRARIES kr_ilqr_optimizer
 CATKIN_DEPENDS kr_planning_msgs roscpp visualization_msgs altro
 DEPENDS # we should put fmt and eigen here
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
# include
  include
  ${catkin_INCLUDE_DIRS}
)

add_library(quaternion_utils src/quaternion_utils.cpp)
target_link_libraries(quaternion_utils
  PUBLIC
  Eigen3::Eigen
  )

add_library(finitediff 
  src/finitediff.cpp
) 
target_link_libraries(finitediff
  PUBLIC
  Eigen3::Eigen
  )

add_library(altro_test_utils src/test_utils.cpp)
target_link_libraries(altro_test_utils
  PUBLIC
  finitediff
  quaternion_utils
  ${catkin_LIBRARIES}
)
## Declare a C++ library
add_library(quad_mpc include/kr_ilqr_optimizer/quadMPC.hpp)
set_target_properties(quad_mpc PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(quad_mpc
  PUBLIC
  altro_test_utils
  )
target_compile_definitions(altro_test_utils
  PUBLIC
  ALTRO_TEST_DIR="/home/yifei/Documents/optimal_ctrl/altro/test"
)
## Add cmake target dependencies of the library
## as an example, code may need to be generated before libraries
## either from message generation or dynamic reconfigure
# add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Declare a C++ executable
## With catkin_make all packages are built within a single CMake context
## The recommended prefix ensures that target names across packages don't collide

## Rename C++ executable without prefix
## The above recommended prefix causes long target names, the following renames the
## target back to the shorter version for ease of user use
## e.g. "rosrun someones_pkg node" instead of "rosrun someones_pkg someones_pkg_node"
# set_target_properties(${PROJECT_NAME}_node PROPERTIES OUTPUT_NAME node PREFIX "")

## Add cmake target dependencies of the executable
## same as for the library above
# add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
add_library(kr_ilqr_optimizer src/spline_trajectory_sampler.cpp)
target_link_libraries(kr_ilqr_optimizer 
      PUBLIC quad_mpc)
## Specify libraries to link a library or executable target against
add_executable(sampler_pub src/run_sampler.cpp)
target_link_libraries(sampler_pub 
  PUBLIC  kr_ilqr_optimizer) # WHY DOES IT NEED FMT!!!
# 


# target_link_libraries(kr_ilqr_optimizer PRIVATE )

#############
## Install ##
#############

# all install targets should use catkin DESTINATION variables
# See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
# catkin_install_python(PROGRAMS
#   scripts/my_python_script
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark executables for installation
## See http://docs.ros.org/melodic/api/catkin/html/howto/format1/building_executables.html
# install(TARGETS ${PROJECT_NAME}_node
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark libraries for installation
## See http://docs.ros.org/melodic/api/catkin/html/howto/format1/building_libraries.html
# install(TARGETS ${PROJECT_NAME}
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
# )

## Mark cpp header files for installation
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
#   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
# )

## Mark other files for installation (e.g. launch and bag files, etc.)
# install(FILES
#   # myfile1
#   # myfile2
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )

#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_kr_ilqr_optimizer.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
