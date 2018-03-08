cmake_minimum_required(VERSION 2.8.12)
project(ocam_cpp)

option(VERBOSE "More verbose output while configuring/building" OFF)
##=======================================================================
## CMAKE_MODULE_PATH
##=======================================================================
## Make the scripts available in the 'cmake' directory available for the
## 'include()' and 'find_package()' commands.
#include(${PROJECT_SOURCE_DIR}/../../cmake/CMakeModulePath.cmake)

##=======================================================================
## Dependencies
##=======================================================================

##-----------------------------------------------------------------------
## catkin
##-----------------------------------------------------------------------
find_package(catkin REQUIRED COMPONENTS roscpp)

##=======================================================================
## Catkin Setup
##=======================================================================
#catkin_package(CATKIN_DEPENDS  camera_calibration
#               INCLUDE_DIRS include
#               LIBRARIES ${PROJECT_NAME}_lib
#              )

##------------------------------------------------------------------------
## libnavvis
##------------------------------------------------------------------------
find_package(NAVVIS REQUIRED)
list(APPEND INCLUDE_DIRS ${NAVVIS_INCLUDE_DIRS})
list(APPEND LIBRARY_DIRS ${NAVVIS_LIBRARY_DIRS})
list(APPEND LIBRARIES    ${NAVVIS_LIBRARIES})

##------------------------------------------------------------------------------
## Boost
##------------------------------------------------------------------------------
find_package(Boost REQUIRED COMPONENTS filesystem program_options thread QUIET)
list(APPEND INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
list(APPEND LIBRARY_DIRS ${Boost_LIBRARY_DIRS})
list(APPEND LIBRARIES    ${Boost_LIBRARIES})

##------------------------------------------------------------------------------
## Eigen
##------------------------------------------------------------------------------
find_package(Eigen3 REQUIRED QUIET)
list(APPEND INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
message(STATUS "EIGEN3_VERSION: ${EIGEN3_VERSION}")

##------------------------------------------------------------------------
## ceres
##------------------------------------------------------------------------
find_package(Ceres REQUIRED)
list(APPEND INCLUDE_DIRS ${CERES_INCLUDES})
list(APPEND LIBRARIES    ${CERES_LIBRARIES})

##========================================================================
## Set compiler flags
##========================================================================
include(DefaultCompilerFlags)
include(UseCxx11)
include(UseOpenMP)

##========================================================================
## Build
##========================================================================
SET(CMAKE_INSTALL_PREFIX "/opt/NavVis/")

##=============================================================================
## Build
##=============================================================================

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(src include ${INCLUDE_DIRS} ${catkin_INCLUDE_DIRS})

link_directories(${LIBRARY_DIRS} ${catkin_LIBRARY_DIRS})

add_executable(ocam src/main.cpp )
target_link_libraries(ocam ${LIBRARIES})
add_dependencies(ocam ${catkin_EXPORTED_TARGETS})
