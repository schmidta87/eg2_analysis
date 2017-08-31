# CMakeLists for the eg2 analysis repository

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

# Project name
project (eg2)

# ADD SUBDIRECTORIES

include_directories(${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(acc_corr acc_corr.cpp)
target_link_libraries(acc_corr ${ROOT_LIBRARIES})

