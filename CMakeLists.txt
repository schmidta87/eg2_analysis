# CMakeLists for the eg2 analysis repository

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

# Project name
project (eg2)

# ADD SUBDIRECTORIES

include_directories(${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )

# Add executable for acceptance corrections
add_executable(acc_corr acc_corr.cpp)
target_link_libraries(acc_corr ${ROOT_LIBRARIES})

# Add executable for generating pseudo data
add_executable(gen_pseudo gen_pseudo.cpp)
target_link_libraries(gen_pseudo ${ROOT_LIBRARIES})

# Add executable for printing the tree from mcmc
add_executable(print_tree print_tree.cpp)
target_link_libraries(print_tree ${ROOT_LIBRARIES})

# Add executable for processing burn-in and autocorrelation
add_executable(process_tree process_tree.cpp)
target_link_libraries(process_tree ${ROOT_LIBRARIES})

# Add executable for doing basic fits
add_executable(fit_basic fit_basic.cpp)
target_link_libraries(fit_basic ${ROOT_LIBRARIES})

