# CMakeLists for the eg2 analysis repository

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")

#OpenMP
ENABLE_LANGUAGE(CXX)
OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
 FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  ENDIF()
ENDIF()

# ROOT
find_package(ROOT REQUIRED)

# Project name
project (eg2)

# ADD SUBDIRECTORIES
add_subdirectory(helpers)
add_subdirectory(acceptance)
add_subdirectory(fiducials)
add_subdirectory(generator)

include_directories(${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )

# Add executable for acceptance corrections
add_executable(mcmc_kde mcmc_kde.cpp)
target_link_libraries(mcmc_kde ${ROOT_LIBRARIES} Helpers)
add_executable(mcmc_bin mcmc_bin.cpp)
target_link_libraries(mcmc_bin ${ROOT_LIBRARIES} AccMap Helpers Fid)

# Add executable for generating pseudo data
add_executable(gen_pseudo gen_pseudo.cpp)
target_link_libraries(gen_pseudo ${ROOT_LIBRARIES} AccMap)

# Add executable for printing the tree from mcmc
add_executable(print_tree print_tree.cpp)
target_link_libraries(print_tree ${ROOT_LIBRARIES})

# Add executable for processing burn-in and autocorrelation
add_executable(process_tree process_tree.cpp)
target_link_libraries(process_tree ${ROOT_LIBRARIES})

# Add executable for doing basic fits
add_executable(fit_basic fit_basic.cpp)
target_link_libraries(fit_basic ${ROOT_LIBRARIES} Fid)

# Add executable for calculating the pp/p correction
add_executable(find_correction find_correction.cpp)
target_link_libraries(find_correction ${ROOT_LIBRARIES} AccMap)

# Add executable for calculating the pp/p correction
add_executable(find_corr_single find_corr_single.cpp)
target_link_libraries(find_corr_single ${ROOT_LIBRARIES} AccMap)

# Add executable for calculating the pp/p ratio from data
add_executable(pp_to_p pp_to_p.cpp)
target_link_libraries(pp_to_p ${ROOT_LIBRARIES} Fid)
