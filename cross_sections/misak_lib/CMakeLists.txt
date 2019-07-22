# Set up project
PROJECT(Misak_eN)
SET(LIBNAME "Misak_eN")
enable_language(Fortran)
include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

add_library(${LIBNAME} wrapper.cpp sigma_en_bound.f )

