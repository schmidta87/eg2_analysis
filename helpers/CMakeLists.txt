# Set project name
PROJECT(Helpers)
SET(LIBNAME "Helpers")

include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

# Fiducial Function Library
add_library(${LIBNAME} helpers.cpp)
target_link_libraries(${LIBNAME} ${ROOT_LIBRARIES} )
