# Set project name
PROJECT(AccMap)
SET(LIBNAME "Fid")

include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

# Acceptance Map Library
add_library(${LIBNAME} fiducials.cpp)
target_link_libraries(${LIBNAME} ${ROOT_LIBRARIES} )


# Query map executable
add_executable(fid_test fid_test.cpp)
target_link_libraries(fid_test ${ROOT_LIBRARIES} Fid)
