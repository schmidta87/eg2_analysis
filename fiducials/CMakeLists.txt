# Set project name
PROJECT(Fiducials)
SET(LIBNAME "Fid")

include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

# Fiducial Function Library
add_library(${LIBNAME} fiducials.cpp)
target_link_libraries(${LIBNAME} ${ROOT_LIBRARIES} Helpers)


# Test fiducials executable
add_executable(fid_test fid_test.cpp)
target_link_libraries(fid_test ${ROOT_LIBRARIES} Fid)


# Show fiducial cut results on epp events
add_executable(fid_recoil_test fid_recoil_test.cpp)
target_link_libraries(fid_recoil_test ${ROOT_LIBRARIES} Fid)


# Show fiducial cut results on electrons
add_executable(fid_e_test fid_e_test.cpp)
target_link_libraries(fid_e_test ${ROOT_LIBRARIES} Fid)


# Show fiducial cut results on protons
add_executable(fid_p_test fid_p_test.cpp)
target_link_libraries(fid_p_test ${ROOT_LIBRARIES} Fid)
