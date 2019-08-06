# Directory for Axel's programs

# Include global headers
include_directories( ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

# Program for testing how many recoils are mistaken as leads
add_executable(mistaken_recoils mistaken_recoils.cpp )
target_link_libraries(mistaken_recoils ${ROOT_LIBRARIES} Fid AccMap)

# Program for testing LC weights on a2
add_executable(lc_a2_test lc_a2_test.cpp )
target_link_libraries(lc_a2_test ${ROOT_LIBRARIES})
