# Directory for Axel's programs

# Include global headers
include_directories( ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

# Program for testing how many recoils are mistaken as leads
add_executable(mistaken_recoils mistaken_recoils.cpp )
target_link_libraries(mistaken_recoils ${ROOT_LIBRARIES} Fid AccMap)

