
# Include global headers
include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

# Generator 
add_executable(gen_weight gen_weight.cpp Nuclear_Info.cpp Nuclear_Info.h)
target_link_libraries(gen_weight ${ROOT_LIBRARIES})

# Event Simulator
add_executable(simulator simulator.cpp)
target_link_libraries(simulator ${ROOT_LIBRARIES} AccMap)
