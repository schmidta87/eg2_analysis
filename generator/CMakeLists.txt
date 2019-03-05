
# Include global headers
include_directories( ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

# Generator 
add_executable(gen_weight gen_weight.cpp Nuclear_Info.cpp )
target_link_libraries(gen_weight ${ROOT_LIBRARIES} Cross_Sections)

# Event Simulator
add_executable(simulator simulator.cpp)
target_link_libraries(simulator ${ROOT_LIBRARIES} AccMap Fid)

# ep_cut
add_executable(ep_cut ep_cut.cpp)
target_link_libraries(ep_cut ${ROOT_LIBRARIES})

# Make histograms
add_executable(make_hists make_hists.cpp )
target_link_libraries(make_hists ${ROOT_LIBRARIES} Fid Cross_Sections AccMap)

# Make ratios (when combining more than one histogram file)
add_executable(make_ratios make_ratios.cpp)
target_link_libraries(make_ratios ${ROOT_LIBRARIES})

add_subdirectory(Jackson)

add_subdirectory(Andrew)
