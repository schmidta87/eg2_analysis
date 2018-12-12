
# Include global headers
include_directories( ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

# Generator 
add_executable(gen_weight gen_weight.cpp Nuclear_Info.cpp )
target_link_libraries(gen_weight ${ROOT_LIBRARIES} Cross_Sections)

# pp/np ratio cut
add_executable(pp_np_cut pp_np_cut.cpp)
target_link_libraries(pp_np_cut ${ROOT_LIBRARIES} AccMap Fid)

# Event Simulator
add_executable(simulator simulator.cpp)
target_link_libraries(simulator ${ROOT_LIBRARIES} AccMap Fid)

# Make histograms
add_executable(make_hists make_hists.cpp )
target_link_libraries(make_hists ${ROOT_LIBRARIES} Fid)

# pp_np_prec
add_executable(pp_np_prec pp_np_prec.cpp )
target_link_libraries(pp_np_prec ${ROOT_LIBRARIES} Fid)



