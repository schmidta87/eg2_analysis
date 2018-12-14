
# Include global headers
include_directories( ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

# pp/np ratio cutvy
add_executable(pp_np_cut pp_np_cut.cpp)
target_link_libraries(pp_np_cut ${ROOT_LIBRARIES} AccMap Fid Cross_Sections)

# pp_np_prec
add_executable(prec_hists prec_hists.cpp )
target_link_libraries(prec_hists ${ROOT_LIBRARIES} Fid)



