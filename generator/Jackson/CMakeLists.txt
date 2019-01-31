
# Include global headers
include_directories( ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

# pp/np ratio cut
add_executable(pp_np_cut pp_np_cut.cpp)
target_link_libraries(pp_np_cut ${ROOT_LIBRARIES} AccMap Fid Cross_Sections)

# pn/p ratio cut
add_executable(pn_p_cut pn_p_cut.cpp)
target_link_libraries(pn_p_cut ${ROOT_LIBRARIES} AccMap Fid Cross_Sections)

# pp_np_prec
add_executable(prec_hists prec_hists.cpp )
target_link_libraries(prec_hists ${ROOT_LIBRARIES} Fid)

# pn_p
add_executable(pn_p_hists pn_p_hists.cpp )
target_link_libraries(pn_p_hists ${ROOT_LIBRARIES} Fid)

# finder
add_executable(gapHists gapHists.cpp )
target_link_libraries(gapHists ${ROOT_LIBRARIES} Fid)

