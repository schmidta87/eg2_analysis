# Include global headers
include_directories( ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

# Make histograms
add_executable(TGAE2T TGAE2T.cpp )
target_link_libraries(TGAE2T ${ROOT_LIBRARIES})

# Make histograms
add_executable(H2T H2T.cpp )
target_link_libraries(H2T ${ROOT_LIBRARIES})

# Make file for Emiss
add_executable(Emiss Emiss.cpp )
target_link_libraries(Emiss ${ROOT_LIBRARIES})
