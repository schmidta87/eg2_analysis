
# Include global headers
include_directories( ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

# Electron energy corrections
add_executable(electron_E_corr electron_E_corr.cpp)
target_link_libraries(electron_E_corr ${ROOT_LIBRARIES})






