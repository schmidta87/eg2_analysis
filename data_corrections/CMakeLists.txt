
# Include global headers
include_directories( ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

# Corrector
add_executable(corrector corrector.cpp)
target_link_libraries(corrector ${ROOT_LIBRARIES})






