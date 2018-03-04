# Set project name
PROJECT(AccMap)
SET(LIBNAME "AccMap")

include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

# Acceptance Map Library
add_library(${LIBNAME} AccMap.cpp)
target_link_libraries(${LIBNAME} ${ROOT_LIBRARIES} )


# Query map executable
add_executable(map_test map_test.cpp)
target_link_libraries(map_test ${ROOT_LIBRARIES} AccMap)
