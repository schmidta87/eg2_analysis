# Set project name
PROJECT(Cross_Sections)
SET(LIBNAME "Cross_Sections")

include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

# Misak Lib
add_subdirectory(misak_lib)

# Cross Sections Library
add_library(${LIBNAME} Cross_Sections.cpp)
target_link_libraries(${LIBNAME} ${ROOT_LIBRARIES} Helpers Misak_eN)

# Cross Section test executable
add_executable(cs_query cs_query.cpp)
target_link_libraries(cs_query ${ROOT_LIBRARIES} Cross_Sections )

