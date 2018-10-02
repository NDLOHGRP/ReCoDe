cmake_minimum_required (VERSION 2.6)
project (ReCoDe)
# The version number.
set (ReCoDe_VERSION_MAJOR 1)
set (ReCoDe_VERSION_MINOR 0)

cmake_policy(SET CMP0015 NEW)
cmake_policy(SET CMP0014 NEW)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "src/ReCoDeConfig.h.in"
  "bin/ReCoDeConfig.h"
  )
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories(include)
link_directories(lib)
 
add_executable(ReCoDe src/recode.cpp)

target_link_libraries (ReCoDe zlib)