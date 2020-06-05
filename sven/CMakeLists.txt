# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (epr CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)
find_package (BZip2)

set (CMAKE_PREFIX_PATH "/home/infri/develop/seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "/home/infri/develop/seqan/include")

enable_testing()
# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)
find_package (OpenMP REQUIRED)
find_package(benchmark REQUIRED)
# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

#add_subdirectory(benchmark)
# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})



#include (ExternalProject)
#set (GBENCHMARK_CMAKE_ARGS "-DBENCHMARK_ENABLE_TESTING=false;-DCMAKE_BUILD_TYPE=Release;-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}")
#ExternalProject_Add (
#    gbenchmark_project
#    PREFIX gbenchmark_project
#    GIT_REPOSITORY "https://github.com/google/benchmark.git"
#    GIT_TAG "v1.4.1"
#    CMAKE_ARGS ${GBENCHMARK_CMAKE_ARGS}
#    UPDATE_DISCONNECTED yes
#)

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
#set (CMAKE_CXX_FLAGS "-std=c++14 -O0 ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -I/home/sven/devel/sdsl-lite/include")
#-fsanitize=undefined
set (CMAKE_CXX_FLAGS "-fuse-ld=gold -fopenmp -march=native -std=c++17 -pthread -DSEQAN_DISABLE_VERSION_CHECK=YES ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.

add_executable (benchmark_epr benchmark_epr.cpp)
#add_dependencies(benchmark_epr gbenchmark_project)
#set_source_files_properties(benchmark_epr.cpp PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_BINARY_DIR}/include")
#target_link_libraries (benchmark_epr "${SEQAN_LIBRARIES} -L${CMAKE_CURRENT_BINARY_DIR}/lib -lbenchmark")
target_link_libraries (benchmark_epr "${SEQAN_LIBRARIES}")
target_link_libraries (benchmark_epr benchmark::benchmark)




#add_executable (epr epr.cpp)
#target_link_libraries (epr ${SEQAN_LIBRARIES})
#add_executable (benchmark_epr benchmark_epr.cpp)
#target_link_libraries (benchmark_epr ${SEQAN_LIBRARIES})

