# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (EPR-SeqAn2 CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

include (../epr-test.cmake)

epr_require_benchmark()

# ----------------------------------------------------------------------------
# Add targets
# ----------------------------------------------------------------------------

# Add executable and link against SeqAn dependencies.
add_executable (benchmark_epr benchmark_epr.cpp)
target_link_libraries (benchmark_epr seqan::epr::test::performance)

#add_executable (epr epr.cpp)
#target_link_libraries (epr ${SEQAN_LIBRARIES})
#add_executable (benchmark_epr benchmark_epr.cpp)
#target_link_libraries (benchmark_epr ${SEQAN_LIBRARIES})

