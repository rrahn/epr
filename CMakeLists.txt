project(EPR CXX)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)

enable_testing()
find_package(GTest REQUIRED)
find_package(benchmark REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories("${CMAKE_SOURCE_DIR}/include")

add_executable(rank_test rank_test.cpp)
target_link_libraries(rank_test ${GTEST_MAIN_LIBRARIES})
target_link_libraries(rank_test -lgtest)
target_link_libraries(rank_test -pthread)
add_test(NAME rank_test_test COMMAND rank_test)


if (CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif (CMAKE_BUILD_TYPE MATCHES Release)

add_executable(rank_benchmark rank_benchmark.cpp)
target_link_libraries(rank_benchmark benchmark::benchmark)

add_executable(epr_benchmark epr_benchmark.cpp)
target_link_libraries(epr_benchmark benchmark::benchmark)

add_executable(access_benchmark access_benchmark.cpp)
target_link_libraries(access_benchmark benchmark::benchmark)
