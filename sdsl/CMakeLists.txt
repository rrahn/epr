project(EPR-SDSL CXX)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)

include (../epr-test.cmake)
include (GoogleTest OPTIONAL)

epr_require_test()
epr_require_benchmark()

include_directories("${CMAKE_SOURCE_DIR}/../include")

add_executable(rank_test rank_test.cpp)
target_link_libraries (rank_test epr::test::unit)
add_test(NAME rank_test_test COMMAND rank_test)

add_executable(rank_benchmark rank_benchmark.cpp)
target_link_libraries(rank_benchmark epr::test::performance)

add_executable(rank_profile rank_profile.cpp)
target_link_libraries(rank_profile epr::test::performance)

add_executable(epr_benchmark epr_benchmark.cpp)
target_link_libraries(epr_benchmark epr::test::performance)

add_executable(access_benchmark access_benchmark.cpp)
target_link_libraries(access_benchmark epr::test::performance)
