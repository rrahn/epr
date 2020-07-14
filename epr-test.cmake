cmake_minimum_required(VERSION 3.16)

# Include the sdsl from the submodule.
find_path (EPR_SUBMODULES_DIR NAMES submodule/sdsl-lite HINTS "${CMAKE_SOURCE_DIR}/..")

message (STATUS "Looking for submodules")
# Find all submodules and add the include directories.
if (EPR_SUBMODULES_DIR)
    set (SDSL_LITE_SUBMODULE_DIR "${EPR_SUBMODULES_DIR}/submodule/sdsl-lite")
    message(STATUS "   ... adding sdsl submodule: ${SDSL_LITE_SUBMODULE_DIR}")
else ()
    message (FATAL_ERROR "Could not find the submodule directory.")
endif ()

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)
find_package (BZip2)

# Load the SeqAn module and fail if not found.
set (SEQAN_BASE_DIRECTORY "/Users/rrahn/Development/seqan/seqan") # Set to the path on your system
set (SEQAN_INCLUDE_PATH "${SEQAN_BASE_DIRECTORY}/include")
find_package (OpenMP REQUIRED)
find_package (SeqAn REQUIRED PATHS "${SEQAN_BASE_DIRECTORY}/util/" NO_CMAKE_PATH)

set (EPR_BENCHMARK_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/benchmark")
set (EPR_TEST_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/googletest")

# needed for add_library (seqan3::test::* INTERFACE IMPORTED)
# see cmake bug https://gitlab.kitware.com/cmake/cmake/issues/15052
file(MAKE_DIRECTORY ${EPR_BENCHMARK_CLONE_DIR}/include/)
file(MAKE_DIRECTORY ${EPR_TEST_CLONE_DIR}/googletest/include/)

# ----------------------------------------------------------------------------
# Interface targets for the different test modules in seqan3.
# ----------------------------------------------------------------------------

# epr::test exposes a base set of required flags, includes, definitions and
# libraries which are in common for **all** epr tests
add_library (epr_test INTERFACE)
target_compile_options (epr_test INTERFACE "")
target_link_libraries (epr_test INTERFACE "pthread")
target_include_directories (epr_test INTERFACE "${EPR_TEST_INCLUDE_DIR}" "${SDSL_LITE_SUBMODULE_DIR}/include")
add_library (epr::test ALIAS epr_test)

# epr::test::performance specifies required flags, includes and libraries
# needed for performance test cases
add_library (epr_test_performance INTERFACE)
target_link_libraries (epr_test_performance INTERFACE "epr::test" "gbenchmark")
target_include_directories (epr_test_performance INTERFACE "epr::test" "${EPR_BENCHMARK_CLONE_DIR}/include/")
add_library (epr::test::performance ALIAS epr_test_performance)

# epr::test::unit specifies required flags, includes and libraries
# needed for the unit test cases
add_library (epr_test_unit INTERFACE)
target_link_libraries (epr_test_unit INTERFACE "epr::test" "gtest_main" "gtest")
target_include_directories (epr_test_unit INTERFACE "epr::test" "${EPR_TEST_CLONE_DIR}/googletest/include/")
add_library (epr::test::unit ALIAS epr_test_unit)

# SeqAn build system might add a leading whitespace wich confuses cmake and escaps the content in "" which breaks the compile command.
# seqan::epr::test::performance specifies required flags, includes and libraries
# needed for the SeqAn peformance test cases
string(STRIP "${SEQAN_CXX_FLAGS}" SEQAN_CXX_FLAGS_)
add_library (seqan_epr_test_performance INTERFACE)
target_compile_options (seqan_epr_test_performance INTERFACE "-DSEQAN_DISABLE_VERSION_CHECK=YES" "${SEQAN_DEFINITIONS}" "${SEQAN_CXX_FLAGS_}")
target_link_libraries (seqan_epr_test_performance INTERFACE "epr::test::performance" "${SEQAN_LIBRARIES}")
target_include_directories (seqan_epr_test_performance INTERFACE "epr::test::performance" "${SEQAN_INCLUDE_DIRS}")
add_library (seqan::epr::test::performance ALIAS seqan_epr_test_performance)

# ----------------------------------------------------------------------------
# Commonly shared options for external projects.
# ----------------------------------------------------------------------------

set (EPR_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND EPR_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND EPR_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND EPR_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND EPR_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

macro (epr_require_benchmark)
    enable_testing ()

    set (gbenchmark_project_args ${EPR_EXTERNAL_PROJECT_CMAKE_ARGS})
    list (APPEND gbenchmark_project_args "-DBENCHMARK_ENABLE_TESTING=false")
    # list (APPEND gbenchmark_project_args "-DBENCHMARK_ENABLE_LTO=true")

    # force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`
    list (APPEND gbenchmark_project_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

    set (gbenchmark_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}benchmark${CMAKE_STATIC_LIBRARY_SUFFIX}")

    include (ExternalProject)
    ExternalProject_Add (
        gbenchmark_project
        PREFIX gbenchmark_project
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        GIT_TAG "v1.5.0"
        SOURCE_DIR "${EPR_BENCHMARK_CLONE_DIR}"
        CMAKE_ARGS "${gbenchmark_project_args}"
        BUILD_BYPRODUCTS "${gbenchmark_path}"
    )
    unset (gbenchmark_project_args)

    add_library (gbenchmark STATIC IMPORTED)
    add_dependencies (gbenchmark gbenchmark_project)
    set_target_properties (gbenchmark PROPERTIES IMPORTED_LOCATION "${gbenchmark_path}")
    set_property (TARGET gbenchmark APPEND PROPERTY INTERFACE_LINK_LIBRARIES "pthread")

    unset (gbenchmark_path)
endmacro ()

macro (epr_require_test)
    enable_testing ()

    set (gtest_project_args ${EPR_EXTERNAL_PROJECT_CMAKE_ARGS})
    list (APPEND gtest_project_args "-DBUILD_GMOCK=0")

    # force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`
    list (APPEND gtest_project_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

    # google sets CMAKE_DEBUG_POSTFIX = "d"
    set (gtest_main_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (gtest_main_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_maind${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif ()

    set (gtest_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (gtest_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtestd${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif ()

    include (ExternalProject)
    ExternalProject_Add (
        gtest_project
        PREFIX gtest_project
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        # we currently have warnings that were introduced in
        # 03867b5389516a0f185af52672cf5472fa0c159c, which are still available
        # in "release-1.8.1", see https://github.com/google/googletest/issues/1419
        GIT_TAG "release-1.10.0"
        SOURCE_DIR "${EPR_TEST_CLONE_DIR}"
        CMAKE_ARGS "${gtest_project_args}"
        BUILD_BYPRODUCTS "${gtest_main_path}" "${gtest_path}"
    )
    unset (gtest_project_args)

    add_library (gtest_main STATIC IMPORTED)
    add_dependencies (gtest_main gtest_project)
    set_target_properties (gtest_main PROPERTIES IMPORTED_LOCATION "${gtest_main_path}")

    add_library (gtest STATIC IMPORTED)
    add_dependencies (gtest gtest_main)
    set_target_properties (gtest PROPERTIES IMPORTED_LOCATION "${gtest_path}")
    set_property (TARGET gtest APPEND PROPERTY INTERFACE_LINK_LIBRARIES "pthread")

    unset(gtest_main_path)
    unset(gtest_path)
endmacro ()
