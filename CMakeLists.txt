cmake_minimum_required(VERSION 3.16)

project(bpt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_library(bpt
    src/qp_solver.cpp
    src/bicycle.cpp
    src/path.cpp
    src/ltv_mpc.cpp
    src/pure_pursuit.cpp
    src/simulator.cpp
)
target_include_directories(bpt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(bpt PRIVATE -Wall -Wextra -Wpedantic -Werror -Wshadow)

add_executable(bpt_benchmark src/benchmark.cpp)
target_link_libraries(bpt_benchmark PRIVATE bpt)
target_compile_options(bpt_benchmark PRIVATE -Wall -Wextra -Wpedantic -Werror -Wshadow)

include(FetchContent)
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(bpt_tests
    tests/test_bicycle.cpp
    tests/test_qp.cpp
    tests/test_constraints.cpp
    tests/test_closed_loop.cpp
    tests/test_path.cpp
    tests/test_pure_pursuit.cpp
)
target_link_libraries(bpt_tests PRIVATE bpt GTest::gtest_main)
target_compile_options(bpt_tests PRIVATE -Wall -Wextra -Wpedantic -Werror -Wshadow)
include(GoogleTest)
gtest_discover_tests(bpt_tests)
