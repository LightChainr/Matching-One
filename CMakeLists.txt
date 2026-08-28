cmake_minimum_required(VERSION 3.16)
project(matching_one_exact_torus LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_compile_options(-O3 -march=native -Wall -Wextra -Wpedantic)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

add_executable(exact_torus_matching src/exact_torus_matching.cpp)
target_include_directories(exact_torus_matching PRIVATE src)
target_link_libraries(exact_torus_matching PRIVATE OpenMP::OpenMP_CXX)
target_compile_options(exact_torus_matching PRIVATE -fopenmp)

add_executable(test_torus_connectivity tests/test_torus_connectivity.cpp)
target_include_directories(test_torus_connectivity PRIVATE src)
target_link_libraries(test_torus_connectivity PRIVATE OpenMP::OpenMP_CXX)
target_compile_options(test_torus_connectivity PRIVATE -fopenmp)

enable_testing()
add_test(NAME torus_connectivity COMMAND test_torus_connectivity)

add_test(
  NAME exact_L2_L3
  COMMAND exact_torus_matching --L 2 --L 3 --threads 8 --outdir ${CMAKE_SOURCE_DIR}/results/issue-7 --dump-configs
)
set_tests_properties(exact_L2_L3 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=8;OMP_PROC_BIND=close;OMP_PLACES=cores")

add_test(
  NAME exact_L4
  COMMAND exact_torus_matching --L 4 --threads 8 --outdir ${CMAKE_SOURCE_DIR}/results/issue-7
)
set_tests_properties(exact_L4 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=8;OMP_PROC_BIND=close;OMP_PLACES=cores")

add_test(
  NAME exact_L5
  COMMAND exact_torus_matching --L 5 --threads 8 --outdir ${CMAKE_SOURCE_DIR}/results/issue-7
)
set_tests_properties(exact_L5 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=8;OMP_PROC_BIND=close;OMP_PLACES=cores")
