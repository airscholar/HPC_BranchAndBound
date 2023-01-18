cmake_minimum_required(VERSION 3.23)
project(HPC_BNB)
#include_directories(/opt/homebrew/Cellar/open-mpi/4.1.4_2/include)
#
## Add MPI Package to Project
#find_package(MPI REQUIRED)
#
#set(CMAKE_CXX_STANDARD 14)
#
#add_executable(HPC_BNB main.cpp Scatter.cpp)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(HPC_BNB main.cpp sample.cpp test_sample.cpp branch_n_bound_parallel.cpp adjust_sample.cpp another_sample.cpp optimised_sample.cpp branch_and_bound.cpp bnb_parallel_v1.cpp example.cpp)

target_link_libraries(HPC_BNB ${MPI_C_LIBRARIES})
