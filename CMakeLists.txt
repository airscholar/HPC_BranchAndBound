cmake_minimum_required(VERSION 3.23)
project(HPC_BNB)
include_directories(/opt/homebrew/Cellar/open-mpi/4.1.4_2/include)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(HPC_BNB wsp_bruteforce.cpp  branch_n_bound_serial.cpp branch_n_bound_parallel.cpp)

target_link_libraries(HPC_BNB ${MPI_C_LIBRARIES})
