cmake_minimum_required(VERSION 3.20)

project(dem2d)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# MPI
find_package(MPI REQUIRED)

set(SOURCE
  src/main.cpp
  src/particle.cpp
  src/grid.cpp
  src/contact.cpp
  src/boundary.cpp
  src/world.cpp
  src/output.cpp
)

add_executable(dem2d ${SOURCE})
target_link_libraries(dem2d MPI::MPI_CXX)
