cmake_minimum_required(VERSION "3.14")

project(EQUILIBRATE LANGUAGES Fortran C VERSION "0.1.0")

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
include(cmake/CPM.cmake)

option(BUILD_EXECUTABLES "if ON, then will build the
Fortran executables" ON)

add_subdirectory(src)
add_subdirectory(test)