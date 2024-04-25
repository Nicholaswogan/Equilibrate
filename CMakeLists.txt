cmake_minimum_required(VERSION "3.14")

project(EQUICHEM LANGUAGES Fortran C VERSION "0.1.0")

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
include(cmake/CPM.cmake)

add_subdirectory(src)