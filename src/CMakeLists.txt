
CPMAddPackage(
  NAME fortran-yaml-c
  VERSION 0.2.4
  GITHUB_REPOSITORY "Nicholaswogan/fortran-yaml-c"
  GIT_TAG "v0.2.4"
  EXCLUDE_FROM_ALL ON
)

add_library(equilibrate
  equilibrate_yaml.f90
  equilibrate_const.f90
  equilibrate_cea.f90
  equilibrate.f90
)

target_link_libraries(equilibrate
  fortran-yaml-c
)

# This is important for making Python wrapper
set_property(TARGET equilibrate PROPERTY POSITION_INDEPENDENT_CODE ON)

# Add special compiler flags to Clima
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(equilibrate PRIVATE -Wunused -Wimplicit-interface -fimplicit-none)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(equilibrate PRIVATE -fcheck=all,no-array-temps)
  endif()
endif()