
set(TESTS test_equilibrate)

if (BUILD_EXECUTABLES)
  foreach(TEST ${TESTS})
    add_executable(${TEST} "${TEST}.f90")
    target_link_libraries(${TEST} equilibrate)
    target_include_directories(${TEST} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
  endforeach()
endif()