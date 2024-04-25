
set(TESTS test_chemequi)

if (BUILD_EXECUTABLES)
  foreach(TEST ${TESTS})
    add_executable(${TEST} "${TEST}.f90")
    target_link_libraries(${TEST} equichem)
    target_include_directories(${TEST} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
  endforeach()
endif()