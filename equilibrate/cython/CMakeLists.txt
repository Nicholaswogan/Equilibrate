
add_cython_target(_equilibrate _equilibrate.pyx)

if(WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMS_WIN64")
endif()

add_library(_equilibrate MODULE ${_equilibrate})

python_extension_module(_equilibrate)

target_link_libraries(_equilibrate equilibrate_c_api)
target_include_directories(_equilibrate PUBLIC ${NumPy_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET _equilibrate PROPERTY LINKER_LANGUAGE Fortran)