
add_library(equilibrate_c_api 
  equilibrate_c_api.f90 
)
target_link_libraries(equilibrate_c_api equilibrate)
