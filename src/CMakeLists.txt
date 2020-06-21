file(GLOB sources *.f90)

add_executable(test_quaternions ${sources})
install (TARGETS test_quaternions RUNTIME DESTINATION bin)
      
