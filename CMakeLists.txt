cmake_minimum_required(VERSION 2.6.0)	
	
## Project information ##
project(sdsl CXX C)
ADD_DEFINITIONS(
      -std=c++0x # Or -std=c++0x
      -march=native
      -g
       )

add_executable(fastbitvector_test fastbitvector_test.cpp fastbitvector.cpp)
target_link_libraries(fastbitvector_test
  libgtest.a
  libgtest_main.a
  pthread
  )
add_custom_command(TARGET fastbitvector_test POST_BUILD
    COMMAND fastbitvector_test || true)

add_executable(fastbitvector_benchmark fastbitvector_benchmark.cpp fastbitvector.cpp)
set_target_properties(fastbitvector_benchmark PROPERTIES COMPILE_FLAGS "-O3 -DNDEBUG")
