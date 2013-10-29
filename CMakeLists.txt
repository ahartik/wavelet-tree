cmake_minimum_required(VERSION 2.6.0)	
	
## Project information ##
project(sdsl CXX C)
ADD_DEFINITIONS(
      -std=c++0x # Or -std=c++0x
      -g
       )

add_executable(test test.cpp)
target_link_libraries(test
  libgtest.a
  libgtest_main.a
  pthread
  )

add_custom_command(TARGET test POST_BUILD COMMAND test || true)

ADD_DEFINITIONS(
      -march=native
      -O3
       )

add_executable(benchmark benchmark.cpp)
