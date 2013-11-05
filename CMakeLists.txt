cmake_minimum_required(VERSION 2.6.0)	
	
## Project information ##
project(sdsl CXX C)
ADD_DEFINITIONS(
      -std=c++0x # Or -std=c++0x
      -march=native
      -Wall
      -g
       )

add_executable(fast-bit-vector_test fast-bit-vector_test.cpp fast-bit-vector.cpp)
target_link_libraries(fast-bit-vector_test
  libgtest.a
  libgtest_main.a
  pthread
  )
add_custom_command(TARGET fast-bit-vector_test POST_BUILD
    COMMAND fast-bit-vector_test || true)

add_executable(bit-vector_test bit-vector_test.cpp fast-bit-vector.cpp)
target_link_libraries(bit-vector_test
  libgtest.a
  libgtest_main.a
  pthread
  )
add_custom_command(TARGET bit-vector_test POST_BUILD
    COMMAND bit-vector_test || true)

add_executable(fast-bit-vector_benchmark fast-bit-vector_benchmark.cpp fast-bit-vector.cpp)
set_target_properties(fast-bit-vector_benchmark PROPERTIES COMPILE_FLAGS "-O3 -DNDEBUG -march=native")

add_executable(wavelet_test wavelet_test.cpp fast-bit-vector.cpp)
target_link_libraries(wavelet_test
  libgtest.a
  libgtest_main.a
  pthread
  )
add_custom_command(TARGET wavelet_test POST_BUILD
    COMMAND wavelet_test || true)

add_executable(wavelet_benchmark wavelet_benchmark.cpp fast-bit-vector.cpp)
set_target_properties(wavelet_benchmark PROPERTIES COMPILE_FLAGS "-O3 -DNDEBUG -march=native")
