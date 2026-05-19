cmake_minimum_required(VERSION 3.16)

add_executable(bench_controllers bench_controllers.cpp)
target_link_libraries(bench_controllers PRIVATE controller_toolbox)
target_compile_options(bench_controllers PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O3 -march=native>
    $<$<CXX_COMPILER_ID:MSVC>:/O2 /arch:AVX2>
)
