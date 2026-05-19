cmake_minimum_required(VERSION 3.16)
project(ControllerToolbox CXX)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(lib)
add_subdirectory(examples)
add_subdirectory(tests)
add_subdirectory(scripts)
add_subdirectory(case-study)