cmake_minimum_required(VERSION 3.16)

add_executable(test_controllers test_controllers.cpp)
target_link_libraries(test_controllers PRIVATE controller_toolbox)

add_executable(test_tuners_extended test_tuners_extended.cpp)
target_link_libraries(test_tuners_extended PRIVATE controller_toolbox)

add_executable(test_integration test_integration.cpp)
target_link_libraries(test_integration PRIVATE controller_toolbox)
target_compile_features(test_controllers PRIVATE cxx_std_17)

# Register with CTest so "ctest" runs the tests
enable_testing()
add_test(NAME controller_tests COMMAND test_controllers)
