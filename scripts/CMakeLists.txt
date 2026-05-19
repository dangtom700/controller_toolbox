cmake_minimum_required(VERSION 3.16)

# tune_all: outputs tuned_params.txt with optimal parameters
add_executable(tune_all tune_all.cpp)
target_link_libraries(tune_all PRIVATE controller_toolbox)
target_compile_features(tune_all PRIVATE cxx_std_17)

# simulate_all: closed-loop simulation → sim_*.csv + sim_summary.csv
add_executable(simulate_all simulate_all.cpp)
target_link_libraries(simulate_all PRIVATE controller_toolbox)
target_compile_features(simulate_all PRIVATE cxx_std_17)

# realtime_all: wall-clock-paced real-time runner → rt_*.csv + rt_summary.csv
add_executable(realtime_all realtime_all.cpp)
target_link_libraries(realtime_all PRIVATE controller_toolbox)
target_compile_features(realtime_all PRIVATE cxx_std_17)
# realtime_all uses std::this_thread::sleep_until; link pthreads on Linux/macOS
if(NOT WIN32)
    target_link_libraries(realtime_all PRIVATE pthread)
endif()
