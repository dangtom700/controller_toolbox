cmake_minimum_required(VERSION 3.16)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

add_library(controller_toolbox STATIC
    PlantModel.cpp
    DiscretePID.cpp
    DiscreteLQR.cpp
    DiscreteMPC.cpp
    ExtremumSeeker.cpp
    SmithPredictor.cpp
    DiscreteLeadLag.cpp
    DiscreteSMC.cpp
    DiscreteADRC.cpp
    KalmanFilter.cpp
    DiscreteLQG.cpp
    ControllerTuner.cpp
    ControllerStack.cpp
    TunerSuite.cpp
    MetricsAnalyzer.cpp
    SystemAnalysis.cpp
    # HAL simulation adapters
    hal/SimPlant.cpp
)

# Expose lib/ as the public include root so consumers write:
#   #include "ControllerToolbox.h"
target_include_directories(controller_toolbox PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(controller_toolbox PUBLIC Eigen3::Eigen)
target_compile_features(controller_toolbox PUBLIC cxx_std_17)
