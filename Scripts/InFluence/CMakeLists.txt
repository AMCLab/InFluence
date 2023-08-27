cmake_minimum_required(VERSION 3.12)
project(InFluence VERSION 1.0 LANGUAGES CXX)

# Find Python and pybind11
find_package(pybind11 REQUIRED)

# Add your C++ executable or library source files
add_library(InFluence MODULE main.cpp helper_functions.cpp electron_generator.cpp constants.cpp MCloop.cpp)

# Link the C++ code with pybind11 and Python
target_link_libraries(InFluence PRIVATE pybind11::module)

# Set the C++ standard (e.g., C++11)
target_compile_features(InFluence PRIVATE cxx_std_11)
