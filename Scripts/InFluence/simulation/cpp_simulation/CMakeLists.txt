cmake_minimum_required(VERSION 3.5)
project(InFluence)

# Set the C++ standard to C++11 or higher
set(CMAKE_CXX_STANDARD 11)

# Set the name of the Python module
set(PYTHON_MODULE_NAME "InFluence")

# Add your source files here
set(SRCS
    main.cpp
    helper_functions.cpp
    electron_generator.cpp
    constants.cpp
    MCloop.cpp
    wrapper.cpp
)

# Find the Python library and include directories
find_package(Python3 REQUIRED COMPONENTS Development)

# Add the include directories for Python headers
include_directories(${Python3_INCLUDE_DIRS})

# Create the target for your shared library (Python module)
add_library(${PYTHON_MODULE_NAME} MODULE ${SRCS})

# Link against Python library
target_link_libraries(${PYTHON_MODULE_NAME} PRIVATE Python3::Python)

# Set output name of the shared library (Python module)
set_target_properties(${PYTHON_MODULE_NAME} PROPERTIES OUTPUT_NAME ${PYTHON_MODULE_NAME})

# If you need to add additional compile options, you can add them here
# For example, to add specific compiler flags:
# target_compile_options(${PYTHON_MODULE_NAME} PRIVATE -Wall -Wextra -fPIC)

# If you need to link against additional libraries, you can add them here
# For example, to link against some external library:


# target_link_libraries(${PYTHON_MODULE_NAME} PRIVATE YourExternalLibrary)
