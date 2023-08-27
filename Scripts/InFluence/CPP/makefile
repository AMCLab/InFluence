CXX = g++
CXXFLAGS = -Wall -Wextra -fPIC

# Get the include flags from pybind11
PYBIND11_INCLUDES = $(shell python3 -m pybind11 --includes)

# Set the name of the Python module
PYTHON_MODULE_NAME = InFluence

SRCS = main.cpp helper_functions.cpp electron_generator.cpp constants.cpp MCloop.cpp wrapper.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = $(PYTHON_MODULE_NAME).so

# Linking libraries
LDLIBS = -shared

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDLIBS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(PYBIND11_INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
