#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Include this header for automatic conversion

#include "MCloop.h"

namespace py = pybind11;

// Wrap the RunMCScatteringSimulation function
py::list runMCScatteringSimulationWrapper(
    py::list pixels,
    double E_i,
    double ProbeDiameter,
    double MinimumEnergy,
    double dE_threshold,
    int perfect_image_0,
    int perfect_image_1,
    double Density,
    double t_counting
) {
    // Convert Python list of lists to C++ vector of vectors explicitly
    std::vector<std::vector<double>> pixels_cpp;
    for (const auto& row : pixels) {
        // Convert each row (Python list) to a C++ vector explicitly
        pixels_cpp.emplace_back(py::cast<std::vector<double>>(row));
    }

    // Call the original C++ function
    std::vector<std::vector<double>> result = RunMCScatteringSimulation(
        pixels_cpp,
        E_i,
        ProbeDiameter,
        MinimumEnergy,
        dE_threshold,
        perfect_image_0,
        perfect_image_1,
        Density,
        t_counting
    );

    // Convert the result back to a Python list of lists
    py::list result_py;
    for (const auto& row : result) {
        result_py.append(py::cast(row));
    }

    return result_py;
}

PYBIND11_MODULE(InFluence, m) {
    m.doc() = "Your module's docstring";  // Optional module docstring

    // Bind the runMCScatteringSimulationWrapper function
    m.def("RunMCScatteringSimulation", &runMCScatteringSimulationWrapper, "A description of your function");
}
