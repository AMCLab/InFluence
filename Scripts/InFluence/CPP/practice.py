import InFluence

# Example usage
pixels = [
    [1, 1, 1, 0.1, 0.1, 0.1],
    # Add more pixel data as needed
]

E_i = 100.0  # Example initial energy value
ProbeDiameter = 0.1  # Example ProbeDiameter value
MinimumEnergy = 1.0  # Example MinimumEnergy value
dE_threshold = 0.1  # Example dE_threshold value
perfect_image_0 = 10  # Example perfect_image_0 value
perfect_image_1 = 10  # Example perfect_image_1 value
Density = 2.0  # Example Density value
t_counting = 100.0  # Example t_counting value

result = InFluence.RunMCScatteringSimulation(
    pixels, E_i, ProbeDiameter, MinimumEnergy, dE_threshold,
    perfect_image_0, perfect_image_1, Density, t_counting
)

# The 'result' variable now holds the output from your C++ function
