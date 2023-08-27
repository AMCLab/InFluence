import math
import timeit
import numpy as np

def evaluate_cross_section(E, Z):
    u = math.log10(8 * E * Z ** -1.33)
    return (4.7 * 10 ** -18) * (Z ** 1.33 + 0.032 * Z ** 2) / (
                (E + (E ** 0.5) * 0.0155 * Z ** 1.33) * (1 - math.exp(-1*u**2) * 0.02 * Z ** 0.5))

def constants(Z):
    CrossSectionNumorator = (4.7 * 10 ** -18) * (Z ** 1.33 + 0.032 * Z ** 2)
    CrossSectionLogArgMultiplier = 8 * Z ** -1.33
    CrossSectionDenominatorA = 0.0155 * Z ** 1.33
    CrossSectionDenominatorB = 0.02 * Z ** 0.5
    return CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB

def evaluate_cross_section_opt(CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, E):
    LogArg = math.log10(E * CrossSectionLogArgMultiplier)
    CrossSection = CrossSectionNumorator / ((E + (E ** 0.5) * CrossSectionDenominatorA) * (1 - math.exp(-1 * LogArg ** 2) * CrossSectionDenominatorB))
    return CrossSection

Z = 16
CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB = constants(Z)



# Additional functions
def retrieve_pi():
    val = np.sin(0.25)

local_sin = np.sin(0.25)  # This should be outside the function, so it's pre-computed.
def use_assigned_pi():
    val = local_sin

# Benchmarking for original and optimized versions
N_values = [1000, 10000, 100000, 1000000]
original_times = []
opt_times = []

for N in N_values:
    orig_time = timeit.timeit(lambda: evaluate_cross_section(200, Z), number=N)
    opt_time = timeit.timeit(lambda: evaluate_cross_section_opt(CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, 200), number=N)

    original_times.append(orig_time)
    opt_times.append(opt_time)

# Print the results in table format for original and optimized versions
print("N-values    | Evaluate Cross Section Function | Opt Function")
print("------------|---------------------------------|--------------")
for N, orig_time, opt_time in zip(N_values, original_times, opt_times):
    print(f"{N:<12}| {orig_time:<17.5f} | {opt_time:.5f}")
print("\n")

# Benchmarking for retrieve_pi and use_assigned_pi functions
retrieve_times = []
assigned_times = []

for N in N_values:
    retrieve_time = timeit.timeit(retrieve_pi, number=N)
    assigned_time = timeit.timeit(use_assigned_pi, number=N)

    retrieve_times.append(retrieve_time)
    assigned_times.append(assigned_time)

# Print the results in table format for retrieve_pi and use_assigned_pi
print("N-values    | Retrieve Pi       | Use Local Pi")
print("------------|-------------------|---------------")
for N, retrieve_time, assigned_time in zip(N_values, retrieve_times, assigned_times):
    print(f"{N:<12}| {retrieve_time:<17.5f} | {assigned_time:.5f}")


local_sin = math.sin
local_cos = math.cos
local_sqrt = math.sqrt

# Your new function using local_sin and local_cos
def evaluate_direction_cosine_a_local(phi, psi, cosineX, cosineY, cosineZ):
    alpha = local_sin(psi) * local_sin(phi)
    beta = local_sin(phi) * local_cos(psi)
    gamma = local_cos(phi)
    cos_1 = cosineZ
    sin_1 = local_sqrt(1 - cosineZ**2)
    cos_2 = cosineY / sin_1
    sin_2 = cosineX / sin_1
    return alpha * cos_2 + sin_2 * (beta * cos_1 + gamma * sin_1)

# Your new function using numpy
def evaluate_direction_cosine_a_numpy(phi, psi, cosineX, cosineY, cosineZ):
    alpha = math.sin(psi) * math.sin(phi)
    beta = math.sin(phi) * math.cos(psi)
    gamma = math.cos(phi)
    cos_1 = cosineZ
    sin_1 = math.sqrt(1 - cosineZ**2)
    cos_2 = cosineY / sin_1
    sin_2 = cosineX / sin_1
    return alpha * cos_2 + sin_2 * (beta * cos_1 + gamma * sin_1)

# Benchmarking for the direction cosine functions
local_times = []
numpy_times = []

for N in N_values:
    local_time = timeit.timeit(lambda: evaluate_direction_cosine_a_local(0.2, 0.2, 0.2, 0.2, 0.2), number=N)
    numpy_time = timeit.timeit(lambda: evaluate_direction_cosine_a_numpy(0.2, 0.2, 0.2, 0.2, 0.2), number=N)
    
    local_times.append(local_time)
    numpy_times.append(numpy_time)

# Print the results in table format for the direction cosine functions
print("N-values    | Evaluate Diretion Cosine Function with Local Math Function  | Evaluate Diretion Cosine Function with Math Function")
print("------------|-------------------------------------------------------------|------------------------------------------------------")
for N, local_time, numpy_time in zip(N_values, local_times, numpy_times):
    print(f"{N:<12}| {local_time:<17.5f} | {numpy_time:.5f}")
