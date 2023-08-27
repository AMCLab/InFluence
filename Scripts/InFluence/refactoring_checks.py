import math
import timeit




#Original
def evaluate_cross_section(E, Z):

    u = math.log10(8 * E * Z ** -1.33)
    return (4.7 * 10 ** -18) * (Z ** 1.33 + 0.032 * Z ** 2) / (
                (E + (E ** 0.5) * 0.0155 * Z ** 1.33) * (1 - math.exp(-1*u**2) * 0.02 * Z ** 0.5))


#Optimized
def constants(Z):

    CrossSectionNumorator = (4.7 * 10 ** -18) * (Z ** 1.33 + 0.032 * Z ** 2)
    CrossSectionLogArgMultiplier = 8*Z ** -1.33
    CrossSectionDenominatorA = 0.0155 * Z ** 1.33
    CrossSectionDenominatorB = 0.02 * Z ** 0.5

    return CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB

    
def evaluate_cross_section_opt(CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB):
    
    LogArg = math.log10(E*CrossSectionLogArgMultiplier)

    CrossSection = CrossSectionNumorator / ((E + (E ** 0.5) * CrossSectionDenominatorA) * (1 - math.exp(-1*LogArg**2) * CrossSectionDenominatorB))

    return CrossSection


Z = 2

CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB = constants(Z)




for E in range(1,100):
    original_result = evaluate_cross_section(E, Z)
    optimized_result = evaluate_cross_section_opt(CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB)
    assert math.isclose(original_result, optimized_result, rel_tol=1e-9)






import math
import timeit

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


Z = 2
CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB = constants(Z)

# Benchmarking the original version
original_time = timeit.timeit(lambda: evaluate_cross_section(1000, Z), number=10000000)

# Benchmarking the optimized version
optimized_time = timeit.timeit(lambda: evaluate_cross_section_opt(CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, 1000), number=10000000)

print("Original Version Time:", original_time)
print("Optimized Version Time:", optimized_time)

local_sin = np.sin

def retrieve_pi():
    val = np.sin(0.25)
    

def use_assigned_pi():
    val = local_sin

retrieve_pi_time = timeit.timeit(retrieve_pi, number=1)
use_assigned_pi_time = timeit.timeit(use_assigned_pi, number=1)

print("Retrieve math.pi inside loop:", retrieve_pi_time)
print("Use assigned pi outside loop:", use_assigned_pi_time)




#######################################################################################################################



