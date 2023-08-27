import Scripts.InFluence.HelperFunctions as new
import LegacyHelperFunctions as legacy
from parameters import params

import pytest
import numpy as np
from numba import jit
import math



def test_rounding_functions():

    test_inputs = [(2.5, 0), (10.4, 0), (-10.5555, 0), (-10.49, 0), (0.545, 2), (10.555, 2)]

    for x, decimal in test_inputs:

        result_1 = new.custom_round_down(x, decimal)
        result_2 = legacy.round_half_down(x, decimal)


        assert result_1 == result_2, f"For input ({x}, {decimal}), custom_round_down got {result_1} but round_half_down got {result_2}"



def test_evaluate_alpha():


    test_inputs = [(10, 2), (20, 3), (30, 4), (40, 5), (1e-6, 2e-6), (1e-7, 2e-7), (1e-8, 2e-18)]

    for E, Z in test_inputs:

        AlphaMultiplier = (3.4*10**-3)*(Z**0.67) 

        result_1 = new.evaluate_alpha(E, AlphaMultiplier)
        result_2 = legacy.evaluate_alpha(E, Z)

        assert result_1 == result_2, f"For input ({E}, {Z}), evaluate_alpha in new got {result_1} but evaluate_alpha in legacy got {result_2}"


def test_evaluate_path_length():


    test_inputs = [(10, 2, 3, 4), (20, 3, 5, 6), (30, 4, 7, 8), (40, 5, 9, 10), (1e-6, 2e-6, 3e-6, 4e-6), (1e-7, 2e-7, 5e-7, 6e-7), (1e-8, 2e-8, 7e-8, 8e-8)]

    for A, N, p, cross in test_inputs:


        PathLengthMultiplier = A/(N * p)

        result_1 = new.evaluate_path_length(cross, PathLengthMultiplier)
        result_2 = legacy.evaluate_path_length(A, N, p, cross)

      
        assert result_1 == result_2, f"For input ({A}, {N}, {p}, {cross}), evaluate_path_length in new got {result_1} but evaluate_path_length in legacy got {result_2}"


def test_evaluate_phi():

    test_inputs = [(0.5, 2), (0.4, 3), (0.3, 4), (0.2, 5), (1e-7, 2), (1e-8, 3), (1e-9, 4), (1e-10, 5)]

    for RandomNum, alpha in test_inputs:

        result_1 = new.evaluate_phi(RandomNum, alpha)
        result_2 = legacy.evaluate_phi(RandomNum, alpha)

        assert result_1 == result_2, f"For input ({RandomNum}, {alpha}), evaluate_phi in new got {result_1} but evaluate_phi in legacy got {result_2}"


def test_evaluate_pho():

    test_inputs = [0.1, 0.2, 0.3, 0.4, 0.5, 1e-7, 1e-8, 1e-9, 1e-10]

    for RandomNum in test_inputs:

        result_1 = new.evaluate_pho(RandomNum, np.pi)
        result_2 = legacy.evaluate_pho(RandomNum)

        assert result_1 == result_2, f"For input ({RandomNum}), evaluate_pho in new got {result_1} but evaluate_pho in legacy got {result_2}"

def test_evaluate_step():

    test_inputs = [(1, 0.5), (2, 0.4), (3, 0.3), (4, 0.2), (5, 0.1), (10, 1e-7), (20, 1e-8), (30, 1e-9), (40, 1e-10)]

    for path_length, RandomStep in test_inputs:

        result_1 = new.evaluate_step(path_length, RandomStep)
        result_2 = legacy.evaluate_step(path_length, RandomStep)

        assert result_1 == result_2, f"For input ({path_length}, {RandomStep}), evaluate_step in new got {result_1} but evaluate_step in legacy got {result_2}"

def test_evaluate_cross_section():

    test_inputs = [(10, 2), (20, 3), (30, 4), (40, 5), (50, 6), (100, 7), (1000, 8), (1e-5, 1), (1e-6, 2), (1e-7, 3)]

    for E, Z in test_inputs:
        CrossSectionNumorator = (4.7 * 10 ** -18) * (Z ** 1.33 + 0.032 * Z  ** 2)
        CrossSectionLogArgMultiplier = 8 * Z  ** -1.33
        CrossSectionDenominatorA = 0.0155 * Z ** 1.33
        CrossSectionDenominatorB = 0.02 * Z  ** 0.5

        result_1 = new.evaluate_cross_section_opt(E, CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB)
        result_2 = legacy.evaluate_cross_section(E, Z)

        assert np.isclose(result_1, result_2, atol=1e-10), f"For input ({E}, {Z}), evaluate_cross_section_opt in new got {result_1} but evaluate_cross_section in legacy got {result_2}"


def test_evaluate_direction_cosine_a():

    test_inputs = [(0.1, 0.2, 0.3, 0.4, 0.5),
                (0.5, 0.4, 0.3, 0.2, 0.8),  # cz or cosineZ is 0.8, not 1.0
                (0.6, 0.7, 0.8, 0.9, 0.99),  # cz or cosineZ is 0.99, not 1.0
                (1e-6, 1e-7, 1e-8, 1e-9, 0.5),  # cz or cosineZ is 0.5, not 1.0
                (np.pi, np.pi/2, np.pi/3, np.pi/4, np.pi/5)]  # cz or cosineZ is np.pi/5, not 1.0

    for phi, psi, cx, cy, cz in test_inputs:

        # Run the functions with the test inputs.
        result_1 = new.evaluate_direction_cosine_a(phi, psi, cx, cy, cz)
        result_2 = legacy.evaluate_direction_cosine_a(phi, psi, cx, cy, cz)

        # Check if the results are equal to within a tolerance.
        assert math.isclose(result_1, result_2, rel_tol=1e-9), f"For input ({phi}, {psi}, {cx}, {cy}, {cz}), evaluate_direction_cosine_a in new got {result_1} but evaluate_direction_cosine_a in legacy got {result_2}"


def test_evaluate_direction_cosine_b():
    # Generating random sets of inputs
    test_inputs = [(0.1, 0.2, 0.3, 0.4, 0.5), 
                (0.6, 0.7, 0.8, 0.9, 0.8), 
                (0.005, 0.01, 0.015, 0.02, 0.99),
                (0.0001, 0.0002, 0.0003, 0.0004, 0.5)]

    for phi, psi, cx, cy, cz in test_inputs:
        result_1 = new.evaluate_direction_cosine_b(phi, psi, cx, cy, cz)
        result_2 = legacy.evaluate_direction_cosine_b(phi, psi, cx, cy, cz)

        assert math.isclose(result_1, result_2, rel_tol=1e-9), f"For input ({phi}, {psi}, {cx}, {cy}, {cz}), evaluate_direction_cosine_b in new got {result_1} but evaluate_direction_cosine_b in legacy got {result_2}"

def test_evaluate_direction_cosine_c():

    test_inputs = [(0.1, 0.2, 0.3), 
                   (0.6, 0.7, 0.8), 
                   (0.005, 0.01, 0.02),
                   (0.0001, 0.0002, 0.0003)]

    for phi, psi, cz in test_inputs:
        result_1 = new.evaluate_direction_cosine_c(phi, psi, cz)
        result_2 = legacy.evaluate_direction_cosine_c(phi, psi, cz)

        assert math.isclose(result_1, result_2, rel_tol=1e-9), f"For input ({phi}, {psi}, {cz}), evaluate_direction_cosine_c in new got {result_1} but evaluate_direction_cosine_c in legacy got {result_2}"


def test_evaluate_energy_loss_rate():
    # Generating random sets of inputs
    test_inputs = [(0.1, 0.2, 0.3), 
                   (0.6, 0.7, 0.8), 
                   (0.005, 0.01, 0.02),
                   (0.0001, 0.0002, 0.0003)]

    for E, Z, A in test_inputs:
        # Compute EnergyLossMultiplierA and EnergyLossMultiplierB for the new function
        ProtonNum = Z
        AtomicMass = A
        EnergyLossMultiplierA = -78500*ProtonNum/AtomicMass
        EnergyLossMultiplierB = (9.76*ProtonNum + 58.5/ProtonNum**0.19)*10**-3

        result_1 = new.evaluate_energy_loss_rate(E, EnergyLossMultiplierA, EnergyLossMultiplierB)
        result_2 = legacy.evaluate_energy_loss_rate(E, Z, A)

        assert math.isclose(result_1, result_2, rel_tol=1e-9), f"For input ({E}, {Z}, {A}), evaluate_energy_loss_rate in new got {result_1} but evaluate_energy_loss_rate in legacy got {result_2}"



