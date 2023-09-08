import random
import timeit

# Setting up the lists
list_1d = [i for i in range(100)]
list_2d = [[i, i + 1] for i in range(0, 100, 2)]  # This creates pairs of numbers in sublists

def test_1d_list():
    item = random.randint(0, 9999)
    if item in list_1d:
        pass

def test_2d_list():
    sublist = [random.randint(0, 9998), random.randint(0, 9998)]  # Adjusted to create pairs
    if sublist in list_2d:
        pass

Ns = [1000, 10000, 100000, 1000000]
results = []

# Loop through values of N and measure time for both 1D and 2D lists
for N in Ns:
    time_1d = timeit.timeit(test_1d_list, number=N)
    time_2d = timeit.timeit(test_2d_list, number=N)
    results.append((N, time_1d, time_2d))

# Printing the results in table format
print("N        | 1D List Time | 2D List Time")
print("--------------------------------------")
for N, time_1d, time_2d in results:
    print(f"{N:<8} | {time_1d:<12.6f} | {time_2d:.6f}")
