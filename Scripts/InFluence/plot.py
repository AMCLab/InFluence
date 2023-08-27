

import matplotlib.pyplot as plt

# Sample data (you can replace these with your own data)
x = [1, 2, 3, 4, 5]
y = [2, 4, 6, 8, 10]

# Create the plot
plt.plot(x, y, marker='o', linestyle='-', color='b', label='Data points')

# Add labels and title
plt.xlabel('x')
plt.ylabel('y')
plt.title('Plot of y vs. x')

# Add a legend
plt.legend()

# Show the plot
plt.grid(True)
plt.show()