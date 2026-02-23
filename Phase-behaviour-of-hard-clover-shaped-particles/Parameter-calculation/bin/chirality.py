import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Folder path containing the data files
folder_path = r'C:\Users\Nathan\Desktop\Nova pasta\carlos\OrderAndStructure\bin\Calculations'

# List of chirality file names to include
chirality_files = [
    'chirality_20250508_161505.dat'
]

# Construct full file paths
file_paths = [os.path.join(folder_path, file) for file in chirality_files]

# Colors for each file's plot
colors = ['blue', 'red', 'green', 'orange', 'purple']  # Extend this list if needed

# Loop through each file and plot its PDF
for i, file_path in enumerate(file_paths):
    # Load the data from the file
    data = np.loadtxt(file_path, delimiter=',', skiprows=1, usecols=1)
    
    # Perform kernel density estimation
    kde = gaussian_kde(data)
    x_vals = np.linspace(min(data), max(data), 1000)
    pdf_vals = kde(x_vals)
    
    # Plot the PDF
    plt.plot(x_vals, pdf_vals, color=colors[i % len(colors)], label=f'File {i + 1} PDF')

# Add plot details
plt.title('Probability Density Function')
plt.xlabel('χ')
plt.ylabel('Probability Density')
plt.grid(True)
plt.legend()

# Save the plot to a PDF file
plt.savefig(os.path.join(folder_path, 'combined_probability_density_function.pdf'))

# Show the plot
plt.show()