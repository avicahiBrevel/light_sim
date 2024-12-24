import numpy as np
import matplotlib.pyplot as plt

# Define the grid size and resolution
x_min, x_max = -10, 10
y_min, y_max = -10, 10
resolution = 500

# Create a grid of points
x = np.linspace(x_min, x_max, resolution)
y = np.linspace(y_min, y_max, resolution)
x_grid, y_grid = np.meshgrid(x, y)

# Define the positions of the two point lights
light1 = np.array([-2, 0])
light2 = np.array([2, 0])

# Define a function to calculate light intensity at a point
def light_intensity(point, light_pos, intensity=1):
    distance = np.sqrt((point[0] - light_pos[0])**2 + (point[1] - light_pos[1])**2)
    return intensity / (distance**2 + 1e-6)  # Avoid division by zero

# Calculate the total intensity at each point on the grid
intensity_grid = (
    light_intensity((x_grid, y_grid), light1) +
    light_intensity((x_grid, y_grid), light2)
)

# Plot the color map
plt.figure(figsize=(8, 6))
plt.imshow(
    intensity_grid, 
    extent=[x_min, x_max, y_min, y_max], 
    origin='lower', 
    cmap='plasma',
    interpolation='bilinear'
)

plt.colorbar(label='Light Intensity')
plt.scatter([light1[0], light2[0]], [light1[1], light2[1]], color='white', label='Light Sources')
plt.title('Light Intensity Map for Two Adjacent Point Lights')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.legend()
plt.grid(False)
plt.show()