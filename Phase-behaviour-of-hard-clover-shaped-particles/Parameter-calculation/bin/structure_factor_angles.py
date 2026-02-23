import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import colormaps
import os

# Get the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

file_name = 'structure_factor_20251122_200446'
file_path = os.path.join(current_dir, file_name + '.dat')

# Create the 2D dictionary with indices from -600 to 600
Z = {x: {} for x in range(-600, 601)}

with open(file_path, 'r') as f:
    header = f.readline()  # skip "x,y,S(k)"

    for line in f:
        if not line.strip():
            continue

        parts = line.split(',')
        if len(parts) < 3:
            continue

        x = int(parts[0])
        y = int(parts[1])
        z = float(parts[2])

        # Store z directly using x and y as dictionary keys
        Z[x][y] = z

for x in Z:
    for y in Z[x]:
        Z[x][y] = max(Z[x][y] ** 3.5, 0.01)

def max_in_window(Z, xmin, xmax, ymin, ymax):
    max_value = 0.0
    max_x = None
    max_y = None

    # x varies in the given window
    for x in range(xmin, xmax + 1):
        # y varies in the given window
        for y in range(ymin, ymax + 1):
            value = Z[x][y]
            if value > max_value:
                max_value = value
                max_x = x
                max_y = y

    return max_value, max_x, max_y

def max_in_window(Z, xmin, xmax, ymin, ymax):
    max_value = -1.0
    max_x = None
    max_y = None

    for x in range(xmin, xmax + 1):
        for y in range(ymin, ymax + 1):
            value = Z[x][y]
            if value > max_value:
                max_value = value
                max_x = x
                max_y = y

    return max_value, max_x, max_y

windows = [
    (-100, -50, 100, 200),
    (60, 180, 90, 190),
    (120, 220, -50, 50),
    (0, 120, -200, -100),
    (-150, -50, -200, -90),
    (-220, -140, -50, 50),
]

hex_vertices = []   # will store (x, y)

for (xmin, xmax, ymin, ymax) in windows:
    value, x, y = max_in_window(Z, xmin, xmax, ymin, ymax)
    hex_vertices.append((x, y))
    print(f"Window {xmin},{xmax},{ymin},{ymax} -> max at (x={x}, y={y}) value={value}")

def polygon_angles(points):
    n = len(points)
    angles = []

    for i in range(n):
        p_prev = np.array(points[i-1])
        p_curr = np.array(points[i])
        p_next = np.array(points[(i+1) % n])

        v1 = p_prev - p_curr
        v2 = p_next - p_curr

        dot = np.dot(v1, v2)
        n1 = np.linalg.norm(v1)
        n2 = np.linalg.norm(v2)

        cos_theta = dot / (n1 * n2)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

        angle = np.degrees(np.arccos(cos_theta))
        angles.append(angle)

    return angles

# Compute angles
angles = polygon_angles(hex_vertices)

# Print angles and norms from (0,0)
for i, (vertex, ang) in enumerate(zip(hex_vertices, angles)):
    norm = np.linalg.norm(vertex)
    print(f"Vertex {i+1}: x={vertex[0]}, y={vertex[1]}, angle={ang:.3f} deg, ||vertex||={norm:.3f}")

# Convert Z into NumPy array
xs = range(-600, 601)
ys = range(-600, 601)

data = np.zeros((len(xs), len(ys)))
for i, x in enumerate(xs):
    for j, y in enumerate(ys):
        data[i, j] = Z[x][y]

# --- Create circular mask ---
N = len(xs)
X, Y = np.meshgrid(xs, ys, indexing='ij')
R = np.sqrt(X**2 + Y**2)
radius = max(xs)           # 600
mask = R > radius          # outside circle

# Apply mask using alpha channel
# alpha=1 inside circle, alpha=0 outside
alpha = (~mask).astype(float)

plt.figure(figsize=(12, 12), dpi=600)

plt.imshow(
    data.T,
    origin='lower',
    extent=[min(xs), max(xs), min(ys), max(ys)],
    cmap='inferno',
    norm=LogNorm(),
    interpolation='bicubic',
    alpha=alpha.T           # mask applied here
)

# Remove axes, ticks, frame
plt.axis('off')

# Ensure square aspect ratio
plt.gca().set_aspect('equal')

plt.tight_layout(pad=0)

# --- Save PNG ---
plt.savefig(
    os.path.join(current_dir, f"{file_name}_hires.png"),
    dpi=600,
    bbox_inches='tight',
    pad_inches=0,
    transparent=True
)

# --- Save SVG ---
plt.savefig(
    os.path.join(current_dir, f"{file_name}.svg"),
    bbox_inches='tight',
    pad_inches=0,
    transparent=True
)

plt.show()