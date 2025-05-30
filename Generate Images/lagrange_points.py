import matplotlib.pyplot as plt
import numpy as np

# Define positions in DU (non-dimensional units)
earth_pos = (0, 0)
moon_pos = (1, 0)
L1 = (0.8369, 0)
L2 = (1.1557, 0)
L3 = (-1.005, 0)
L4 = (0.4878, 0.8660)
L5 = (0.4878, -0.8660)

# Triangle connecting L4-L1-L5
triangle_x = [L4[0], L1[0], L5[0], L4[0]]
triangle_y = [L4[1], L1[1], L5[1], L4[1]]

# Plot setup
plt.figure(figsize=(6, 6))
plt.plot(triangle_x, triangle_y, 'b-', label="Lagrange Triangle")
plt.plot([earth_pos[0], moon_pos[0]], [earth_pos[1], moon_pos[1]], 'b-')  # Earth-Moon line

# Plot Earth and Moon
plt.plot(*earth_pos, 'go', markersize=10, label='Earth')
plt.plot(*moon_pos, 'ko', markersize=5, label='Moon')
plt.text(earth_pos[0] - 0.1, earth_pos[1] + 0.05, 'Earth', fontsize=10, color='green')
plt.text(moon_pos[0] + 0.05, moon_pos[1] + 0.05, 'Moon', fontsize=10, color='black')

# Plot Lagrange points
lagrange_points = {'L1': L1, 'L2': L2, 'L3': L3, 'L4': L4, 'L5': L5}
for label, pos in lagrange_points.items():
    plt.plot(*pos, 'r+', markersize=12, markeredgewidth=2)
    plt.text(pos[0] + 0.03, pos[1] + 0.03, label, fontsize=10, color='red')

# Add angle annotation
theta_text_pos = (0.25, 0.1)
plt.text(*theta_text_pos, r'$\theta = 60^\circ$', fontsize=12)

# Axis settings
plt.xlabel('$x (l^*)$')
plt.ylabel('$y (l^*)$')
plt.axis('equal')
plt.grid(True)
plt.tight_layout()

# Show plot
plt.show()
