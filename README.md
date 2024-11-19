# Gravitational Simulation of Colliding Galaxies

## Overview
This project simulates the gravitational interactions and eventual collision of galaxies using a **Barnes-Hut algorithm** for efficient computation. It models up to **2,000 particles**(until parallization fully implemented), representing stars or any other stellar objects, with realistic physics to study the dynamics of galactic evolution and interaction.
---

## Features
- **Barnes-Hut Algorithm**: Efficient gravitational force calculations using an octree to approximate far-field forces.
- **Dynamic Galaxy Collisions**: Simulates two or more galaxies interacting and merging due to gravitational forces.
- **Realistic Physics**:
  - Gravity based on Newtonian mechanics.
  - Accurate computation of center of mass and total mass for tree nodes.
- **Scalability**: Designed to handle up to 10,000 objects.
- **3D Visualization**: Real-time galaxy interactions with adjustable camera controls.
- **Customizable Parameters**:
  - Initial positions, velocities, and masses of celestial objects.
  - Adjustable time step for precision and performance.
- **Real Time Updating** Octree Visualization

---

## Requirements

### Software
- Unity 2021.3 or higher.
- .NET Framework 4.7 or higher (for C# scripting).

### Hardware
- Minimum 8GB RAM and a dedicated GPU recommended for smooth visualization.

---

## How It Works

### Initialization
- Objects are initialized with positions, masses, and velocities.
- The octree is constructed based on object positions.

### Force Calculation
- The Barnes-Hut algorithm divides space into a hierarchical octree structure.
- Nodes far away are treated as single points to approximate gravitational forces.

### Update Loop
- Forces are calculated for each object.
- Positions and velocities are updated using Verlet integration.
---
## Current Features in Progress
- Multithreading as a benchmark before fully working on the GPU
---
##Future Additions
- Fully move to GPU for all computations
- Instead of directly visualizing the particles, return txt files containing coordinates that I can make a new scene to loop through these files for much more computationally complex visualizations eg. Raymarching.
