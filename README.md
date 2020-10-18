# Raytracer

A toy raytracer for learning.

## Features
- A monte-carlo raytracer
- Diffuse, Specular, Mirror materials
- Sphere, plane, box geometries
- Multithreading via openMP

## TODO
### Misc. Features
- BDPT
- Import .obj files
- Skybox
- Movable camera preview to get good camera angles and positioning before rendering
- Spectral rendering
- alias tables for O(1) sampling

### Bugs
- Specular refractive index fix

### Raycasting improvements
- ~~Implement AABB~~
- BVH acceleration (see PBRT book)
- Optimization for shadow rays (early termination by stopping on any found intersection between 2 points)
- Improve self intersection handling (see raytracing gems)

### Variance reduction
- ~~Russian roulette~~
- ~~Hemisphere sampling~~
- Multiple importance sampling
- Metropolis-Hastings (Metropolis Light Transport or Energy redistribution path tracing)

### CUDA
- Use small kernels over a singular megakernel (see jacco blog)
- Data oriented design (DOD) such that all the attributes of a certain types are arrays
- CUDA SIMD intrinsics for even more speed

### Misc. Optimization
- Test Vec4
