# Raytracer

A toy raytracer for learning.

## Features
- A monte-carlo raytracer
- Diffuse, Specular, Mirror materials
- Sphere, plane geometries
- Multithreading via openMP

## TODO
### Misc. Optimization
- Test Vec4
- Investigate using plane reflection to generate rotated hemisphere over rotation matrices. Branching can be avoided by being clever
- Precomputing rotation matrix since the second Vec is a constant (0,0,1)

### Raycasting improvements
- Implement AABB
- BVH acceleration (see PBRT book)
- Optimization for shadow rays (early termination by stopping on any found intersection between 2 points)
- Improve self intersection handling (see raytracing gems)

### Variance reduction
- Russian roulette
- Metropolis-Hastings
- Multiple importance sampling

### CUDA
- Use small kernels over a singular megakernel (see jacco blog)
- Data oriented design (DOD) such that all materials and geometric primitives are be processed in their own array
- CUDA SIMD intrinsics for even more speed

### Other
- Spectral rendering
- Movable prerender cameraview to get good angles of the scene