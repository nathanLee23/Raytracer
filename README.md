# Raytracer

A WIP toy raytracer for learning.

## Samples
Diffuse material:  
<img src="samples/sample1.png" height="600">

Glass material:  
<img src="samples/sample2.png" height="600">

Multiple importance sampling in a slightly more complicated scene:  
<img src="samples/sample3.png" height="600">

## Features
- A monte-carlo pathtracer
- Diffuse, Specular, Mirror materials
- OBJ loading
- GUI to watch in progress render
- Multithreaded rendering using openMP

## TODO
### Features
- ~~Import .obj files~~
- Import cubemaps
- Interactive and movable camera
- Spectral rendering
- Scene editing

### Bugs
- Specular refractive index fix

### Raycasting improvements
- ~~Implement AABB~~
- BVH acceleration (see PBRT book)
- Optimization for shadow rays (early termination by stopping on any found intersection between 2 points)

### Variance reduction
- ~~Russian roulette~~
  - ~~Brightness based russian roulette~~
- ~~Hemisphere sampling~~
- ~~Multiple importance sampling~~
- Low discrepancy sequences

### CUDA
- Use small kernels over a singular megakernel (see jacco blog)
- Data oriented design (DOD) such that all the attributes of a certain types are arrays

### Misc. Optimization
- Test Vec4

## Dependancies
You'll have to install Intel Embree to be able to compile and run this at the moment
