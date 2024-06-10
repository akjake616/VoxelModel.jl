# VoxelModel.jl

`VoxelModel.jl` is a Julia module for creating, manipulating, and visualizing 3D voxel geometries. It provides functionalities to create geometric shapes, add or delete them from a voxel grid, and visualize the grid using [`PlotlyJS.jl`](https://github.com/JuliaPlots/PlotlyJS.jl). This project is dedicated to (hopefully) easy model creation in FDTD/PSTD simulations.

<p align="center">
  <img alt="VoxelModel.jl" src="./media/illus.png" width="50%" height="auto" />
</p>

## Installation

To install `VoxelModel.jl`, use the following command in the Julia REPL:

```julia
using Pkg
Pkg.add("VoxelModel")
```

## Learn by Examples

Please go to the examples folder for quick understanding of how to use the module:

```bash
cd examples
julia ex_basics.jl
```

This example code illustrates the basic functionalities of the module. To see how rotation works, please run the fllowing script in the examples folder:

```bash
julia ex_rotation.jl
```

To create more complicated geometries, please refer to the following example which demonstrates the creation of a sphere hollowd by cylinders:

```bash
julia ex_advanced.jl
```

## Usage

Please refer to the [user manual](./docs/MANUAL.md) in the docs folder.




