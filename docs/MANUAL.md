# Manual for VoxelModel.jl

## Overview

The `Geometry` struct represents a geometric shape with a list of positions `pos`, an integer index `index` and a unique `ID`:

```julia
mutable struct Geometry
    pos::Vector{Vector{Float64}}
    index::Int
    const ID::Int
end
```

The `index` value is analogous to material index in FDTD/PSTD simulations. In `VoxelModel.jl`, they are used to set the colors of the geometries. In order to set the color, one can edit the dict `colorDict`, for example:

```julia
colorDict[4] = "pink"
 ```

sets the grids with index value `4` as pink voxels. If the key is not found, a random color is assigned. It is noted that key `0` is used for geometry deletion (please refer to  `ex_advanced.jl` in the examples folder).

Normally one does not need to deal with the `Geometry` struct directly. To add geometries, use API functions which starts with `craete_` to build geometries. `trans!` and `rot!` are used for translations and rotations of the created geometry.

The `Voxels` struct represents the voxel space as a result of the prsented geometries. The struct is composed of `grid`, which stores the grid array with the integer index of the geometries (last added), with customizable grid spacing `dl` and start position `start` (default of `shift[]` is `true`):

 ```julia
Base.@kwdef mutable struct Voxels
    grid::Array{Int, 3} = zeros(1, 1, 1)
    dl::Vector{Float64} = [1.0, 1.0, 1.0]
    start::Vector{Float64} = [shift[] * 0.5, shift[] * 0.5, shift[] * 0.5]
end
 ```

If one needs to add extra traces on the voxel plot, the PlotlyJS Plot is exported as `canvas` which can be used for further modifications. 

## API

### reset voxel space

```julia
reset_voxel()
```

Reset the full voxel space.

#### Arguments
- None

___

### export voxel

```julia
export_voxel()
```

Exports the current voxel as a `Voxels` struct.


#### Returns
- `voxel::Voxels`: The current voxel space.

___


### save voxel

```julia
save_voxel(fileName::String)
```
save voxel in JLD format. 

#### Arguments
- `fileName::String`: file name of the JLD file

___

### load voxel

```julia
load_voxel(fileName::String)
```

load voxel in JLD format. This will reset the current voxel space.


#### Arguments
- `fileName::String`: file name of the JLD file

___

### plot voxel

```julia
plot_voxel(addRef::Bool=true)
```

Plots the voxel space. If `addRef=false`, the reference axes will not be added. Call this function if the plot window is accidentally closed.

#### Arguments
- `addRef::Bool=true`: Boolean value to specify whether to add the reference axes to the plot.

___

### reset referecnce axes

```julia
reset_ref(b::Bool)
```

Toggles the display of the reference axes at the origin. The default state is `true` (axes visible).

#### Arguments
- `b::Bool`: Boolean value to set the visibility of the reference axes.

___

### reset shift 

```julia
reset_shift(b::Bool)
```

Sets the `shift[]` parameter to the specified boolean value `b`. `shift[]` means whether to shift the center with half grid spacing. For example, if the spacing is `1`, the grid center will be at `-0.5, 0.5, 1.5, ...`. The default of `shift[]` is `true`. 

#### Arguments
- `b::Bool`: Boolean value to set the `shift[]` parameter.

___

### reset_dl

```julia
reset_dl(dl::Vector{<:Real})
```

Updates the grid spacing to the specified vector `dl`.

#### Arguments
- `dl::Vector{<:Real}`: A vector containing the new grid spacing values.

___

### create_cuboid

```julia
create_cuboid(origin::Vector{<:Real}, dim::Vector{<:Real}, ind::Int=1, mode::String="corner", fac::Real=2)
```

Creates a cuboid with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The origin point of the cuboid.
- `dim::Vector{<:Real}`: The dimensions of the cuboid.
- `ind::Int=1`: The color index of the cuboid.
- `mode::String="corner"`: The mode specifying the cuboid's origin ("corner" or "center").
- `fac::Real=2`: The interior densified factor according to the grid spacing.

___

### create_cube

```julia
creat_cube(origin::Vector{<:Real}, dim::Real, ind::Int=1, mode="corner", fac::Real=2)
```

Creates a cube with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The origin point of the cube.
- `dim::Real`: The side length of the cube.
- `ind::Int=1`: The color index of the cube.
- `mode::String="corner"`: The mode specifying the cube's origin ("corner" or "center").
- `fac::Real=2`: The interior densified factor according to the grid spacing.

___

### create_sphere

```julia
create_sphere(origin::Vector{<:Real}, radius::Real, ind::Int=1, fac::Real=2)
```

Creates a sphere with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The origin point of the sphere.
- `radius::Real`: The radius of the sphere.
- `ind::Int=1`: The color index of the sphere.
- `fac::Real=2`: The interior densified factor according to the grid spacing.

___

### create_ellipsoid

```julia
create_ellipsoid(origin::Vector{<:Real}, par::Vector{<:Real}, ind::Int=1, fac::Real=2)
```

Creates an ellipsoid with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The origin point of the ellipsoid.
- `par::Vector{<:Real}`: The lengths of the semi-axes.
- `ind::Int=1`: The color index of the ellipsoid.
- `fac::Real=2`: The interior densified factor according to the grid spacing.

___

### create_cylinder

```julia
create_cylinder(origin::Vector{<:Real}, radius::Real, height::Real, ind::Int=1, fac::Real=2)
```

Creates a cylinder with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The base origin point of the cylinder.
- `radius::Real`: The radius of the cylinder.
- `height::Real`: The height of the cylinder.
- `ind::Int=1`: The color index of the cylinder.
- `fac::Real=2`: The interior densified factor according to the grid spacing.

___

### trans!

```julia
trans!(geo::Geometry, dl::Vector{<:Real})
```

Translates the geometry by the specified vector `dl`.

#### Arguments
- `geo::Geometry`: The geometry to be translated.
- `dl::Vector{<:Real}`: The translation vector.

___

### rot!

```julia
rot!(geo::Geometry, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])
```

Rotates the geometry by the specified angle `ang` around the axis `axis` and origin `origin`.

#### Arguments
- `geo::Geometry`: The geometry to be rotated.
- `ang::Real`: The rotation angle.
- `axis::Vector{<:Real}`: The rotation axis.
- `origin::Vector{<:Real}=[0]`: The rotation origin. Defaults to the center of the geometry if not specified.

___

### clear_geom (single geometry)

```julia
clear_geom(geo::Geometry)
```

Removes the specified geometry from the voxel space.

#### Arguments
- `geo::Geometry`: The geometry to be removed.

___

### clear_geom (multiple geometries)

```julia
clear_geom(geoList::Vector{Geometry})
```

Removes the specified list of geometries from the voxel space.

#### Arguments
- `geoList::Vector{Geometry}`: The list of geometries to be removed.

___


### export_grid

```julia
export_grid()
```

Exports the grid array filled with color indexes. Note that when geometries overlap, the index of the last-added geometry is used.

#### Returns
- `Array{Int}`: The grid array with color indexes.

___