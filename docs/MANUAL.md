# Manual for VoxelModel.jl

## Overview

It is best for the user to understand the usage by the examples provided in the exa,ples folder. To undersatnd the underlying design, we first introduce the `Geometry` struct. `Geometry` represents a geometric shape with a list of positions `pos`, an integer index `index` and a unique `ID`:

```julia
mutable struct Geometry
    pos::Vector{Vector{Float64}}
    index::Int
    const ID::Int
end
```
The positions `pos` are actually a set of equally spaced dense grid points which is used for further voxel approximations. The `index` value is analogous to material index in FDTD/PSTD simulations. In `VoxelModel.jl`, they are used to set the colors of the geometries. In order to set the color, one can edit the dict `colorDict`, for example:

```julia
colorDict[4] = "pink"
 ```

sets the grids with index value `4` as pink voxels. If the key is not found, a random color is assigned. It is noted that key `0` is used for geometry deletion (please refer to  `ex_advanced.jl` in the examples folder). `ID` is a unique ID number (used for efficient identification of the geometries in the voxel space.

Normally one does not need to deal with the `Geometry` struct directly. To add geometries, use API functions which starts with `craete_` to build geometries. `trans!` and `rot!` are used for translations and rotations of the created geometry. It is noted that there is a keyword `render` (default`=false`) determining whether to plot the voxel space according to the creation/operation of the geometries. Due to the time to first plot (TTFP) problem, it is recommended to complete all creations/operations first and then call `plot_voxel()` to see the reult.

The `Voxels` struct represents the voxel space as a result of the present geometries. The struct is composed of `grid`, which stores the grid array with the integer index of the geometries (only the last added geometry is stored in `grid` if geometry collision occurs), with customizable grid spacing `dl` and start position `start` (default of `shift[]` is `true`):

 ```julia
Base.@kwdef mutable struct Voxels
    grid::Array{Int, 3} = zeros(1, 1, 1)
    dl::Vector{Float64} = [1.0, 1.0, 1.0]
    start::Vector{Float64} = [shift[] * 0.5, shift[] * 0.5, shift[] * 0.5]
end
 ```

`shift[]` is a boolean value representing whether to shift the center of each voxel with half grid spacing (one can set the value of `shift[]` using th API function `reset_shift()`). For example, if the spacing is `1`, the voxel center will be at `-0.5, 0.5, 1.5, ...` if `shift[] = true`, and the corner of each voxel will be at  `-1, 0, 1,...` etc. It is noted that when crating a gemetries, one should be cautious about grid shifts. For the default setting (`shift[] = true` and `dl = [1, 1, 1]`), consider the following example:

```julia
c1 = create_cube([0, 0, 0], 1, "corner") # OK - creates a voxl centered at [0.5, 0.5, 0.5] with side length = 1
c2 = create_cube([0, 0, 0], 1, "center") # !! - creates a voxl centered at [0.0, 0.0, 0.0] with side length = 2
```
`c1` is created as expected, but `c2` is not. This is because, to preserve the center position of the cube (which does not lie on the grid centers), the size is scaled (rounded) to 2 in this case.

If one needs to add extra traces on the voxel plot, the PlotlyJS Plot is exported as `canvas` which can be used for further modifications. 

## API

### reset voxel space

```julia
reset_voxel(;render=false)
```

Reset the full voxel space.

#### Arguments
- None

#### Keywords
- `render=false`: real-time rendering for creation/operation.
___

### export voxel

```julia
export_voxel()
```

Exports a copy of the current voxel.


#### Returns
- `voxelCopy::Voxels`: The copy of the current voxel space.

___


### save voxel

```julia
save_voxel(fileName::String)
```
save voxel in JLD format. 

#### Arguments
- `fileName::String`: file name of the JLD file.

___

### load voxel

```julia
load_voxel(fileName::String)
```

load voxel in JLD format. This will reset the current voxel space.


#### Arguments
- `fileName::String`: file name of the JLD file.

___

### plot voxel

```julia
plot_voxel(addRef::Bool=true)
```

Plots the voxel space. If `addRef=false`, the reference axes will not be added. Call this function if the plot window is accidentally closed.

#### Arguments
- `addRef::Bool=true`: Boolean value to specify whether to add the reference axes to the plot.

___

### assign voxel

```julia
assign_voxel(grid::Array{Int, 3}, dl::Vector{<:Real}=[1.0, 1.0, 1.0], start::Vector{<:Real}=[0, 0, 0])
```


Assign grid to voxel space.
Updates: One can also set `dl` as a real number for equal spacings.
    
#### Arguments
- `grid::Array{Int, 3}`: Interger grid array.
- `dl::Vector{<:Real}=[1.0, 1.0, 1.0]`: grid spacings.
- `start::Vector{<:Real}=[shift[] * 0.5, shift[] * 0.5, shift[] * 0.5]`: start point.

___

### reset referecnce axes

```julia
reset_ref(b::Bool, len::Real=refLen[])
```

Toggles the display of the reference axes at the origin. The default state is `true` (axes visible).Since the gird space is shifted, the voxel space will be reseted accordingly. Use this function before adding geomtries to the voxel space.

#### Arguments
- `b::Bool`: Boolean value to set the visibility of the reference axes.
- `len::Float64=refLen[]`: reference length of the axes. The default is the minimum of the grid spacings.

___

### reset shift 

```julia
reset_shift(b::Bool)
```

Sets the `shift[]` parameter to the specified boolean value `b`. 
!! Since the gird space is changed, the voxel space will be reseted accordingly. Use this function before adding geomtries to the voxel space.

#### Arguments
- `b::Bool`: Boolean value to set the `shift[]` parameter.

___

### reset grid spacings

```julia
reset_dl(dl::Vector{<:Real})
```

Resets the grid spacing to the specified vector `dl`. 
!! Since the gird space is changed, the voxel space will be reseted accordingly. Use this function before adding geomtries to the voxel space.
Updates: One can also set `dl` as a real number for equal spacings.


#### Arguments
- `dl::Vector{<:Real}`: A vector containing the new grid spacing values.

___

### reset start

```julia
reset_start(start::Vector{<:Real})
```

reset the start point of the voxel space. 
!! Note that th spart point is ceiled to the neasrest grid center.
    
#### Arguments
- `start::Vector{<:Real}`: A vector indicating the spart point.

___

### create cuboid

```julia
create_cuboid(origin::Vector{<:Real}, dim::Vector{<:Real}, ind::Int=1, mode::String="corner", fac::Real=2; render=false)
```

Creates a cuboid with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The origin point of the cuboid.
- `dim::Vector{<:Real}`: The dimensions of the cuboid.
- `ind::Int=1`: The color index of the cuboid.
- `mode::String="corner"`: The mode specifying the cuboid's origin ("corner" or "center").
- `fac::Real=2`: The interior densified factor according to the grid spacing.

#### Keywords
- `render=false`: real-time rendering for creation/operation.

___

### create cube

```julia
creat_cube(origin::Vector{<:Real}, dim::Real, ind::Int=1, mode="corner", fac::Real=2; render=false)
```

Creates a cube with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The origin point of the cube.
- `dim::Real`: The side length of the cube.
- `ind::Int=1`: The color index of the cube.
- `mode::String="corner"`: The mode specifying the cube's origin ("corner" or "center").
- `fac::Real=2`: The interior densified factor according to the grid spacing.

#### Keywords
- `render=false`: real-time rendering for creation/operation.

___

### create sphere

```julia
create_sphere(origin::Vector{<:Real}, radius::Real, ind::Int=1, fac::Real=2; render=false)
```

Creates a sphere with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The origin point of the sphere.
- `radius::Real`: The radius of the sphere.
- `ind::Int=1`: The color index of the sphere.
- `fac::Real=2`: The interior densified factor according to the grid spacing.

#### Keywords
- `render=false`: real-time rendering for creation/operation.

___

### create ellipsoid

```julia
create_ellipsoid(origin::Vector{<:Real}, par::Vector{<:Real}, ind::Int=1, fac::Real=2; render=false)
```

Creates an ellipsoid with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The origin point of the ellipsoid.
- `par::Vector{<:Real}`: The lengths of the semi-axes.
- `ind::Int=1`: The color index of the ellipsoid.
- `fac::Real=2`: The interior densified factor according to the grid spacing.

#### Keywords
- `render=false`: real-time rendering for creation/operation.

___

### create cylinder

```julia
create_cylinder(origin::Vector{<:Real}, radius::Real, height::Real, ind::Int=1, fac::Real=2; render=false)
```

Creates a cylinder with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The base origin point of the cylinder.
- `radius::Real`: The radius of the cylinder.
- `height::Real`: The height of the cylinder.
- `ind::Int=1`: The color index of the cylinder.
- `fac::Real=2`: The interior densified factor according to the grid spacing.

#### Keywords
- `render=false`: real-time rendering for creation/operation.

___

### translation

```julia
trans!(geo::Geometry, dl::Vector{<:Real}; render=false)
```

Translates the geometry by the specified vector `dl`.

#### Arguments
- `geo::Geometry`: The geometry to be translated.
- `dl::Vector{<:Real}`: The translation vector.

#### Keywords
- `render=false`: real-time rendering for creation/operation.

___

### rotation

```julia
rot!(geo::Geometry, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0]; render=false)
```

Rotates the geometry by the specified angle `ang` around the axis `axis` and origin `origin`.

#### Arguments
- `geo::Geometry`: The geometry to be rotated.
- `ang::Real`: The rotation angle.
- `axis::Vector{<:Real}`: The rotation axis.
- `origin::Vector{<:Real}=[0]`: The rotation origin. Defaults to the center of the geometry if not specified.

#### Keywords
- `render=false`: real-time rendering for creation/operation.

___

### clear geometry (single geometry)

```julia
clear_geom(geo::Geometry; render=false)
```

Removes the specified geometry from the voxel space.

#### Arguments
- `geo::Geometry`: The geometry to be removed.

#### Keywords
- `render=false`: real-time rendering for creation/operation.

___

### clear geometry (multiple geometries)

```julia
clear_geom(geoList::Vector{Geometry}; render=false)
```

Removes the specified list of geometries from the voxel space.

#### Arguments
- `geoList::Vector{Geometry}`: The list of geometries to be removed.

#### Keywords
- `render=false`: real-time rendering for creation/operation.

___


### export grid

```julia
export_grid()
```

Exports the grid array filled with color indexes. Note that when geometries overlap, the index of the last-added geometry is used.

#### Returns
- `Array{Int}`: The grid array with color indexes.

___