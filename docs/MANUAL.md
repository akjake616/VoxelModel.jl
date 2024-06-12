# Manual for VoxelModel.jl

## Overview

The source folder is composed the main module `VoxelModel.jl`, the definition of color index file `color_dict.jl` and other supplementary files. In the modeule, the `Geometry` struct represents a geometric shape with a list of positions `pos`, a color index `index` and a unique `ID`:

```julia
mutable struct Geometry
    pos::Vector{Vector{<:Real}}
    index::Int
    const ID::Int
end
```

 The `Voxels` struct represents the voxel grid array `gridID` whcih stores the geometry ID, with customizable grid spacing `dl` and start position `start`:

 ```julia
Base.@kwdef mutable struct Voxels
    gridID::Array = []
    dl::Vector{<:Real} = [1, 1, 1]
    start::Vector{<:Real} = [shift_half[][] * 1 / 2, shift_half[][] * 1 / 2, shift_half[][] * 1 / 2]
end
 ```

 In order to set the color of the geometry, one can edit the dict `colorDict`. The following is the `color_dict.jl` file:

 ```julia
const colorDict = Dict{Int, String}()

colorDict[1] = "blue"
colorDict[2] = "green"
colorDict[3] = "red"
 ```

The default index 1, 2, 3 are mapped to color blue, green and red, respectively. One can add new color mapping (or modify the default), for example:

```julia
colorDict[4] = "pink"
 ```

If one needs to add extra traces on the voxel plot, the PlotlyJS Plot is exported as `canvas` which can be used for further modifications. 

## APIs

### reset_voxel

```julia
reset_voxel()
```

Reset the full voxel space.

#### Arguments
- None

___

### reset_ref

```julia
reset_ref(b::Bool)
```

Toggles the display of the reference axes at the origin. The default state is `true` (axes visible).

#### Arguments
- `b::Bool`: Boolean value to set the visibility of the reference axes.

___

### reset_shift

```julia
reset_shift(b::Bool)
```

Sets the `shift_half[]` parameter to the specified boolean value `b`. `shift_half[]` means whether to shift the center with half grid spacing. For example, if the spacing is `1`, the grid center will be at `-0.5, 0.5, 1.5, ...`. The default of `shift_half[]` is `true`. 

#### Arguments
- `b::Bool`: Boolean value to set the `shift_half[]` parameter.

___

### reset_dl

```julia
reset_dl(dl::Vector{<:Real})
```

Updates the grid spacing to the specified vector `dl`.

#### Arguments
- `dl::Vector{<:Real}`: A vector containing the new grid spacing values.

___

### create_cube

```julia
create_cube(origin::Vector{<:Real}, dim::Vector{<:Real}, ind::Int=1, mode::String="corner", fac::Real=2)
```

Creates a cuboid with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The origin point of the cuboid.
- `dim::Vector{<:Real}`: The dimensions of the cuboid.
- `ind::Int=1`: The color index of the cuboid.
- `mode::String="corner"`: The mode specifying the cuboid's origin ("corner" or "center").
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

### create_ellip

```julia
create_ellip(origin::Vector{<:Real}, par::Vector{<:Real}, ind::Int=1, fac::Real=2)
```

Creates an ellipsoid with the specified parameters.

#### Arguments
- `origin::Vector{<:Real}`: The origin point of the ellipsoid.
- `par::Vector{<:Real}`: The lengths of the semi-axes.
- `ind::Int=1`: The color index of the ellipsoid.
- `fac::Real=2`: The interior densified factor according to the grid spacing.

___

### create_cylin

```julia
create_cylin(origin::Vector{<:Real}, radius::Real, height::Real, ind::Int=1, fac::Real=2)
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

### plot_voxel

```julia
plot_voxel(addRef::Bool=true)
```

Plots the voxel space. If `addRef=false`, the reference axes will not be added. Call this function if the plot window is accidentally closed.

#### Arguments
- `addRef::Bool=true`: Boolean value to specify whether to add the reference axes to the plot.

___

### export_voxel

```julia
export_voxel()
```

Exports the current voxel space as a `Voxels` struct.

#### Returns
- `Voxels`: The current voxel space.

___

### export_grid

```julia
export_grid()
```

Exports the grid array filled with color indexes. Note that when geometries overlap, the index of the last-added geometry is used.

#### Returns
- `Array{Int}`: The grid array with color indexes.

___