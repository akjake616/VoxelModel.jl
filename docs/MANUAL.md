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
    start::Vector{<:Real} = [shift_half * 1 / 2, shift_half * 1 / 2, shift_half * 1 / 2]
end
 ```

## APIs

- ### Global Settings
___
```julia
reset_voxel()
```

Reset the full voxel space. 

___
```julia
reset_ref(b::Bool)
```

Reset whether to show the refenrnce axes at the origin. The default is on. 

___
```julia
reset_shift(b::Bool)
```

Set `shift_half` to `true` or `false`. `shift_half` denotes whether to shift the grid origin half apart from the grid spacing. For example, if the grid spacing is 1, the orign of the grids are at the positions of 0.5, 1.5, 2.5... etc. It is set to `true` by default. Noted that the call of this function will reset the grid space.

___
```julia
reset_dl(dl::Vector{<: Real})
```

Set the grid spacing to `dl` (a 3 element vector). The initial grid spacing is `[1, 1, 1]`. Noted that the call of this function will reset tWhe grid space.
___
- ### Core Functionalities
___
```julia
export_space()
```

Export the the voxel model.

___
```julia
export_grid()
```

Export the grid array filled with color indexes. Note that when geometries coincide, index of the last-added geometry is taken. 

___
```julia
create_cube(origin::Vector{<:Real}, dimension::Vector{<:Real}, ind::Int=1, mode="corner", fac::Real=2)
```

 Create a cuboid with origin at `origin` and size `dimension` (color index `ind`). `mode="corner"` is the default option; otherwise, one can choose `mode="center"`, which specifies the origin as the center of the cuboid.  `fac` denotes the interior densified factor according to the grid spacing. It is suggested that one take `fac=2` if rotation is needed, otherwise choose `fac=1`.

___
```julia
create_sphere(origin::Vector{<:Real}, radius::Real, ind::Int=1, fac::Real=2)
```

Create a sphere with origin at `origin` and a radius equals to `radius` (color index `ind`).

___
```julia
create_ellip(origin::Vector{<: Real}, dimension::Vector{<: Real}, ind::Int=1, fac::Real=2)
```

Create an ellipsoid with origin at `origin` and the length of the semi-axes `dimension` (color index `ind`).

___
```julia
create_cylin(origin::Vector{<:Real}, radius::Real, height::Real, ind::Int=1, fac::Real=2)
```

Create a cylinder with `radius` and `height` from the base `origin` (color index `ind`).

___
```julia
translate!(geo::Geometry, dl::Vector{<: Real})
```

Translate geometry with `dl`.

___
```julia
rotate!(geo::Geometry, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])
```

Rotate geometry with angle = `ang` with respect to the axis `axis` according to the origin `origin`. If origin is not specified, then the rotation is conducted according to the center of the geometry.

___
```julia
clear_geom(geo::Geometry)
clear_geom(geoList::Vector{Geometry})
```

Clear geometry from the voxel space.

___
 ```julia
plot_voxel(addRef::Bool=true)
```

Plot the voxel space. If `addRef=false` the reference axes will not be added. Call this function if the plot window is closed accidentally. 
___