#region API
"""
    reset_voxel()

    Reset the full voxel space. 
"""
function reset_voxel()
    global space = Voxels()
end

"""
    reset_ref(b::Bool)

    Toggles the display of the reference axes at the origin. The default state is `true` (axes visible).
    
    # Arguments
    - `b::Bool`: Boolean value to set the visibility of the reference axes.
"""
function reset_ref(b::Bool)
    refAxis[] = b
    plot_voxel(refAxis[])
end

"""
    reset_shift(b::Bool)

    Sets the `shift[]` parameter to the specified boolean value `b`.
    
    # Arguments
    - `b::Bool`: Boolean value to set the `shift[]` parameter.
"""
function reset_shift(b::Bool)
    shift[] = b
    global space = Voxels()
end

"""
    reset_dl(dl::Vector{<:Real})

    Updates the grid spacing to the specified vector `dl`.
    
    # Arguments
    - `dl::Vector{<:Real}`: A vector containing the new grid spacing values.
"""
function reset_dl(dl::Vector{<:Real})
    @assert length(dl) == 3
    start = [shift[] * 1 / 2 * dl[1], shift[] * 1 / 2 * dl[2], shift[] * 1 / 2 * dl[3]]
    global space = Voxels([], dl, start)
end

"""
    create_cuboid(origin::Vector{<:Real}, dim::Vector{<:Real}, ind::Int=1, mode::String="corner", fac::Real=2)

    Creates a cuboid with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The origin point of the cuboid.
    - `dim::Vector{<:Real}`: The dimensions of the cuboid.
    - `ind::Int=1`: The color index of the cuboid.
    - `mode::String="corner"`: The mode specifying the cuboid's origin ("corner" or "center").
    - `fac::Real=2`: The interior densified factor according to the grid spacing.
"""
function create_cuboid(origin::Vector{<:Real}, dim::Vector{<:Real}, ind::Int=1, mode="corner", fac::Real=2)
    @assert length(origin) == 3
    @assert length(dim) == 3
    @assert fac > 0

    dx = space.dl[1]
    dy = space.dl[2]
    dz = space.dl[3]

    sx = Int(_round((dim[1]-2*dx/fac) / (dx/fac))) + 1
    sy = Int(_round((dim[2]-2*dy/fac) / (dy/fac))) + 1
    sz = Int(_round((dim[3]-2*dz/fac) / (dz/fac))) + 1

    pos = []

    if mode == "center"
        xs = origin[1] - dim[1]/2 + dx / fac
        ys = origin[2] - dim[2]/2 + dy / fac
        zs = origin[3] - dim[3]/2 + dz / fac
    else
        xs = origin[1] + dx / fac
        ys = origin[2] + dy / fac
        zs = origin[3] + dz / fac
    end

    for i in 1:sx, j in 1:sy, k in 1:sz
        push!(pos, [xs + (i - 1) * dx/fac, ys + (j - 1) * dy/fac, zs + (k - 1) * dz/fac])
    end

    idCount[] += 1
    geo = Geometry(pos, ind, idCount[])
    idDict[][idCount[]] = ind
    _add_geom(geo)

    _plot_voxel(refAxis[])

    return geo
end

"""
    create_sphere(origin::Vector{<:Real}, radius::Real, ind::Int=1, fac::Real=2)

    Creates a sphere with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The origin point of the sphere.
    - `radius::Real`: The radius of the sphere.
    - `ind::Int=1`: The color index of the sphere.
    - `fac::Real=2`: The interior densified factor according to the grid spacing.
"""
function create_sphere(origin::Vector{<:Real}, radius::Real, ind::Int=1, fac::Real=2)
    @assert length(origin) == 3
    @assert fac > 0

    dx = space.dl[1]
    dy = space.dl[2]
    dz = space.dl[3]

    sr = maximum([ceil(Int, radius / (dx/fac)), ceil(Int, radius / (dy/fac)), ceil(Int, radius / (dz/fac))])

    pos = []
    for i in -sr:sr, j in -sr:sr, k in -sr:sr

        x1 = origin[1] + i * (dx/fac) 
        y1 = origin[2] + j * (dy/fac) 
        z1 = origin[3] + k * (dz/fac) 

        if ((x1 - origin[1])./(radius-dx/fac))^2 + ((y1 - origin[2])./(radius-dy/fac))^2 + ((z1 - origin[3])./(radius-dz/fac))^2 < 1
            push!(pos, [x1, y1, z1])
        end
    end

    idCount[] += 1
    geo = Geometry(pos, ind, idCount[])
    idDict[][idCount[]] = ind
    _add_geom(geo)

    _plot_voxel(refAxis[])

    return geo
end

"""
    create_ellip(origin::Vector{<:Real}, par::Vector{<:Real}, ind::Int=1, fac::Real=2)

    Creates an ellipsoid with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The origin point of the ellipsoid.
    - `par::Vector{<:Real}`: The lengths of the semi-axes.
    - `ind::Int=1`: The color index of the ellipsoid.
    - `fac::Real=2`: The interior densified factor according to the grid spacing.
"""
function create_ellip(origin::Vector{<:Real}, par::Vector{<:Real}, ind::Int=1, fac::Real=2)
    @assert length(origin) == 3
    @assert length(par) == 3
    @assert fac > 0

    dx = space.dl[1]
    dy = space.dl[2]
    dz = space.dl[3]

    sa = ceil(Int, par[1] / (dx/fac))
    sb = ceil(Int, par[2] / (dy/fac))
    sc = ceil(Int, par[3] / (dz/fac))

    pos = []
    for i in -sa:sa, j in -sb:sb, k in -sc:sc
        x1 = origin[1] + i * (dx/fac) 
        y1 = origin[2] + j * (dy/fac) 
        z1 = origin[3] + k * (dz/fac) 

        if ((x1 - origin[1]) / (par[1]-dx/fac))^2 + ((y1 - origin[2]) / (par[2]-dy/fac))^2 + ((z1 - origin[3]) / (par[3]-dz/fac))^2 < 1
            push!(pos, [x1, y1, z1])
        end
    end

    idCount[] += 1
    geo = Geometry(pos, ind, idCount[])
    idDict[][idCount[]] = ind
    _add_geom(geo)

    _plot_voxel(refAxis[])

    return geo
end

"""
    create_cylin(origin::Vector{<:Real}, radius::Real, height::Real, ind::Int=1, fac::Real=2)

    Creates a cylinder with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The base origin point of the cylinder.
    - `radius::Real`: The radius of the cylinder.
    - `height::Real`: The height of the cylinder.
    - `ind::Int=1`: The color index of the cylinder.
    - `fac::Real=2`: The interior densified factor according to the grid spacing.
"""
function create_cylin(origin::Vector{<:Real}, radius::Real, height::Real, ind::Int=1, fac::Real=2)
    @assert length(origin) == 3
    @assert fac > 0

    dx = space.dl[1]
    dy = space.dl[2]
    dz = space.dl[3]

    sr = maximum([ceil(Int, radius / (dx/fac)), ceil(Int, radius / (dy/fac))])
    sz = Int(_round((height-2*dz/fac) / (dz/fac))) + 1
    pos = []
    for i in -sr:sr, j in -sr:sr, k in 1:sz

        x1 = origin[1] + i * (dx/fac) 
        y1 = origin[2] + j * (dy/fac) 
        z1 = origin[3] + k * (dz/fac) 

        if ((x1 - origin[1])./(radius-dx/fac))^2 + ((y1 - origin[2])./(radius-dy/fac))^2  < 1
            push!(pos, [x1, y1, z1])
        end
    end

    idCount[] += 1
    geo = Geometry(pos, ind, idCount[])
    idDict[][idCount[]] = ind
    _add_geom(geo)

    _plot_voxel(refAxis[])

    return geo
end

"""
    trans!(geo::Geometry, dl::Vector{<:Real})

    Translates the geometry by the specified vector `dl`.
    
    # Arguments
    - `geo::Geometry`: The geometry to be translated.
    - `dl::Vector{<:Real}`: The translation vector.
"""
function trans!(geo::Geometry, dl::Vector{<:Real})
    @assert length(dl) == 3

    _del_geom(geo)

    for n in eachindex(geo.pos)
        geo.pos[n][1] += dl[1]
        geo.pos[n][2] += dl[2]
        geo.pos[n][1] += dl[3]
    end

    _add_geom(geo)

    _plot_voxel(refAxis[])
end

"""
    rot!(geo::Geometry, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])

    Rotates the geometry by the specified angle `ang` around the axis `axis` and origin `origin`.
    
    # Arguments
    - `geo::Geometry`: The geometry to be rotated.
    - `ang::Real`: The rotation angle.
    - `axis::Vector{<:Real}`: The rotation axis.
    - `origin::Vector{<:Real}=[0]`: The rotation origin. Defaults to the center of the geometry if not specified.
"""
function rot!(geo::Geometry, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])
    @assert length(axis) == 3

    _del_geom(geo)

    axis = axis ./ norm(axis)
    vrot = similar(geo.pos)
    
    if origin == [0] # rotation center set at the geometry center
        origin = sum(geo.pos) ./ length(geo.pos)
    else
        @assert length(origin) == 3
    end

    for n in eachindex(vrot)
        v = (geo.pos[n] .- origin)
        vrot[n] = cosd(ang) * v + sind(ang) * cross(axis, v) + (1-cosd(ang)) * dot(axis, v) * axis
        geo.pos[n] = vrot[n] .+ origin
    end

    _add_geom(geo)

    _plot_voxel(refAxis[])
end

"""
    clear_geom(geo::Geometry)

    Removes the specified geometry from the voxel space.
    
    # Arguments
    - `geo::Geometry`: The geometry to be removed.
"""
function clear_geom(geo::Geometry)

    _del_geom(geo)
    geo = nothing
    _plot_voxel(refAxis[])
end

"""
    clear_geom(geoList::Vector{Geometry})

    Removes the specified list of geometries from the voxel space.
    
    # Arguments
    - `geoList::Vector{Geometry}`: The list of geometries to be removed.
"""
function clear_geom(geoList::Vector{Geometry})

    for i in eachindex(geoList)
        _del_geom(geoList[i])
        geoList[i] = nothing
    end
end

"""
    plot_voxel(addRef::Bool=true)

    Plots the voxel space. If `addRef=false`, the reference axes will not be added. Call this function if the plot window is accidentally closed.
    
    # Arguments
    - `addRef::Bool=true`: Boolean value to specify whether to add the reference axes to the plot.
"""
function plot_voxel(addRef::Bool=true)
    
    global canvas = plot([mesh3d(x=0, y=0, z=0)], blank_layout())
    display(canvas)
    
    _plot_voxel(addRef)
end

"""
    export_voxel()

    Exports the current voxel space as a `Voxels` struct.
    
    # Returns
    - `Voxels`: The current voxel space.
"""
function export_voxel()
    return space
end

"""
    export_grid()

    Exports the grid array filled with color indexes. Note that when geometries overlap, the index of the last-added geometry is used.
    
    # Returns
    - `Array{Int}`: The grid array with color indexes.
"""
function export_grid()
    grid = zeros(Int, size(space.gridID))
    for i in eachindex(grid)
        if space.gridID[i] == []
            grid[i] = 0
        else
            grid[i] = idDict[][space.gridID[i][end]]
        end
    end
    return grid
end
#endregion

