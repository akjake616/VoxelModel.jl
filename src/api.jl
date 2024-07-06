#region API
"""
    reset_voxel()

    Reset the full voxel voxel. 
"""
function reset_voxel()
    global voxel = Voxels()
    global gridID = []
    idCount[] = 0
    empty!(idDict)
    return nothing
end

"""
    export_voxel()

    Exports the current voxel as a `Voxels` struct.
    
    # Returns
    - `voxel::Voxels`: The current voxel space.
"""
function export_voxel()
    return voxel
end

"""
    save_voxel(fileName::String)

    save voxel in JLD format. 
    
    # Arguments
    - `fileName::String`: file name of the JLD file
"""
function save_voxel(fileName::String)
    save(fileName, "voxel", voxel)
    return nothing
end

"""
    load_voxel(fileName::String)

    load voxel in JLD format. This will reset the current voxel space.
    
    # Arguments
    - `fileName::String`: file name of the JLD file
"""
function load_voxel(fileName::String)
    global voxel = load(fileName, "voxel")
    empty!(idDict)
    idCount[] = 0
    grid_ind = sort(unique(voxel.grid))
    filter!(x -> x != 0, grid_ind)
    
    for ind in grid_ind
        idCount[] += 1
        idDict[idCount[]] = ind
    end
    nx = size(voxel.grid, 1)
    ny = size(voxel.grid, 2)
    nz = size(voxel.grid, 3)
    global gridID = Array{Vector}(undef, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        gridID[i, j, k] = []
        if voxel.grid[i, j, k] != 0
            ind = findfirst(x -> x .== voxel.grid[i, j, k], collect(values(idDict)))
            push!(gridID[i, j, k], collect(keys(idDict))[ind])
        end
    end
    
    return nothing
end


"""
    plot_voxel(addRef::Bool=true)

    Plots the voxel voxel. If `addRef=false`, the reference axes will not be added. Call this function if the plot window is accidentally closed.
    
    # Arguments
    - `addRef::Bool=true`: Boolean value to specify whether to add the reference axes to the plot.
"""
function plot_voxel(addRef::Bool=true)
    
    global canvas = plot([mesh3d(x=0, y=0, z=0)], blank_layout())
    display(canvas)
    
    _plot_voxel(gridID, addRef)
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
    global voxel = Voxels()
    return nothing
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
    global voxel = Voxels(zeros(1, 1, 1), dl, start)
    return nothing
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

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

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
    idDict[idCount[]] = ind
    _add_geom(geo, gridID)

    _plot_voxel(gridID, refAxis[])

    return geo
end

"""
    create_cube(origin::Vector{<:Real}, dim::Real, ind::Int=1, mode::String="corner", fac::Real=2)

    Creates a cube with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The origin point of the cube.
    - `dim::Real`: The side length of the cube.
    - `ind::Int=1`: The color index of the cube.
    - `mode::String="corner"`: The mode specifying the cube's origin ("corner" or "center").
    - `fac::Real=2`: The interior densified factor according to the grid spacing.
"""
function create_cube(origin::Vector{<:Real}, dim::Real, ind::Int=1, mode="corner", fac::Real=2)
    @assert dim > 0
    return create_cuboid(origin, [dim, dim, dim], ind, mode, fac)
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

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

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
    idDict[idCount[]] = ind
    _add_geom(geo, gridID)

    _plot_voxel(gridID, refAxis[])

    return geo
end

"""
    create_ellipsoid(origin::Vector{<:Real}, par::Vector{<:Real}, ind::Int=1, fac::Real=2)

    Creates an ellipsoid with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The origin point of the ellipsoid.
    - `par::Vector{<:Real}`: The lengths of the semi-axes.
    - `ind::Int=1`: The color index of the ellipsoid.
    - `fac::Real=2`: The interior densified factor according to the grid spacing.
"""
function create_ellipsoid(origin::Vector{<:Real}, par::Vector{<:Real}, ind::Int=1, fac::Real=2)
    @assert length(origin) == 3
    @assert length(par) == 3
    @assert fac > 0

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

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
    idDict[idCount[]] = ind
    _add_geom(geo, gridID)

    _plot_voxel(gridID, refAxis[])

    return geo
end

"""
    create_cylinder(origin::Vector{<:Real}, radius::Real, height::Real, ind::Int=1, fac::Real=2)

    Creates a cylinder with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The base origin point of the cylinder.
    - `radius::Real`: The radius of the cylinder.
    - `height::Real`: The height of the cylinder.
    - `ind::Int=1`: The color index of the cylinder.
    - `fac::Real=2`: The interior densified factor according to the grid spacing.
"""
function create_cylinder(origin::Vector{<:Real}, radius::Real, height::Real, ind::Int=1, fac::Real=2)
    @assert length(origin) == 3
    @assert fac > 0

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

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
    idDict[idCount[]] = ind
    _add_geom(geo, gridID)

    _plot_voxel(gridID, refAxis[])

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

    _del_geom(geo, gridID)

    for n in eachindex(geo.pos)
        geo.pos[n][1] += dl[1]
        geo.pos[n][2] += dl[2]
        geo.pos[n][1] += dl[3]
    end

    _add_geom(geo, gridID)

    _plot_voxel(gridID, refAxis[])
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

    _del_geom(geo, gridID)

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

    _add_geom(geo, gridID)

    _plot_voxel(gridID, refAxis[])
end

"""
    clear_geom(geo::Geometry)

    Removes the specified geometry from the voxel voxel.
    
    # Arguments
    - `geo::Geometry`: The geometry to be removed.
"""
function clear_geom(geo::Geometry)

    _del_geom(geo, gridID)
    geo = nothing
    _plot_voxel(gridID, refAxis[])
end

"""
    clear_geom(geoList::Vector{Geometry})

    Removes the specified list of geometries from the voxel voxel.
    
    # Arguments
    - `geoList::Vector{Geometry}`: The list of geometries to be removed.
"""
function clear_geom(geoList::Vector{Geometry})

    for i in eachindex(geoList)
        _del_geom(geoList[i], gridID)
        geoList[i] = nothing
    end
end

"""
    export_grid()

    Exports the grid array filled with color indexes. Note that when geometries overlap, the index of the last-added geometry is used.
    
    # Returns
    - `Array{Int}`: The grid array with color indexes.
"""
function export_grid()
    grid = zeros(Int, size(gridID))
    for i in eachindex(grid)
        if gridID[i] == []
            grid[i] = 0
        else
            grid[i] = idDict[gridID[i][end]]
        end
    end
    return grid
end


#endregion

