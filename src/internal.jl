##region internal functionalities
function _plot_voxel(gridID::Array{Vector}, addRef::Bool=true)
    if isnothing(canvas)
        global canvas = plot([mesh3d(x=0, y=0, z=0)], blank_layout())
        display(canvas)
    else
        react!(canvas, [mesh3d(x=0, y=0, z=0)], blank_layout())
    end
    if addRef
        add_ref_axes!(canvas, [0, 0, 0], refLen[])
    end

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

    nx = size(gridID, 1)
    ny = size(gridID, 2)
    nz = size(gridID, 3)

    xmin = voxel.start[1]
    ymin = voxel.start[2]
    zmin = voxel.start[3]

    id_list = sort(unique(gridID))
    filter!(x -> x != [], id_list)

    @all pts1 pts2 pts3 pts4 pts5 pts6 pts7 pts8 = fill(0.0, 3)
    @all r g b = 0.0
    for ind in eachindex(id_list)
        if idDict[id_list[ind][end]] != 0
            ptsArray = []
            for i in 1:nx, j in 1:ny, k in 1:nz
                if gridID[i, j, k] == id_list[ind]
                    pts1 = [(i - 1.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 1.5) * dz + zmin]
                    pts2 = [(i - 0.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 1.5) * dz + zmin]
                    pts3 = [(i - 0.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 1.5) * dz + zmin]
                    pts4 = [(i - 1.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 1.5) * dz + zmin]
                    pts5 = [(i - 1.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 0.5) * dz + zmin]
                    pts6 = [(i - 0.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 0.5) * dz + zmin]
                    pts7 = [(i - 0.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 0.5) * dz + zmin]
                    pts8 = [(i - 1.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 0.5) * dz + zmin]
    
                    if i == 1 || gridID[i-1, j, k] == [] || idDict[gridID[i-1, j, k][end]] == 0
                        push!(ptsArray, pts1)
                        push!(ptsArray, pts4)
                        push!(ptsArray, pts8)
                        push!(ptsArray, pts5)
                    end
                    if i == nx || gridID[i+1, j, k] == [] || idDict[gridID[i+1, j, k][end]] == 0
                        push!(ptsArray, pts2)
                        push!(ptsArray, pts3)
                        push!(ptsArray, pts7)
                        push!(ptsArray, pts6)
                    end
                    if j == 1 || gridID[i, j-1, k] == [] || idDict[gridID[i, j-1, k][end]] == 0
                        push!(ptsArray, pts1)
                        push!(ptsArray, pts2)
                        push!(ptsArray, pts6)
                        push!(ptsArray, pts5)
                    end
                    if j == ny || gridID[i, j+1, k] == [] || idDict[gridID[i, j+1, k][end]] == 0
                        push!(ptsArray, pts4)
                        push!(ptsArray, pts3)
                        push!(ptsArray, pts7)
                        push!(ptsArray, pts8)
                    end
                    if k == 1 || gridID[i, j, k-1] == [] || idDict[gridID[i, j, k-1][end]] == 0
                        push!(ptsArray, pts1)
                        push!(ptsArray, pts2)
                        push!(ptsArray, pts3)
                        push!(ptsArray, pts4)
                    end
                    if k == nz || gridID[i, j, k+1] == [] || idDict[gridID[i, j, k+1][end]] == 0
                        push!(ptsArray, pts5)
                        push!(ptsArray, pts6)
                        push!(ptsArray, pts7)
                        push!(ptsArray, pts8)
                    end
                end
            end
            if !haskey(colorDict, idDict[id_list[ind][end]])
                @all r g b = round(Int, rand() * 255)
                colorDict[idDict[id_list[ind][end]]] = "rgb($r, $g, $b)"
            end
            voxel_obj = polygons(ptsArray, 4, colorDict[idDict[id_list[ind][end]]])
    
            addtraces!(canvas, voxel_obj)
            sleep(0.1)
        end
    end
end

function _add_geom(geo::Geometry, gridID::Array{Vector})
    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

    np = length(geo.pos)
    xrange = getindex.(geo.pos, 1)
    yrange = getindex.(geo.pos, 2)
    zrange = getindex.(geo.pos, 3)

    if isempty(gridID)
        @all ngx ngy ngz = 1
    else
        ngx = size(gridID, 1)
        ngy = size(gridID, 2)
        ngz = size(gridID, 3)
    end

    xmin = _round((minimum([minimum(xrange), voxel.start[1]]) - shift[]*dx/2)/dx)*dx + shift[]*dx/2
    ymin = _round((minimum([minimum(yrange), voxel.start[2]]) - shift[]*dy/2)/dy)*dy + shift[]*dy/2
    zmin = _round((minimum([minimum(zrange), voxel.start[3]]) - shift[]*dz/2)/dz)*dz + shift[]*dz/2

    xmax = _round((maximum([maximum(xrange), voxel.start[1] + (ngx - 1) * dx]) - shift[]*dx/2)/dx)*dx + shift[]*dx/2
    ymax = _round((maximum([maximum(yrange), voxel.start[2] + (ngy - 1) * dy]) - shift[]*dy/2)/dy)*dy + shift[]*dy/2
    zmax = _round((maximum([maximum(zrange), voxel.start[3] + (ngz - 1) * dz]) - shift[]*dz/2)/dz)*dz + shift[]*dz/2

    x = collect(xmin:dx:xmax)
    y = collect(ymin:dy:ymax)
    z = collect(zmin:dz:zmax)

    nx = round(Int, (xmax - xmin) / dx + 1)
    ny = round(Int, (ymax - ymin) / dy + 1)
    nz = round(Int, (zmax - zmin) / dz + 1)
    
    gridID_new = Array{Vector}(undef, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        gridID_new[i, j, k] = []
    end

    if !isempty(gridID)
        for i in 1:ngx, j in 1:ngy, k in 1:ngz
            gridID_new[findfirst(x .== voxel.start[1])+(i-1), findfirst(y .== voxel.start[2])+(j-1), findfirst(z .== voxel.start[3])+(k-1)] = gridID[i, j, k]
        end
    end

    for n in 1:np
        indx = _find_nearest(x, geo.pos[n][1])
        indy = _find_nearest(y, geo.pos[n][2])
        indz = _find_nearest(z, geo.pos[n][3])
        for i in eachindex(indx), j in eachindex(indy), k in eachindex(indz)
            if !(geo.ID in gridID_new[indx[i], indy[j], indz[k]])
                push!(gridID_new[indx[i], indy[j], indz[k]], geo.ID)
            end
        end
    end

    global gridID = gridID_new

    voxel.start[1] = xmin
    voxel.start[2] = ymin
    voxel.start[3] = zmin
    
    voxel.grid = export_grid()
    
    return nothing
end

function _add_geom(geoList::Vector{Geometry}, gridID::Array{Vector})
    for i in eachindex(geoList)
        _add_geom(geoList[i], gridID)
    end
end

function _del_geom(geo::Geometry, gridID::Array{Vector}, trim::Bool=true)
    np = length(geo.pos)

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

    ngx = size(gridID, 1)
    ngy = size(gridID, 2)
    ngz = size(gridID, 3)

    x = collect(voxel.start[1]:dx:(ngx-1)*dx+voxel.start[1])
    y = collect(voxel.start[2]:dy:(ngy-1)*dy+voxel.start[2])
    z = collect(voxel.start[3]:dz:(ngz-1)*dz+voxel.start[3])

    for n in 1:np
        indx = _find_nearest(x, geo.pos[n][1])
        indy = _find_nearest(y, geo.pos[n][2])
        indz = _find_nearest(z, geo.pos[n][3])
        for i in eachindex(indx), j in eachindex(indy), k in eachindex(indz)
            filter!(e -> e != geo.ID, gridID[indx[i], indy[j], indz[k]])
        end
    end

    if trim
        x1 = 1
        x2 = ngx
        y1 = 1
        y2 = ngy
        z1 = 1
        z2 = ngz

        for i in 1:ngx
            if unique(gridID[i, :, :]) != []
                x1 = i
                break
            end
        end
        for i in ngx:1
            if unique(gridID[i, :, :]) != []
                x2 = i
                break
            end
        end
        for i in 1:ngy
            if unique(gridID[:, i, :]) != []
                y1 = i
                break
            end
        end
        for i in ngy:1
            if unique(gridID[:, i, :]) != []
                y2 = i
                break
            end
        end
        for i in 1:ngz
            if unique(gridID[:, :, i]) != []
                z1 = i
                break
            end
        end
        for i in ngz:1
            if unique(gridID[:, :, i]) != []
                z2 = i
                break
            end
        end

        gridID = gridID[x1:x2, y1:y2, z1:z2]
        voxel.start[1] = voxel.start[1] + (x1 - 1) * dx
        voxel.start[2] = voxel.start[2] + (y1 - 1) * dy
        voxel.start[3] = voxel.start[3] + (z1 - 1) * dz
        
        voxel.grid = export_grid()
        
        return nothing
    end
end

function _del_geom(geoList::Vector{Geometry}, gridID::Array{Vector}, trim::Bool=true)
    for i in eachindex(geoList)
        _del_geom(geoList[i], gridID, trim)
    end
end

function _reset_gridID()
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
end
function _find_nearest(ary::Array{<:Number}, ele::Number)
    tmp = abs.(ary .- ele)
    val, ind = findmin(tmp)
    ind = findall(tmp .== val)
    return ind
end

function _round(num::Number)
    if num - floor(num) == 0.5
        if num < 0
            return floor(num)     
        else
            return ceil(num)
        end
    else
        return round(num)
    end
end
#endregion