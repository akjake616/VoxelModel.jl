##region inner functionalitis
function _plot_voxel(addRef::Bool=true)

    if isnothing(canvas)
        global canvas = plot([mesh3d(x=0, y=0, z=0)], layout)
        display(canvas)
    else
        react!(canvas, [mesh3d(x=0, y=0, z=0)], layout)
    end
    if addRef
        add_ref_axes(canvas, [0, 0, 0], 1)
    end

    dx = space.dl[1]
    dy = space.dl[2]
    dz = space.dl[3]

    nx = size(space.gridID, 1)
    ny = size(space.gridID, 2)
    nz = size(space.gridID, 3)

    xmin = space.start[1]
    ymin = space.start[2]
    zmin = space.start[3]

    id_index = sort(unique(space.gridID))
    filter!(x -> x != [], id_index)

    for ind in eachindex(id_index)

        ptsArray = []
        for i in 1:nx, j in 1:ny, k in 1:nz

            if space.gridID[i, j, k] == id_index[ind]

                pts1 = [(i - 1.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 1.5) * dz + zmin]
                pts2 = [(i - 0.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 1.5) * dz + zmin]
                pts3 = [(i - 0.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 1.5) * dz + zmin]
                pts4 = [(i - 1.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 1.5) * dz + zmin]
                pts5 = [(i - 1.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 0.5) * dz + zmin]
                pts6 = [(i - 0.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 0.5) * dz + zmin]
                pts7 = [(i - 0.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 0.5) * dz + zmin]
                pts8 = [(i - 1.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 0.5) * dz + zmin]

                if i == 1 || space.gridID[i-1, j, k] == []
                    push!(ptsArray, pts1)
                    push!(ptsArray, pts4)
                    push!(ptsArray, pts8)
                    push!(ptsArray, pts5)
                end
                if i == nx || space.gridID[i+1, j, k] == []
                    push!(ptsArray, pts2)
                    push!(ptsArray, pts3)
                    push!(ptsArray, pts7)
                    push!(ptsArray, pts6)
                end

                if j == 1 || space.gridID[i, j-1, k] == []
                    push!(ptsArray, pts1)
                    push!(ptsArray, pts2)
                    push!(ptsArray, pts6)
                    push!(ptsArray, pts5)
                end
                if j == ny || space.gridID[i, j+1, k] == []
                    push!(ptsArray, pts4)
                    push!(ptsArray, pts3)
                    push!(ptsArray, pts7)
                    push!(ptsArray, pts8)
                end

                if k == 1 || space.gridID[i, j, k-1] == []
                    push!(ptsArray, pts1)
                    push!(ptsArray, pts2)
                    push!(ptsArray, pts3)
                    push!(ptsArray, pts4)
                end
                if k == nz || space.gridID[i, j, k+1] == []
                    push!(ptsArray, pts5)
                    push!(ptsArray, pts6)
                    push!(ptsArray, pts7)
                    push!(ptsArray, pts8)
                end
            end
        end
        voxel_obj = create_mesh(ptsArray, 4, colorDict[idDict[id_index[ind][end]]])

        addtraces!(canvas, voxel_obj)
        sleep(0.01)
    end
end

function _add_geom(geo::Geometry)

    setrounding(BigFloat, RoundToZero)

    dx = space.dl[1]
    dy = space.dl[2]
    dz = space.dl[3]

    np = length(geo.pos)
    xrange = getindex.(geo.pos, 1)
    yrange = getindex.(geo.pos, 2)
    zrange = getindex.(geo.pos, 3)

    if isempty(space.gridID)
        ngx = 1
        ngy = 1
        ngz = 1
    else
        ngx = size(space.gridID, 1)
        ngy = size(space.gridID, 2)
        ngz = size(space.gridID, 3)
    end

    xmin = _round((minimum([minimum(xrange), space.start[1]]) - shift_half*dx/2)/dx)*dx + shift_half*dx/2
    ymin = _round((minimum([minimum(yrange), space.start[2]]) - shift_half*dy/2)/dy)*dy + shift_half*dy/2
    zmin = _round((minimum([minimum(zrange), space.start[3]]) - shift_half*dz/2)/dz)*dz + shift_half*dz/2

    xmax = _round((maximum([maximum(xrange), space.start[1] + (ngx - 1) * dx]) - shift_half*dx/2)/dx)*dx + shift_half*dx/2
    ymax = _round((maximum([maximum(yrange), space.start[2] + (ngy - 1) * dy]) - shift_half*dy/2)/dy)*dy + shift_half*dy/2
    zmax = _round((maximum([maximum(zrange), space.start[3] + (ngz - 1) * dz]) - shift_half*dz/2)/dz)*dz + shift_half*dz/2

    x = collect(xmin:dx:xmax)
    y = collect(ymin:dy:ymax)
    z = collect(zmin:dz:zmax)

    nx = Int((xmax - xmin) / dx + 1)
    ny = Int((ymax - ymin) / dy + 1)
    nz = Int((zmax - zmin) / dz + 1)

    gridID_new = Array{Vector}(undef, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        gridID_new[i, j, k] = []
    end

    if !isempty(space.gridID)
        for i in 1:ngx, j in 1:ngy, k in 1:ngz
            gridID_new[findfirst(x .== space.start[1])+(i-1), findfirst(y .== space.start[2])+(j-1), findfirst(z .== space.start[3])+(k-1)] = space.gridID[i, j, k]
        end
    end

    for n in 1:np
        indx = _find_nearest(x, geo.pos[n][1])
        indy = _find_nearest(y, geo.pos[n][2])
        indz = _find_nearest(z, geo.pos[n][3])
        for i in eachindex(indx), j in eachindex(indy), k in eachindex(indz)
            if !(geo.ID in gridID_new[indx[i], indy[j], indz[k]])
                if idDict[geo.ID] != 0
                    push!(gridID_new[indx[i], indy[j], indz[k]], geo.ID)
                else 
                    gridID_new[indx[i], indy[j], indz[k]] = [];
                end
            end
        end
    end

    space.gridID = gridID_new

    space.start[1] = xmin
    space.start[2] = ymin
    space.start[3] = zmin
end

function _add_geom(geoList::Vector{Geometry})

    for i in eachindex(geoList)
        _add_geom(geoList[i])
    end
end

function _del_geom(geo::Geometry, trim::Bool=true)

    np = length(geo.pos)

    dx = space.dl[1]
    dy = space.dl[2]
    dz = space.dl[3]

    ngx = size(space.gridID, 1)
    ngy = size(space.gridID, 2)
    ngz = size(space.gridID, 3)

    x = collect(space.start[1]:dx:(ngx-1)*dx+space.start[1])
    y = collect(space.start[2]:dy:(ngy-1)*dy+space.start[2])
    z = collect(space.start[3]:dz:(ngz-1)*dz+space.start[3])

    for n in 1:np
        # filter!(e -> e != geo.ID, space.gridID[findfirst(x .== geo.pos[n][1]), findfirst(y .== geo.pos[n][2]), findfirst(z .== geo.pos[n][3])])
        # filter!(e -> e != geo.ID, space.gridID[_find_nearest(x, geo.pos[n][1]), _find_nearest(y, geo.pos[n][2]), _find_nearest(z, geo.pos[n][3])])
        indx = _find_nearest(x, geo.pos[n][1])
        indy = _find_nearest(y, geo.pos[n][2])
        indz = _find_nearest(z, geo.pos[n][3])
        for i in eachindex(indx), j in eachindex(indy), k in eachindex(indz)
            filter!(e -> e != geo.ID, space.gridID[indx[i], indy[j], indz[k]])
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
            if unique(space.gridID[i, :, :]) != []
                x1 = i
                break
            end
        end
        for i in ngx:1
            if unique(space.gridID[i, :, :]) != []
                x2 = i
                break
            end
        end
        for i in 1:ngy
            if unique(space.gridID[:, i, :]) != []
                y1 = i
                break
            end
        end
        for i in ngy:1
            if unique(space.gridID[:, i, :]) != []
                y2 = i
                break
            end
        end
        for i in 1:ngz
            if unique(space.gridID[:, :, i]) != []
                z1 = i
                break
            end
        end
        for i in ngz:1
            if unique(space.gridID[:, :, i]) != []
                z2 = i
                break
            end
        end

        space.gridID = space.gridID[x1:x2, y1:y2, z1:z2]
        space.start[1] = space.start[1] + (x1 - 1) * dx
        space.start[2] = space.start[2] + (y1 - 1) * dy
        space.start[3] = space.start[3] + (z1 - 1) * dz
    end
end

function _del_geom(geoList::Vector{Geometry}, trim::Bool=true)

    for i in eachindex(geoList)
        _del_geom(geoList[i], trim)
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