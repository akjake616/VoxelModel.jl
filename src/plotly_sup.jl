using PlotlyJS
using Combinatorics
using LinearAlgebra
using Combinatorics

include("./meshgrid.jl")


# ******************** Easy 2D Plot ********************

function plot2d(y, x,  ylabel="", xlabel="", range=[0, 0], width=500, height=500, mode="lines")
    trace = scatter(y=y, x=x, mode=mode)
    layout = Layout(template=:plotly_white,
            width=width,
            height=height,
            yaxis=attr(
                title_text=ylabel,
                zeroline=false,
                showline=true,
                mirror=true,
                ticks="outside",
                tick0=minimum(y),
                automargin=true,
            ),
            xaxis=attr(
                title_text=xlabel,
                zeroline=false,
                showline=true,
                mirror=true,
                ticks="outside",
                tick0=minimum(x),
                automargin=true,
                range=[minimum(x), maximum(x)],
            ),
        )
        fig = plot(trace, layout)
        if range[1] !== 0 && range[2] !== 0
            update_yaxes!(
                fig, range=range
            )
        end
        display(fig)
end

# ******************** Export Cube ********************

function cubes(size, o, color, opc=1)
    # create points
    x1 = o[1] - size / 2
    x2 = o[1] + size / 2
    y1 = o[2] - size / 2
    y2 = o[2] + size / 2
    z1 = o[3] - size / 2
    z2 = o[3] + size / 2
    x = [x1, x1, x2, x2, x1, x1, x2, x2]
    y = [y1, y2, y2, y1, y1, y2, y2, y1]
    z = [z1, z1, z1, z1, z2, z2, z2, z2]
    i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2]
    j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3]
    k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6]

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        alphahull=1,
        flatshading=true,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0,
        ),
    )
end

function boxs(sizex, sizey, sizez, o, color, opc=1)
    # create points
    x1 = o[1] - sizex / 2
    x2 = o[1] + sizex / 2
    y1 = o[2] - sizey / 2
    y2 = o[2] + sizey / 2
    z1 = o[3] - sizez / 2
    z2 = o[3] + sizez / 2

    x = [x1, x1, x2, x2, x1, x1, x2, x2]
    y = [y1, y2, y2, y1, y1, y2, y2, y1]
    z = [z1, z1, z1, z1, z2, z2, z2, z2]
    i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2]
    j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3]
    k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6]

    # @infiltrate

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        alphahull=1,
        flatshading=true,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0,
        ),
    )
end

# ******************** Export Square ********************

function squares(size, o, color, mode="z", opc=1)
    if mode == "x"
        x = [o[1], o[1], o[1], o[1]]
        y = [o[2] - size / 2, o[2] + size / 2, o[2] + size / 2, o[2] - size / 2]
        z = [o[3] - size / 2, o[3] - size / 2, o[3] + size / 2, o[3] + size / 2]
        i = [0, 2]
        j = [1, 3]
        k = [2, 0]
    elseif mode == "y"
        x = [o[1] - size / 2, o[1] - size / 2, o[1] + size / 2, o[1] + size / 2]
        y = [o[2], o[2], o[2], o[2]]
        z = [o[3] - size / 2, o[3] + size / 2, o[3] + size / 2, o[3] - size / 2]
        i = [0, 2]
        j = [1, 3]
        k = [2, 0]
    else
        x = [o[1] - size / 2, o[1] - size / 2, o[1] + size / 2, o[1] + size / 2]
        y = [o[2] - size / 2, o[2] + size / 2, o[2] + size / 2, o[2] - size / 2]
        z = [o[3], o[3], o[3], o[3]]
        i = [0, 2]
        j = [1, 3]
        k = [2, 0]
    end

    return mesh3d(x=x, y=y, z=z, i=i, j=j, k=k,
        color=color,
        opacity=opc,
        alphahull=1,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

function rectangles(pts, color, opc=1)
    restore_pts!(pts)
    x = []
    y = []
    z = []
    for i in 1:4
        push!(x, pts[i][1])
        push!(y, pts[i][2])
        push!(z, pts[i][3])
    end

    i = [0, 2]
    j = [1, 3]
    k = [2, 0]
    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        # alphahull=1, 
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

function polygons(pts, color="yellow", opc=1)
    N = length(pts)
    # sort_pts!(pts, thtr, phir)
    sort_pts!(pts)
    x = []
    y = []
    z = []
    for i in eachindex(pts)
        push!(x, pts[i][1])
        push!(y, pts[i][2])
        push!(z, pts[i][3])
    end

    mid = zeros(3)
    for n in eachindex(pts)
        for m in 1:3
            mid[m] += pts[n][m] 
        end
    end
    mid = mid./N
    push!(x, mid[1])
    push!(y, mid[2])
    push!(z, mid[3])

    a = 0:1:length(pts)
    i = []
    j = []
    k = []
    # @infiltrate
    for n = 0:length(a)-2
        push!(i, a[mod(n + 0, length(a)-1)+1])
        push!(j, a[mod(n + 1, length(a)-1)+1])
        push!(k, a[end])
    end

    # @infiltrate

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        alphahull=1,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

function polygons(pts, color, thtr, phir, opc=1)
    N = length(pts)
    sort_pts!(pts, thtr, phir)
    x = []
    y = []
    z = []
    for i in eachindex(pts)
        push!(x, pts[i][1])
        push!(y, pts[i][2])
        push!(z, pts[i][3])
    end

    mid = zeros(3)
    for n in eachindex(pts)
        for m in 1:3
            mid[m] += pts[n][m] 
        end
    end
    mid = mid./N
    push!(x, mid[1])
    push!(y, mid[2])
    push!(z, mid[3])

    a = 0:1:length(pts)
    i = []
    j = []
    k = []
    # @infiltrate
    for n = 0:length(a)-2
        push!(i, a[mod(n + 0, length(a)-1)+1])
        push!(j, a[mod(n + 1, length(a)-1)+1])
        push!(k, a[end])
    end

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        alphahull=1,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

function create_mesh(pts, ng::Int, color="yellow", opc=1)
    N = length(pts)

    x = []
    y = []
    z = []
    i = []
    j = []
    k = []
    Ng = floor(Int, N/ng)
    for p = 1:Ng
        ptsg = []
        for m = 1:ng
            push!(ptsg, pts[(p-1)*ng + m])
        end
        sort_pts!(ptsg)
        for i in eachindex(ptsg)
            push!(x, ptsg[i][1])
            push!(y, ptsg[i][2])
            push!(z, ptsg[i][3])
        end

        mid = zeros(3)
        for n in eachindex(ptsg)
            for m in 1:3
                mid[m] += ptsg[n][m] 
            end
        end
        mid = mid./ng
        push!(x, mid[1])
        push!(y, mid[2])
        push!(z, mid[3])

        a = 0:1:length(ptsg)
        for n = 0:length(a)-2
            push!(i, a[mod(n + 0, length(a)-1)+1] + (p-1)*(ng+1))
            push!(j, a[mod(n + 1, length(a)-1)+1] + (p-1)*(ng+1))
            push!(k, a[end] + (p-1)*(ng+1))
        end
    end

    # @infiltrate

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        alphahull=1,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

function restore_pts!(pts) # restore rectangular for mesh3d
    if length(pts) > 3
        dis = zeros(3)
        for n in 2:4
            dis[n-1] = sum(abs2.(pts[n] .- pts[1]))
        end
        _, ind = findmax(dis)
        # @infiltrate
        if ind != 2
            tmp = pts[3]
            pts[3] = pts[ind+1]
            pts[ind+1] = tmp
        end
        # @infiltrate
    end
end

function sort_pts!(pts, thtr, phir) 
    N = length(pts)
    ang = zeros(length(pts))
    mid = zeros(3)
    for n in eachindex(pts)
        for m in 1:3
            mid[m] += pts[n][m] 
        end
    end
    mid = mid./N
    # Rx = [1  0           0         ;
    #       0  cosd(thtr) -sind(thtr);
    #       0  sind(thtr)  cosd(thtr)]
    Ry = [
        cosd(thtr) 0  -sind(thtr);   
        0  1  0;
        sind(thtr)  0  cosd(thtr)
    ]
    Rz = [
        cosd(phir)  sind(phir) 0;
        -sind(phir)   cosd(phir) 0;
         0            0          1
    ]
    R = Ry * Rz
    pts_rot = similar(pts)
    mid_rot = R * mid

    for n in eachindex(pts)
        pts_rot[n] = R * pts[n] 
        ang[n] = atan(pts_rot[n][2]-mid_rot[2], pts_rot[n][1]-mid_rot[1])
    end
    pts .= pts[sortperm(ang)]
end

function sort_pts!(pts) 
    N = length(pts)
    ang = zeros(length(pts))
    mid = zeros(3)
    for n in eachindex(pts)
        for m in 1:3
            mid[m] += pts[n][m] 
        end
    end
    mid = mid./N
    
    c = collect(combinations(1:N, 3))
    thtr = 0
    phir = 0
    for n in eachindex(c)
        vec = cross(pts[c[n][2]].-pts[c[n][1]], pts[c[n][3]].-pts[c[n][1]])
        if norm(vec) == 0
            continue
        else
            vec = vec./norm(vec)
            thtr = acosd(vec[3])
            phir = atand(vec[2], vec[1])
            break
        end
    end
    
    Ry = [
        cosd(thtr) 0  -sind(thtr);   
        0  1  0;
        sind(thtr)  0  cosd(thtr)
    ]
    Rz = [
        cosd(phir)  sind(phir) 0;
        -sind(phir)   cosd(phir) 0;
         0            0          1
    ]
    R = Ry * Rz
    pts_rot = similar(pts)
    mid_rot = R * mid
    for n in eachindex(pts)
        pts_rot[n] = R * pts[n] 
        ang[n] = atan(pts_rot[n][2]-mid_rot[2], pts_rot[n][1]-mid_rot[1])
    end
    pts .= pts[sortperm(ang)]
end

# ******************** Export Line ********************

function lines(pt1, pt2, color, opc=1, style="")

    x = [pt1[1], pt2[1]]
    y = [pt1[2], pt2[2]]
    z = [pt1[3], pt2[3]]

    return scatter3d(x=x, y=y, z=z,
        mode="lines",
        line=attr(
            color=color,
            dash=style,
        ),
        showlegend=false,
        alphahull=1,
        flatshading=true,
        opacity=opc,
    )
end

# ******************** Export Arrow ********************

function arrows(o, dir, color, opc=1)
    c = cone(x=[o[1] + dir[1] / 2], y=[o[2] + dir[2] / 2], z=[o[3] + dir[3] / 2], u=[dir[1]], v=[dir[2]], w=[dir[3]],
        colorscale=[[0, color], [1, color]],
        opacity=opc,
        showscale=false)

    l = scatter3d(x=[o[1] - dir[1] / 2, o[1] + dir[1] / 2], y=[o[2] - dir[2] / 2, o[2] + dir[2] / 2], z=[o[3] - dir[3] / 2, o[3] + dir[3] / 2],
        line=attr(color=color, width=4),
        mode="lines",
        opacity=opc,
        showlegend=false)

    return [c, l]
end

# ******************** Export Ellipsoid ********************

function ellipsoids(par, o, rotang, color, opc=1, res=25)
    a = par[1]
    b = par[2]
    c = par[3]

    alpha = rotang[1]
    beta = rotang[2]
    gama = rotang[3]

    Rx = [1 0 0;
        0 cosd(alpha) -sind(alpha);
        0 sind(alpha) cosd(alpha)]
    Ry = [cosd(beta) 0 sind(beta);
        0 1 0;
        -sind(beta) 0 cosd(beta)]
    Rz = [cosd(gama) -sind(gama) 0;
        sind(gama) cosd(gama) 0;
        0 0 1]

    R = Rx * Ry * Rz

    # create points
    P, T = meshgrid(
        LinRange(0, 2 * pi, res),
        LinRange(0, pi, res),
    )

    x = sin.(T) .* cos.(P) .* a
    y = sin.(T) .* sin.(P) .* b
    z = cos.(T) .* c
    x = x[:]
    y = y[:]
    z = z[:]

    @inbounds for n in eachindex(x)
        vec = [x[n], y[n], z[n]]
        vec = R * vec + o
        x[n] = vec[1]
        y[n] = vec[2]
        z[n] = vec[3]
    end

    return mesh3d(x=x, y=y, z=z,
        alphahull=0,
        flatshading=true,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=2.0,
            roughness=0.5
        ),
    )
end

# ******************** Add XYZ Axes ********************

function add_ref_axes(plt, o, r)
    cx = cone(x=[r + o[1]], y=[o[2]], z=[o[3]], u=[r / 10], v=[0], w=[0],
        colorscale=[[0, "red"], [1, "red"]],
        showscale=false)
    cy = cone(x=[o[1]], y=[r + o[2]], z=[o[3]], u=[0], v=[r / 10], w=[0],
        colorscale=[[0, "green"], [1, "green"]],
        showscale=false)
    cz = cone(x=[o[1]], y=[o[2]], z=[r + o[3]], u=[0], v=[0], w=[r / 10],
        colorscale=[[0, "blue"], [1, "blue"]],
        showscale=false)
    lx = scatter3d(x=[o[1], r + o[1]], y=[o[2], o[2]], z=[o[3], o[3]],
        line=attr(color="red"),
        mode="lines",
        showlegend=false)
    ly = scatter3d(x=[o[1], o[1]], y=[o[2], r + o[2]], z=[o[3], o[3]],
        line=attr(color="green"),
        mode="lines",
        showlegend=false)
    lz = scatter3d(x=[o[1], o[1]], y=[o[2], o[2]], z=[o[3], r + o[3]],
        line=attr(color="blue"),
        mode="lines",
        showlegend=false)
    addtraces!(plt, cx)
    addtraces!(plt, cy)
    addtraces!(plt, cz)
    addtraces!(plt, lx)
    addtraces!(plt, ly)
    addtraces!(plt, lz)
    relayout!(plt, scene=attr(
        annotations=[
            attr(
                showarrow=false,
                x=o[1] + 1.1 * r,
                y=o[2],
                z=o[3],
                text="x",
                font=attr(color="red")
            ),
            attr(
                showarrow=false,
                x=o[1],
                y=o[2] + 1.1 * r,
                z=o[3],
                text="y",
                font=attr(color="green")
            ),
            attr(
                showarrow=false,
                x=o[1],
                y=o[2],
                z=o[3] + 1.1 * r,
                text="z",
                font=attr(color="blue")
            ),
        ]
    ))
end

# ******************** Turn Off Grids ********************

function no_grids(plt)
    relayout!(plt,
        scene=attr(
            xaxis=attr(
                visible=false,
                showgrid=false
            ),
            yaxis=attr(
                visible=false,
                showgrid=false
            ),
            zaxis=attr(
                visible=false,
                showgrid=false
            ),
        ))
end

function blank_layout()
    layout = Layout(
        scene=attr(
            xaxis=attr(
                visible=false,
                showgrid=false
            ),
            yaxis=attr(
                visible=false,
                showgrid=false
            ),
            zaxis=attr(
                visible=false,
                showgrid=false
            ),
            aspectmode="data",
        ),
    )
    return layout
end