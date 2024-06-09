module VoxelModel

export reset_voxel, reset_shift, reset_dl, 
    export_voxel, export_grid,
    create_cube, create_sphere, create_ellip, create_cylin,
    trans!, rot!,
    clear_geom, plot_voxel,
    colorDict,
    canvas

using PlotlyJS
# using PlotlyGeometries
using LinearAlgebra
using Combinatorics

include("./plotly_sup.jl")
include("./color_dict.jl")

include("./apis.jl")
include("./inner_funcs.jl")


end
