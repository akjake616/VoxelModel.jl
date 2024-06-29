module VoxelModel

export reset_voxel, reset_shift, reset_dl
export create_cuboid, create_sphere, create_ellip, create_cylin
export trans!, rot!
export clear_geom, plot_voxel, export_voxel, export_grid
export colorDict, canvas

using PlotlyJS
using PlotlyGeometries
using LinearAlgebra
using Combinatorics

#region structs
mutable struct Geometry
    pos::Vector{Vector{Float64}}
    index::Int
    const ID::Int
end

Base.@kwdef mutable struct Voxels
    gridID::Array{Vector} = []
    dl::Vector{Float64} = [1.0, 1.0, 1.0]
    start::Vector{Float64} = [shift_half[] * 0.5, shift_half[] * 0.5, shift_half[] * 0.5]
end
#endregion

#region constrols
const shift_half = Ref(true)
const idCount = Ref{Int}(0)
const idDict = Ref(Dict{Int, Int}())
const ref = Ref(true)

global space = Voxels()
global canvas = nothing
#endregion

include("./color_dict.jl")
include("./api.jl")
include("./internal.jl")

end
