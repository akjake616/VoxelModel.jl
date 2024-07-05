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
    start::Vector{Float64} = [shift[] * 0.5, shift[] * 0.5, shift[] * 0.5]
end
#endregion

#region constrols
const shift = Ref(true)
const refAxis = Ref(true)
const idCount = Ref{Int}(0)
const idDict = Ref(Dict{Int, Int}())

global space = Voxels()
global canvas = nothing
const colorDict = Dict{Int, String}()
#endregion

include("api.jl")
include("internal.jl")

end
