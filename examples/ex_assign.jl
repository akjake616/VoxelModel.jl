using Pkg, Revise
Pkg.activate("..")

using VoxelModel

A = zeros(Int, 10, 10, 10)
A[1:2, 1:3, 1:4] .= 1
A[4:7, 4:6, 4:5] .= 2

assign_voxel(A)
plot_voxel()
