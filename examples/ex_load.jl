using Pkg, Revise
Pkg.activate("..")

using VoxelModel

reset_voxel()

load_voxel("hollowed_sphere_voxel.jld")
c1 = create_cuboid([0 ,0, 0], [1, 1, 1], 1)