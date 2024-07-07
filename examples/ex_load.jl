using Pkg, Revise
Pkg.activate("..")

using VoxelModel

load_voxel("hollowed_sphere_voxel.jld")
c1 = create_cube([0 ,0, 0], 2, 1, "center")
