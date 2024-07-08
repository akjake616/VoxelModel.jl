using Pkg, Revise
Pkg.activate(".")

using VoxelModel

reset_dl(3.33102731122781e-8)
create_sphere([0, 0, 0], 5e-7, 6)