## using the modeule

push!(LOAD_PATH, "../src")
using VoxelModel
using Infiltrator

reset_voxel()

s = create_sphere([0, 0, 0], 10, 2)
c1 = create_cylin([0, 0, -10], 4, 20, 0)
c2 = create_cylin([0, 0, -10], 4, 20, 0)
rotate!(c2, 90, [0, 1, 0])
c3 = create_cylin([0, 0, -10], 4, 20, 0)
rotate!(c3, 90, [1, 0, 0])

