using Pkg, Revise
Pkg.activate("..")

using VoxelModel

reset_voxel()

colorDict[4] = "pink"
# create a sphere hollowed by cylinders
s = create_sphere([0, 0, 0], 10, 4)
c1 = create_cylin([0, 0, -10], 4, 20, 0)

c2 = create_cylin([0, 0, -10], 4, 20, 0)
rot!(c2, 90, [0, 1, 0])
c3 = create_cylin([0, 0, -10], 4, 20, 0)
rot!(c3, 90, [1, 0, 0])

c4 = create_cuboid([15, 0, 0], [1, 1, 1], 4)

