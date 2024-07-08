using VoxelModel

reset_voxel(render=true)

colorDict[4] = "pink"
# create a sphere hollowed by cylinders
s = create_sphere([0, 0, 0], 10, 4; render=true)
c1 = create_cylinder([0, 0, -10], 4, 20, 0; render=true)

c2 = create_cylinder([0, 0, -10], 4, 20, 0)
rot!(c2, 90, [0, 1, 0]; render=true)
c3 = create_cylinder([0, 0, -10], 4, 20, 0)
rot!(c3, 90, [1, 0, 0]; render=true)

save_voxel("hollowed_sphere.jld")