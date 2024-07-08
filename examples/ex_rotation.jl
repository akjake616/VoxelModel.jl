using VoxelModel

reset_voxel(render=true)

## create a cuboid with base at [3, 4, 0] and dimensions [1, 2, 3] (color index = 1)
c1 = create_cuboid([0, 0, 0], [30, 20, 10], 1; render=true)

## wait 5 secs
for i in 1:5
    println("rotate cuboid after $(6-i) sec...")
    sleep(1)
end

## rotate 90 deg respect to the x axis ([1, 0, 0]) according to the geometry center
rot!(c1, 90, [1, 0, 0]; render=true)

## wait 5 secs
for i in 1:5
    println("rotate cuboid after $(6-i) sec...")
    sleep(1)
end

## rotate 30 deg respect to the y axis ([0, 1, 0]) according to the origin ([0, 0, 0])
rot!(c1, 30, [0, 1, 0], [0, 0, 0]; render=true)


