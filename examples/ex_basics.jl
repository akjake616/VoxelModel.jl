## using the modeule

push!(LOAD_PATH, "../src")
using VoxelModel
using Infiltrator

reset_voxel()

## create a cube with base at [3, 4, 0] and dimensions [1, 2, 3] (color index = 1)
c1 = create_cube([0, 0, 0], [1, 2, 3], 1)

## create a sphere with center at [0, 0, 8] and radius = 5 (color index = 2)
s1 = create_sphere([0, 0, 8], 5, 2)

## create a ellipsoid with center at [0, 10, 0] and the semi-axes [3, 4, 5] (color index = 3)
e1 = create_ellip([0, 10, 0], [3, 4, 5], 3)

## wait 5 secs
for i in 1:5
    println("translate cube after $(6-i) sec...")
    sleep(1)
end

## translate c1 with [8, 1, 3]
trans!(c1, [8, 1, 3])


## wait 5 secs
for i in 1:5
    println("delete cube after $(6-i) sec...")
    sleep(1)
end

# remove geometries
clear_geom(c1)

# export grid array with geometry indexes
grid = export_grid()
