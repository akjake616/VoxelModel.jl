function meshgrid(x, y)
    sx = size(x, 1)
    sy = size(y, 1)

    T = typeof(x[1])

    X = zeros(T, sy, sx)
    Y = zeros(T, sy, sx)

    @inbounds for i = 1:sx, j = 1:sy
        X[j, i] = x[i]
        Y[j, i] = y[j]
    end

    return (X, Y)
end

function meshgrid(x, y, z)
    sx = size(x, 1)
    sy = size(y, 1)
    sz = size(z, 1)

    T = typeof(x[1])

    X = zeros(T, sy, sx, sz)
    Y = zeros(T, sy, sx, sz)
    Z = zeros(T, sy, sx, sz)

    @inbounds for i = 1:sx, j = 1:sy, k = 1:sz
        X[j, i, k] = x[i]
        Y[j, i, k] = y[j]
        Z[j, i, k] = z[k]
    end

    return (X, Y, Z)
end