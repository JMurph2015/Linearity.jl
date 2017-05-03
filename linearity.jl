#=
    A few useful functions for using linear portions of graphs
    Author: Joseph Murphy
    Date: 05-01-2017
    Version: 0.1.0
=#
function findLinearity(dataset::AbstractArray, target::Real, minSize::Real)
    # objective function f(x) = abs(target - rmse(subset(x)))
    n = size(dataset)[1]
    @show n
    f(x) = abs(target - rmse(continuousSubset(dataset, x)))
    p = n - minSize
    numSubsets = sum((i for i in 1:p+1))
    subsets = 1:numSubsets
    @show subsets
    outputs = f.(subsets)
    best = [1e5, 0]
    idx = 1
    for i in eachindex(outputs)
        j, k = rollingModulo(subsets[i]-1, n)
        sublength = k-j
        score = [outputs[i], sublength]
        if outputs[i] == best[1] && sublength > best[2]
            best = [outputs[i], sublength]
            idx = i
        end
        if outputs[i] < best[1]
            best = [outputs[i], sublength]
            idx = i
        end
    end
    @show best
    @show idx
    return continuousSubset(dataset, subsets[idx])
end

function rollingModulo(i,j)
    # analogous to i % j
    k = j
    while i >= j && j >= 1
        i -= j
        j -= 1
    end
    if j <=0
        return -1, 0
    end
    return i, j
end

function continuousSubset(data, n)
    g(x,n) = rollingModulo(sum(1:n)-x, n)
    i,j = g(n, size(data)[1])
    if i <= -1
        ValueError("Invalid Index")
    end
    return data[(i+1):end-(j-i-1),:]
end

function rmse(data)
    b, m = linreg(data[:,1], data[:,2])
    if b == NaN || m === NaN
        return 1e5
    end
    mse = sum([(data[i,2] - (m*data[i,1]+b))^2 for i in eachindex(data[:,1])])/length(data[:,1])
    rmse = sqrt(mse)
    return rmse
end
