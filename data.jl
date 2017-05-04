using StatsBase
function convertDataToImgArray(data::Array)
    deltas = zeros(size(data)[1], size(data)[2])
    deltaKernel = [0, 1, -1];
    temp1 = conv(data[:,1], deltaKernel)[2:end-1]
    temp2 = conv(data[:,2], deltaKernel)[2:end-1]
    deltas[:,1] = abs.(temp1)
    deltas[:,2] = abs.(temp2)
    for i in 1:size(deltas)[1]
        if deltas[i, 1] == 0.0
            deltas[i,1] = NaN
        end
        if deltas[i, 2] == 0.0
            deltas[i,2] = NaN
        end
    end
    grads, idxs = findmin(deltas, 1)
    normedData = data./grads
    mins, ~ = findmin(normedData, 1)
    maxes, ~ = findmax(normedData, 1)
    histData::Tuple{Vector, Vector} = (normedData[:,1], normedData[:,2])
    numSteps = min.(abs.((maxes.-mins)./grads),1e4)
    @show numSteps
    imgVals = fit(Histogram, histData::Tuple{Vector, Vector} , (linspace(mins[1],maxes[1],trunc(Int64,numSteps[1])), linspace(mins[2],maxes[2],trunc(Int64,numSteps[2]))), closed=:left).weights

    img = imgVals .>= 1.0
    return img, mins, grads
end
