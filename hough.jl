using Iterators
using AxisArrays
using Images, Colors, ColorTypes, PyPlot, StatsBase
#=
function hist(x, y)
    h = fit(Histogram, x, y)
    return h.weights
end
=#
function hist(x::Array{Float64, 2}, y::Range)
    outs = zeros(length(y))
    edges = y
    loops = zeros(5)
    @fastmath @inbounds begin
        for i in eachindex(x)
            if x[i] < edges[2]
                outs[1] += 1
                loops[4] += 1
            elseif x[i] > edges[end-1]
                outs[end] += 1
                loops[5] += 1
            else
                for j in 2:length(y)-1
                    if j == 1 && x[i] < edges[j]
                        outs[j] += 1
                        loops[1] += 1
                    elseif j == length(edges) && x[i] > edges[j]
                        outs[j] += 1
                        loops[2] += 1
                    elseif x[i] >= edges[j] && x[i] < edges[j+1]
                        outs[j] += 1
                        loops[3] += 1
                    end
                end
            end
        end
    end
    return outs
end
#= function hist!(x::Array{Float64}, y::Range, outs::Array{Float64})

    length(y) <= length(outs) || error("Cannot assign more bins than the length of the output.")

    edges = y
    @fastmath @inbounds begin
        for i in eachindex(x)
            if x[i] < edges[2]
                outs[1] += 1
            elseif x[i] > edges[end-1]
                outs[end] += 1
            else
                for j in 2:length(y)-1
                    if j == 1 && x[i] < edges[j]
                        outs[j] += 1
                    elseif j == length(edges) && x[i] > edges[j]
                        outs[j] += 1
                    elseif x[i] >= edges[j] && x[i] < edges[j+1]
                        outs[j] += 1
                    end
                end
            end
        end
    end
end
=#
function hist!(x::Array{Float64}, y::Range, outs::Array{Float64})
    outs = fit(Histogram, x, y, closed=:right).weights
end

function histRef(x, y::Range)
    outs = zeros(length(y))
    edges = y
    outs[1:end] .= sum(x' .< edges, 2)[:]
    @. outs[2:end] = outs[2:end] - outs[1:end-1]
    return outs
end
function testHist()
    temp = true
    for i in 1:2000
        bins = rand()*10:0.1:20
        t3 = rand(100)*15
        temp = temp && histRef(t3,bins) == hist(t3,bins)
    end
    return temp
end
precompile(hist, (Array{Float64, 2}, Range))
#precompile(hist!, (Array{Float64}, Range, Array{Float64}))
function hist(x, y::AbstractArray)
    edges = @. (y[1:end-1]+y[2:end])/2
    outs = zeros(length(y))
    outs[1:end-1] = sum(x' .< edges, 2)
    outs[2:end] = outs[2:end] .- outs[1:end-1]
    outs[end] = length(x) - sum(outs[1:end-1])
    return outs
end
function hough(theImage::AbstractArray{Bool, 2}, thetaSampleFrequency)
    theImage = flipdim(theImage,1);
    width, height = size(theImage);

    rhoLimit = norm([width height]);
    rho = -rhoLimit:1:rhoLimit;
    theta = 0:thetaSampleFrequency:pi;

    numThetas = length(theta);
    houghSpace = zeros(length(rho)-1,numThetas);

    idxs = find(theImage)
    xIndicies, yIndicies = ind2sub(theImage, idxs);
    cartidx = size([xIndicies yIndicies])

    #Preallocate space for the accumulator array
    numEdgePixels = length(xIndicies);
    accumulator = zeros(numEdgePixels,numThetas);

    #Preallocate cosine and sine calculations to increase speed. In
    #addition to precallculating sine and cosine we are also multiplying
    #them by the proper pixel weights such that the rows will be indexed by
    #the pixel number and the columns will be indexed by the thetas.
    #Example: cosine(3,:) is 2*cosine(0 to pi)
    #         cosine(:,1) is (0 to width of image)*cosine(0)
    cosine = @. [i for i in 0:width-1]*cos(theta)'; #Matrix Outerproduct
    sine = @. [i for i in 0:height-1]*sin(theta)'; #Matrix Outerproduct

    for i in 1:numEdgePixels
        accumulator[i,:] = cosine[xIndicies[i],:] + sine[yIndicies[i],:];
    end

    # accumulator[1:numEdgePixels,:] = @. (xIndicies - 1)*cos(theta)' + (yIndicies - 1)*sin(theta)';
    #Scan over the thetas and bin the rhos
    tic()
    for i in 1:numThetas
        houghSpace[:,i] = fit(Histogram, accumulator[:,i], rho, closed=:left).weights
    end
    toc()
    return (theta, rho , houghSpace)
end

function test()
    image = load("Pentagon.png")
    bwImg = convert.(Gray,image)
    edges = canny(bwImg)
    theta, rho, houghSpace, = hough(edges,0.05)
    pcolormesh(Array(theta),Array(rho),Array(houghSpace));
    title("Hough Transform");
    xlabel("Theta (radians)");
    ylabel("Rho (pixels)");
    show()
end
# test()
# y = (p/sin(theta)-x*cot(theta))
