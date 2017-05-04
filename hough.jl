using StatsBase
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
        accumulator[i,:] .= cosine[xIndicies[i],:] .+ sine[yIndicies[i],:];
    end

    # accumulator[1:numEdgePixels,:] = @. (xIndicies - 1)*cos(theta)' + (yIndicies - 1)*sin(theta)';
    #Scan over the thetas and bin the rhos
    for i in 1:numThetas
        houghSpace[:,i] .= fit(Histogram, accumulator[:,i], rho, closed=:left).weights
    end
    return (theta, rho , houghSpace)
end

function weightedHough(data::AbstractArray, numTheta, numRho, thetaLimits=[0, pi])
    numData = size(data)[1]
    maxX = findmax(data[:,1])[1]
    maxY = findmax(data[:,2])[1]
    rhoLimit = norm([maxX maxY])
    theta = linspace(thetaLimits[1],thetaLimits[2],numTheta)
    rho = linspace(-rhoLimit, rhoLimit, numRho)
    accumulator = zeros(numData, numTheta)
    houghSpace = zeros(numTheta-1, numTheta)
    sines = sin.(theta)'
    cosines = cos.(theta)'

    @. accumulator[1:numData,:] = (data[:,1] - 1)*cosines + (data[:,2] - 1)*sines;

    for i in eachindex(theta)
        houghSpace[:,i] .= fit(Histogram, accumulator[:,i], rho, closed=:left).weights
    end

    return (theta, rho , houghSpace)
end

function test()
    using Iterators
    using AxisArrays
    using Images, Colors, ColorTypes, PyPlot
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
