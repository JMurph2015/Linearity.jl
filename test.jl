using ExcelReaders
include("./linearity.jl")
function readData(run::String)
  f = openxl("20161121 Group 7 TensileTestProgram_1.xls")
  indata = readxlsheet(f, "Group 7$run", skipstartrows=3)
  data = transpose(convert(Array{Float64,2},indata))
  #=
  temp1 = "Group 7$(run)!A4:A$top"
  temp2 = "Group 7$(run)!B4:B$top"
  data1 = readxl(f, temp1)
  data2 = readxl(f, temp2)
  outdata1 = convert(Array{Float64,2}, data1)
  outdata2 = convert(Array{Float64,2}, data2)
  data = transpose(cat(2, outdata1, outdata2))
  =#
  return data
end
function extractYoungs(run::String)
  areas = Dict(
  "A" =>  9.393e-6,
  "B" =>  9.731e-6,
  "C" =>  9.3635e-6,
  "D" =>  9.66e-6,
  )
  linearity = Dict(
  "A" => 200:950,
  "B" => 200:935,
  "C" => 250:1070,
  "D" => 600:2500,
  )
  using PyPlot
  rawdata = readData(run)
  data = transpose(cat(2, rawdata[1,:]/178, rawdata[2,:]/areas[run]))
  @show size(data)
  #derivArray = approxDerivative(data)
  #lineardata = transpose(cat(2, data[1,linearity[run]], data[2,linearity[run]], derivArray[linearity[run]]))
  lineardata = transpose(findLinearity(transpose(data),0.01,50))
  @show size(lineardata)
  title("Material $run")
  xlabel("Strain")
  ylabel("Stress (Pa)")
  plot(data[1,:],data[2,:], color="red", linewidth=2.0, linestyle="--")
  plot(lineardata[1,:],lineardata[2,:], color="blue", linewidth=2.0, linestyle=":")
  show()
  #youngs = (sum(lineardata[3,:],1)/size(lineardata,2))[1]
  #printout = "Young's Modulus for Material $run: $youngs"
  #println(printout)
end
