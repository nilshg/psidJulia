################################################################################
#####             PROFILE HETEROGENEITY AND INCOME RISK                    #####
################################################################################
nprocs()==CPU_CORES || addprocs(CPU_CORES-1)
import HDF5, NLopt

@everywhere begin
  include("C:/Users/tew207/My Documents/GitHub/psidJulia/theoretical_varcov.jl")
  include("C:/Users/tew207/My Documents/GitHub/psidJulia/estimate.jl")

  tinit = 1992; tlast = 2009; tmax = tlast - tinit + 1
  data = "C:/Users/tew207/My Documents/GitHub/BHPStools/BHPS_logrnetincdef_1992_2009.h5"
end

methods = [:GN_CRS2_LM, :GN_MLSL, :GN_ISRES, :GN_ESCH]
l = length(methods)
HIPresults = SharedArray(Float64, (l, 6+2*tmax+1))
RIPresults = SharedArray(Float64, (l, 6+2*tmax+1))

@sync @parallel for i = 1:length(methods)
  println("Fitting HIP process with ",string(methods[i]))
  (hipf, hipx, flag) = estimate(data, 1, methods[i], tmax)
  HIPresults[i,1] = hipf
  HIPresults[i,2:end] = hipx
  println(string(methods[i])," exited with ",flag)

  println("Fitting RIP process with ",string(methods[i]))
  (ripf, ripx, flag) = estimate(data, 0, methods[i], tmax)
  RIPresults[i, 1] = ripf
  RIPresults[i, 2:end] = ripx
  println(string(methods[i])," exited with",flag)
end

using DataFrames

results_HIP = DataFrame(Method = ["HIP_Paper"; methods], fmin = zeros(l+1),
  ρ = zeros(l+1), varɛ = zeros(l+1), varη = zeros(l+1), varα = zeros(l+1),
  varβ = zeros(l+1), corrαβ = zeros(l+1) )

results_RIP = DataFrame(Method = ["RIP_Paper"; methods], fmin = zeros(l+1),
  ρ = zeros(l+1), varɛ = zeros(l+1), varη = zeros(l+1), varα = zeros(l+1),
  varβ = zeros(l+1), corrαβ = zeros(l+1) )

paperrip = [0.988, 0.061, 0.015, 0.058]
for i in 1:4
  results_RIP[1, i+2] = paperrip[i]
  for m in 1:l
    results_RIP[1+m,2] = RIPresults[m, 1]
    results_RIP[1+m, i+2] = RIPresults[m, i+1]
  end
end

paperhip = [0.821, 0.047, 0.029, 0.022, 0.00038, -0.23]
for i in 1:6
  results_HIP[1, i+2] = paperhip[i]
  for m in 1:l
    results_HIP[1+m,2] = HIPresults[m, 1]
    results_HIP[1+m, i+2] = HIPresults[m, i+1]
  end
end

println(results_RIP)
println(results_HIP)
writetable("C:/Users/tew207/Desktop/BHPS_netinc_RIP.csv", results_RIP)
writetable("C:/Users/tew207/Desktop/BHPS_netincf_HIP.csv", results_HIP)
