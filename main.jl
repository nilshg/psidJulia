################################################################################
#####             PROFILE HETEROGENEITY AND INCOME RISK                    #####
################################################################################
#
# hip = 1 estimates a 1-by-(6+2T) vector of parameters
# [ρ, varϵ, varη, varα, varβ, covαβ, ϕ_1, ..., ϕ_T, π_1, ..., π_T]
#
# hip = 2 estimates a 1-by-(3+2T) vector of parameters
# [ρ, varϵ, varη, ϕ_1, ..., ϕ_tmax, π_1, ..., π_T]
#
################################################################################
nprocs()==CPU_CORES || addprocs(CPU_CORES-1)

@everywhere begin
using HDF5, Optim, NLopt, DataFrames

  path = "C:/Users/tew207/My Documents/GitHub/psidJulia/"
  include(path*"theoretical_varcov.jl")
  include(path*"objective.jl")

  CovEmp = h5read(path*"output_"*string(tinit)*"_"*string(tlast)*".h5", "Covariances")
  Num = h5read(path*"output_"*string(tinit)*"_"*string(tlast)*".h5", "Observations")

  tinit = 1968; tlast = 2013
  nlag = 10
end

# Invert each cohort's matrix, as HDF5 stores in row-major order
@everywhere for i = 1:size(CovEmp, 3)
  CovEmp[:, :, i] = CovEmp[:, :, i]'
  Num[:, :, i] = Num[:, :, i]'
end

@everywhere begin
  tmax = tlast - tinit + 1
  tik = 0
  for i = 1:size(CovEmp,3)
    if abs(CovEmp[1,1,i]) < 1e-7
      tik = i-1
      break
    end
  end
end

methods = [:GN_CRS2_LM, :GN_MLSL, :GN_ISRES, :GN_ESCH]
l = length(methods)
HIPresults = SharedArray(Float64, (l, 6+2*tmax+1))
RIPresults = SharedArray(Float64, (l, 6+2*tmax-1))

@sync @parallel for i = 1:length(methods)
  println("Now calculating ",string(methods[i]), ", HIP")
  (hipf, hipx, flag) = estimate(CovEmp, Num, 1, methods[i])
  HIPresults[i,1] = hipf
  HIPresults[i,2:end] = hipx

  println("Now calculating ",string(methods[i]), ", RIP")
  (ripf, ripx, flag) = estimate(CovEmp, Num, 0, methods[i])
  RIPresults[i, 1] = ripf
  RIPresults[i, 2:end] = ripx
end

results_HIP = DataFrame(Method = ["HIP_Paper", methods],
  fmin = zeros(l+1), ρ = zeros(l+1), varɛ = zeros(l+1),
  varη = zeros(l+1), varα = zeros(l+1),
  varβ = zeros(l+1), corrαβ = zeros(l+1) )

results_RIP = DataFrame(Method = ["RIP_Paper", methods],
  fmin = zeros(l+1), ρ = zeros(l+1), varɛ = zeros(l+1),
  varη = zeros(l+1), varα = zeros(l+1),
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

println(results_HIP)
println(results_RIP)
