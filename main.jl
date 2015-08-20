nprocs()==CPU_CORES || addprocs(CPU_CORES-1)

@everywhere begin
  using HDF5, Optim, NLopt, DataFrames

  path = "C:/Users/tew207/My Documents/GitHub/psidJulia/"

  CovEmp = h5read(path*"output_1981_1991.h5", "Covariances")
  Num = h5read(path*"output_1981_1991.h5", "Observations")
end

# Invert each cohort's matrix, as HDF5 stores in row-major order
@everywhere for i = 1:size(CovEmp, 3)
  CovEmp[:, :, i] = CovEmp[:, :, i]'
  Num[:, :, i] = Num[:, :, i]'
end

@everywhere begin
tinit = 1980
tlast = 1991
nlag = 10
minyrs = 8
agecell = 4
lastcoh = tlast - minyrs + 1
tmax = tlast - tinit + 1
tik = tlast - tinit - minyrs + int(agecell/2)
end

methods = [:GN_CRS2_LM, :GN_MLSL, :GN_ISRES, :GN_ESCH]
l = length(methods)

HIPresults = SharedArray(Float64, (l, 6+2*tmax+1))
RIPresults = SharedArray(Float64, (l, 6+2*tmax-1))

@sync @parallel for i = 1:length(methods)
  println("Now calculating ",string(methods[i]), ", HIP")
  (hipf, hipx, flag) = estimARMA(CovEmp, Num, 1, methods[i])
  HIPresults[i,1] = hipf
  HIPresults[i,2:end] = hipx
end

@sync @parallel for i = 1:length(methods)
  println("Now calculating ",string(methods[i]), ", RIP")
  (ripf, ripx, flag) = estimARMA(CovEmp, Num, 0, methods[i])
  RIPresults[i, 1] = ripf
  RIPresults[i, 2:end] = ripx
end

results_HIP = DataFrame(Method = ["HIP_Paper", methods],
                    fmin = zeros(l+1),
                    ρ = zeros(l+1), varɛ = zeros(l+1),
                    varη = zeros(l+1), varα = zeros(l+1),
                    varβ = zeros(l+1), corrαβ = zeros(l+1) )

paperhip = [0.821, 0.047, 0.029, 0.022, 0.00038, -0.23]
for i in 1:6
  results_HIP[1, i+2] = paperhip[i]
end

results_RIP = DataFrame(Method = ["RIP_Paper", methods],
                    fmin = zeros(l+1),
                    ρ = zeros(l+1), varɛ = zeros(l+1),
                    varη = zeros(l+1), varα = zeros(l+1),
                    varβ = zeros(l+1), corrαβ = zeros(l+1) )

paperrip = [0.988, 0.061, 0.015, 0.058, 0.0, 0.0]
for i in 1:6
  results_RIP[1, i+2] = paperrip[i]
end

HIPmean = mean(sdata(HIPresults), 1)
RIPmean = mean(sdata(RIPresults), 1)
