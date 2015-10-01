################################################################################
#####             PROFILE HETEROGENEITY AND INCOME RISK                    #####
################################################################################
nprocs()==CPU_CORES || addprocs(CPU_CORES-1)
import HDF5, NLopt
using DataFrames
@everywhere using HDF5

@everywhere begin
  tinit = 1968; tlast = 2013; tmax = tlast - tinit + 1
  data = "C:/Users/tew207/My Documents/GitHub/psidJulia/PSID_1968_2013.h5"
  CovEmp = h5read(data, "Covariances")
  Num = h5read(data, "Observations")
  cmax = size(CovEmp, 3)
  agemax = size(CovEmp, 2)

  # Invert matrices as hd5 stores in column-order
  for i = 1:cmax
    CovEmp[:, :, i] = CovEmp[:, :, i]'
    Num[:, :, i] = Num[:, :, i]'
  end

  # Parameters
  obs_indicator = convert(Array{Int64,3}, Num .> 10)
  hmax = agemax + 2

  # Determine no of lags by maximum distance between nonzero entry in columns
  nlag = 0
  for j = 1:cmax, i = agemax:-1:1
    lagnow = findlast(CovEmp[:,i,j]) - findfirst(CovEmp[:,i,j]) + 1
    (lagnow < nlag) || (nlag = lagnow)
  end

  # Determine new cohorts based on first 0 entry in CovEmp[1,1,:]
  tik = 0
  for i = 1:size(CovEmp,3)
    (abs(CovEmp[1,1,i]) > 1e-7) || (tik = i-1; break)
  end

  println("CovEmp has dimensions $(size(CovEmp)), with $nlag lags and tik $tik")
  include("C:/Users/tew207/My Documents/GitHub/psidJulia/theoretical_varcov.jl")
  include("C:/Users/tew207/My Documents/GitHub/psidJulia/min_dist_estimator.jl")
end

methods = [:GN_CRS2_LM, :GN_MLSL, :GN_ISRES, :GN_ESCH]
l = length(methods)
HIPresults = SharedArray(Float64, (2*l, 6+2*tmax+1))
RIPresults = SharedArray(Float64, (2*l, 6+2*tmax+1))

@sync @parallel for i = 1:length(methods)
  println("Fitting HIP process with ",string(methods[i]))
  (hipf, hipx, flag) = mdest(data, 1, methods[i], tmax)
  HIPresults[2*i-1,1] = hipf
  HIPresults[2*i-1,2:end] = hipx
  println(string(methods[i])," exited with ",flag)

  prt = zeros(12, 6)
  for p = 1:6
    prt[:,p]=hipx[p]; prt[2*p-1,p]+=0.01*hipx[p]; prt[2*p,p]-=0.01*hipx[p]
  end
  prt = hcat(prt, repmat(hipx[7:end]', 12, 1))

  f = obj_jac(hipx, 1); Ω = f*f'
  J = zeros(length(f),6)
  for p = 1:6
    J[:,p] = 0.5*(obj_jac(vec(prt[p,:]), 1) - f)/(0.01*hipx[p])
             + 0.5*(obj_jac(vec(prt[p+1,:]), 1) - f)/(-0.01*hipx[p])
  end
  Σ = inv(J'*J)*J'*Ω*J*inv(J'*J)
  HIPresults[2*i,2:7] = sqrt(abs(diag(Σ)))

  println("Fitting RIP process with ",string(methods[i]))
  (ripf, ripx, flag) = mdest(data, 0, methods[i], tmax)
  RIPresults[2*i-1, 1] = ripf
  RIPresults[2*i-1, 2:end] = ripx
  println(string(methods[i])," exited with",flag)

  prt = zeros(12, 6)
  for p = 1:6
    prt[:,p]=ripx[p]; prt[2*p-1,p]+=0.01*ripx[p]; prt[2*p,p]-=0.01*ripx[p]
  end
  prt = hcat(prt, repmat(ripx[7:end]', 12, 1))

  f = obj_jac(ripx, 0); Ω = f*f'
  J = zeros(length(f),6)
  for p = 1:6
    J[:,p] = 0.5*(obj_jac(vec(prt[p,:]), 0) - f)/(0.01*ripx[p])
             + 0.5*(obj_jac(vec(prt[p+1,:]), 0) - f)/(-0.01*ripx[p])
  end
  Σ = inv(J'*J)*J'*Ω*J*inv(J'*J)
  RIPresults[2*i, 2:7] = sqrt(abs(diag(Σ)))
end

firstcol = Array(AbstractString, (2*l+1))
firstcol[1,1] = "HIP_Paper"
for i = 2:2:2*l
  firstcol[i,1] = string(methods[Int(i/2)])
  firstcol[i+1,1] = ""
end

results_HIP = DataFrame(Method = firstcol, fmin = zeros(2*l+1),
  ρ = zeros(2*l+1), varɛ = zeros(2*l+1), varη = zeros(2*l+1), varα = zeros(2*l+1),
  varβ = zeros(2*l+1), corrαβ = zeros(2*l+1) )

results_RIP = DataFrame(Method = firstcol, fmin = zeros(2*l+1),
  ρ = zeros(2*l+1), varɛ = zeros(2*l+1), varη = zeros(2*l+1), varα = zeros(2*l+1),
  varβ = zeros(2*l+1), corrαβ = zeros(2*l+1) )
results_RIP[1,:Method] = "RIP_Paper"

paperrip = [0.988, 0.061, 0.015, 0.058, 0.0, 0.0]
paperhip = [0.821, 0.047, 0.029, 0.022, 0.00038, -0.23]
for i in 1:6
  results_HIP[1, i+2] = paperhip[i]
  results_RIP[1, i+2] = paperrip[i]
  for m in 1:l
    results_HIP[2*m,2] = HIPresults[2*m-1, 1]
    results_HIP[2*m, i+2] = HIPresults[2*m-1, i+1]
    results_HIP[2*m+1, i+2] = HIPresults[2*m, i+1]
    results_RIP[2*m,2] = RIPresults[2*m-1, 1]
    results_RIP[2*m, i+2] = RIPresults[2*m-1, i+1]
    results_RIP[2*m+1, i+2] = RIPresults[2*m, i+1]
  end
end

println(results_RIP)
println(results_HIP)
writetable("C:/Users/tew207/Desktop/PSID_1968_1986_RIP.csv", results_RIP)
writetable("C:/Users/tew207/Desktop/PSID_1968_1986_RIP_HIP.csv", results_HIP)
