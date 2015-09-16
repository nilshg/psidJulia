using NLopt

function estimate(CovEmp::Array{Float64, 3}, Num::Array{Float64, 3},
  hip::Int64, method::Symbol, tmax=tmax)

  # Parameters
  obs_indicator = convert(Array{Int64,3}, Num .> 10)
  cmax = size(CovEmp,3)
  agemax = size(CovEmp,1); hmax = agemax + 2

  # Determine number of lags by maximum number of entries in CovEmp columns
  nlag = 0
  for j = 1:cmax, i = 1:agemax
    lagnow = sum(CovEmp[:,i,j].>1e-6)
    (lagnow < nlag) || (nlag = lagnow)
  end

  # Determine new cohorts based on first 0 entry in CovEmp[1,1,:]
  tik = 0
  for i = 1:size(CovEmp,3)
    (abs(CovEmp[1,1,i]) > 1e-7) || (tik = i-1; break)
  end

  # start with some initial guess x_0
  π_0 = [linspace(1, 1.2, 13); linspace(1.2, 1.7, 3)]
  if tmax > 16
    π_0 = [π_0; linspace(1.7, 1.7, tmax-16)]
  end

  ϕ_0 = linspace(1,1,tmax)

  # Initial conditions
  x_0 = zeros(6 + 2*tmax)
  x_0[1:6] = [0.8, 0.04, 0.02, 0.02, hip*0.0003, -0.23*hip]
  x_0[7:7+tmax-1] = π_0[1:tmax]
  x_0[7+tmax:end] = ϕ_0
  
  ##############################################################################
  ## Objective, depending on observed and theoretical covariance structure

  function objective(params::Vector{Float64}, grad, CovEmp=CovEmp,
    obs_indicator=obs_indicator, weight=Num, tmax=tmax, nlag=nlag, cmax=cmax,
    tik=tik, hip=hip)

    ρ = params[1]; varϵ = params[2]; varη = params[3]
    varα = params[4]; varβ = params[5]; corrαβ = params[6]

    if (varα >= 0.0) && (varβ >= 0.0)
      covαβ = corrαβ*sqrt(abs(varα*varβ))
    else
      covαβ = 0.0
    end

    π = zeros(tmax)
    ϕ = zeros(tmax)

    π[1] = 1
    π[2:tmax] = params[7:7+tmax-2]
    ϕ[1] = 1
    ϕ[2:tmax-1] = params[7+tmax:7+2*tmax-3]
    ϕ[tmax] = 1

    ϕ .^= 2
    π .^= 2

    # Construct theoretical var-cov matrix (and assert that it has entries in
    # the same places as the empirical one)
    varcov = tvg(varα,varβ,covαβ,varϵ,varη,ρ,ϕ,π,hmax,tmax,nlag,tik,hip)

    if (sum(abs(CovEmp).>1e-8).==sum(abs(varcov).>1e-8)!=length(CovEmp))
      println("WARNING")
      println("Entries in CovEmp: $(sum((abs(CovEmp).>1e-6)))")
      println("Entries in varcov: $(sum((abs(varcov).>1e-6)))")
      for i = 1:size(CovEmp,1), j = 1:size(CovEmp,2), c = 1:size(CovEmp,3)
        if (abs(CovEmp[i,j,c])>1e-6) & (abs(varcov[i,j,c])<=1e-6)
          println("Missing value in varcov[$((i,j,c))]")
          println("CovEmp[$((i,j,c))]=$(CovEmp[i,j,c]), varcov[$((i,j,c))]=$(varcov[i,j,c])")
        end
        if (abs(CovEmp[i,j,c])<1e-6) & (abs(varcov[i,j,c])>=1e-6)
          println("Extra value in varcov[$((i,j,c))]")
          println("CovEmp[$((i,j,c))]=$(CovEmp[i,j,c]), varcov[$((i,j,c))]=$(varcov[i,j,c])")
        end
      end
    end
    @assert sum((abs(CovEmp).>1e-8).==(abs(varcov).>1e-8))==length(CovEmp)

    # matrices to hold observations for
    ∑n = zeros(tmax,tmax)
    # weighted theoretical covariances (from varcov)
    theor = zeros(tmax,tmax)
    # weighted empirical covariances (from CovEmp)
    empir = zeros(tmax,tmax)

    for t1=1:tmax, t2=1:t1
      for coh=1:cmax
        if (t2-tik+coh > 0) && (t1-tik+coh+1 <= agemax-2)
          n = weight[t1-tik+coh, t2-tik+coh, coh]
          ∑n[t1, t2] += n
          theor[t1,t2] += n.*varcov[t1-tik+coh,t2-tik+coh,coh]
          empir[t1,t2] += n.*CovEmp[t1-tik+coh,t2-tik+coh,coh]
        end
      end

      if ∑n[t1,t2] > 0
        theor[t1,t2] = theor[t1,t2]/∑n[t1,t2]
        empir[t1,t2] = empir[t1,t2]/∑n[t1,t2]
      else
        theor[t1,t2] = 0
        empir[t1,t2] = 0
      end
    end

    penalty1 = 10000000.0*(min(0.0, varϵ-0.0005))^2
    penalty2 = 100000.0*(min(0.0, varη - 0.01))^2
    penalty8 = 100000.0*(max(0.0, ρ - 1.2))^2
    penalty3 = 100000.0*(min(0.0, ρ + 0.4))^2

    if (maximum(π) >= 5.0)
        penalty9 = 10000000.0*(maximum(π)-3.0)^2.0
    else
        penalty9 = 0.0
    end

    if (minimum(π) < 0.5)
        penalty10 = 100000.0*(((minimum(π))-0.5)^2.0)
    else
        penalty10 = 0.0
    end

    if (minimum(π) < 0.05)
        penalty11 = 100000.0*(((minimum(π))-0.05)^2.0)
    else
        penalty11 = 0.0
    end

    penaltycorr1 = 10000000.0*(min(corrαβ+1, 0.0)^2.0)
    penaltycorr2 = 10000000.0*(max(corrαβ-1, 0.0)^2.0)

    penaltyVAR = 10000000.0*(min(0.0, varβ)^2.0)

    penaltyβ = 100000.0*((hip-1)*(varβ))^2.0
    penaltycorr3 = 100000.0*((hip-1)*(abs(corrαβ)))^2.0

    penalty = (penalty1 + penalty2 + penalty3 + penalty8
             + penalty9 + penalty10 + penalty11 + penaltyVAR
             + penaltycorr1 + penaltycorr2 + penaltyβ
             + penaltycorr3)

    temp = collect(empir-theor)
    obj = (temp'*temp)[1]*(1+penalty)
  end

  lb = [[ 0.4; 5e-4; 5e-4; 1e-6; 1e-6; -1.]; 0.3*ones(2*(tmax))]
  ub = [[ 1.2;  2.0;  2.0;  0.5;  0.5;  1.];   4*ones(2*(tmax))]

  # Global minimization
  opt_1 = Opt(method, length(x_0))
  if method == :GN_MLSL
    localopt = Opt(:LN_SBPLX, length(x_0))
    local_optimizer!(opt_1, localopt)
  end
  lower_bounds!(opt_1, lb)
  upper_bounds!(opt_1, ub)
  min_objective!(opt_1, objective)
  ftol_abs!(opt_1, 1e-12)
  maxeval!(opt_1, 50000)
  maxtime!(opt_1, 1000)
  (optf_1, optx_1, flag_1) = optimize(opt_1, x_0)

  # Local minimization
  opt_2 = Opt(:LN_SBPLX, length(x_0))
  lower_bounds!(opt_2, lb)
  upper_bounds!(opt_2, ub)
  min_objective!(opt_2, objective)
  ftol_abs!(opt_2, 1e-12)
  maxeval!(opt_2, 50000)
  maxtime!(opt_2, 1000)
  (optf, optx, flag) = optimize(opt_2, optx_1)
end
