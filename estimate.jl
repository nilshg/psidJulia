using HDF5, Optim

path = "C:/Users/tew207/My Documents/GitHub/psidJulia/"

CovN = h5read(path*"output.h5", "Covariances")
Num = h5read(path*"output.h5", "Observations")

# Invert each cohort's matrix, as HDF5 stores in row-major order
for i = 1:size(CovN, 3)
  CovN[:, :, i] = CovN[:, :, i]'
  Num[:, :, i] = Num[:, :, i]'
end

function estimARMA(CovEmp::Array{Float64, 3}, Num::Array{Float64, 3}, hip::Int64)
  # INPUTS: CovEmp (empirical covariances, age*age*cohorts)
  #         Num (# of observations used to calculate covariances in CovEmp)
  # hip = 1 estimates HIP process
  # hip = 0 estimates RIP process

  # Parameters
  lastcoh = 1977; agecell = 4; tmax = 30; nlag = 29;
  AA = 0; choW = 1; Ve = 1; Vn = 1;
  obs_indicator = Num .> 10
  cmax = size(CovEmp,3)
  tik = lastcoh - 1966  # NOTE: this was "Ylastcoh" in original code
  hmax = size(CovEmp,1) + int(agecell/2)

  # start with some initial guess x_0
  mypi1 = linspace(1,1.2,13)
  mypi1 = [mypi1,linspace(1.2,1.7,3)]
  if tmax > 20
    mypi1 = [mypi1,linspace(1.7,1.7,tmax-17)]
  end
  myeps1 = linspace(1,1,tmax-1)

  x_0 = zeros(7 + tmax-2 + length(myeps1))

  if Ve==1 & Vn==1
      x_0[7+2*AA:7+2*AA+tmax-2] = mypi1
      x_0[7+2*AA+tmax-1:7+2*AA+2*tmax-3] = myeps1
  elseif Ve==1 & Vn==0
      x_0[7+2*AA:7+2*AA+tmax-3] = myeps1
  elseif Ve==0 & Vn==1
      x_0[7+2*AA:7+2*AA+tmax-2] = mypi1
  end

  # for m=1:5
  x_0[2]=0.04; x_0[3]=0.02; x_0[5]=0.0005; x_0[6]=-0.5
  x_0[1:6] = [0.80, 0.04, 0.02, 0.02, 0.00025, -0.23]

  ################################################################################
  ## Function to construct theoretical var-cov matrix given parameters of income
  ## process
  ################################################################################

  function theoretical_varcov(varα::Float64, varβ::Float64, covαβ::Float64,
    ρ::Float64, varη::Float64, varϵ::Float64, mypi::Array{Float64,1},
    myeps::Array{Float64,1}, hip::Int64)

    # Calculate CovEmpariances
    ip_part_var = zeros(hmax, tmax)
    varz = zeros(hmax, tmax)
    vary = zeros(hmax, tmax)
    covary = zeros(hmax, nlag, tmax)

    for t = 1:tmax
      varz[1, t] = mypi[t]*varη
      for h = 1:hmax-2
        ip_part_var[h, t] = varα + hip*2*covαβ*h + hip*varβ*h^2
        varz[h, 1] = [[mypi[1]*varη for j = 0:(h-1)]'*[ρ^(2*j) for j = 0:(h-1)]][1]
        if (t > 1) & (h>1)
          varz[h, t] = ρ^2*varz[h-1,t-1] + mypi[t]*varη
        end
        vary[h,t] = ip_part_var[h,t] + varz[h,t] + myeps[t]*varϵ
      end
    end

    # Calculate Autocovariances
    for h = 1:hmax, t = 1:tmax
      for n = 1:min(h-agecell/2+1, t-1, nlag)
        covary[h, n, t] = varα + hip*2*covαβ*(h+n) + hip*varβ*h*(h+n)
          + (ρ^n)*(varz[h, t])
      end
    end

    varcov = zeros(size(CovEmp))

    for coh = 1:cmax, h1 = 1:hmax-(agecell/2), h2 = max(1,h1-nlag+1):h1
      if ((h1+tik-coh >= 1) & ((h1+tik-coh) <= tmax))
        if (h1!=h2)
          varcov[h1,h2,coh] = covary[h1, h1-h2, h1+tik-coh]
        else
          varcov[h1,h2,coh] = vary[h1, h1+tik-coh]
        end
      else
        obs_indicator[h1,h2,coh] = 0
      end
    end

    return varcov
  end

  ################################################################################
  ## Objective function, depending on observed and theoretical covariance
  ## structure
  ################################################################################

  function objective(params::Vector{Float64}, CovEmp=CovEmp,
    obs_indicator=obs_indicator, Num=Num, tmax=tmax, nlag=nlag, agecell=agecell,
    cmax=cmax, tik=tik, AA=AA, choW=choW, Ve=Ve, Vn=Vn, hip=hip)

    # Not sure yet what this is doing
    choW==0 ? weight = obs_indicator : weight = Num

    # Estimate with or without profile heterogeneity
    hip==1 ? qq = 1 : qq = 0

    # maximum experience: T + 2 (WHY?)
    hmax = int(size(CovEmp,1)+(agecell/2))

    # varcov = zeros(T,T,size(CovEmp,3))
    varcov = zeros(size(CovEmp))

    ρ = params[1]; varϵ = params[2]; varη = params[3]
    varα = params[4]; varβ = params[5]; corrαβ = params[6]

    varα*varβ != 0 ? covαβ = corrαβ*sqrt(abs(varα*varβ)) : covαβ = 0

    hvec = Array(Float64, (3, hmax))
    myfi = Array(Float64, hmax)

    if AA==1
      fi1=params[7]
      fi2=params[8]
      fivec = [1.0,fi1,fi2]
      for h = 1:hmax
        hvec[:,h] = [1.0,(h-1),(h-1)^2]
        myfi[h] = [fivec'*hvec[:,h]][1]
      end
    else
      myfi = ones(1, hmax)
    end

    mypi = Array(Float64, tmax)
    myeps = similar(mypi)
    if Ve==1 & Vn==1
      mypi[2:tmax] = params[7+2*AA:7+2*AA+tmax-2]
      mypi[1] = 1
      myeps[1] = 1
      myeps[tmax] = 1
      myeps[2:tmax-1] = params[7+2*AA+tmax-1:7+2*AA+2*tmax-4]
    elseif Ve==1 & Vn==0
      mypi=ones(1,tmax)
      myeps[1] = 1
      myeps[tmax] = 1
      myeps[2:tmax-1] = params[7+2*AA:7+2*AA+tmax-3]
    elseif Ve==0 & Vn==1
      mypi[1] = 1
      mypi[2:tmax] = params[7+2*AA:7+2*AA+tmax-2]
      myeps = ones[1,tmax]
    else
      mypi = ones[1,tmax]
      myeps = ones[1,tmax]
    end

    myeps = myeps.^2
    mypi = mypi.^2

    # Construct theoretical var-cov matrix
    varcov =
      theoretical_varcov(varα, varβ, covαβ, ρ, varη, varϵ, mypi, myeps, hip)

    # matrices to hold observations for
    nobs = zeros(tmax,tmax)
    # weighted theoretical covariances (from varcov)
    theor = zeros(tmax,tmax)
    # weighted empirical covariances (from CovEmp)
    empir = zeros(tmax,tmax)

    for t1=1:tmax, t2=1:t1, coh=1:cmax
      if (t2-tik+coh>=1) & (t1-tik+coh<=hmax-(agecell/2+1))
        if (weight[t1-tik+coh,t2-tik+coh,coh] > 0)
          obs = weight[t1-tik+coh, t2-tik+coh, coh]
          nobs[t1, t2] += obs
          theor[t1,t2] += obs.*varcov[t1-tik+coh,t2-tik+coh,coh]
          empir[t1,t2] += obs.*CovEmp[t1-tik+coh,t2-tik+coh,coh]
        end
      end

      if nobs[t1,t2] >= 1
        theor[t1,t2] = theor[t1,t2]/nobs[t1,t2]
        empir[t1,t2] = empir[t1,t2]/nobs[t1,t2]
      else
        theor[t1,t2] = 0
        empir[t1,t2] = 0
      end
    end

    penalty1 = 10000000*(min(0.0, varϵ-0.0005))^2
    penalty2 = 100000*(min(0.0, varη-0.01))^2
    penalty8 = 100000*(max(0, ρ-1.2))^2
    penalty3 = 100000*((min(0, ρ+0.4))^2)

    if (maximum(mypi)*maximum(myfi)>=5)
      penalty9 = 10000000*(maximum(mypi)*maximum(myfi)-3)^2
    else
      penalty9 = 0
    end

    if (minimum(mypi)<0.5)
      penalty10 = 100000*(((minimum(mypi))-0.5)^2)
    else
      penalty10 = 0
    end

    if (minimum(myfi)<0.05)
      penalty11 = 100000*(((minimum(mypi))-0.05)^2)
    else
      penalty11 = 0
    end

    penaltycorr1 = 10000000*(min(corrαβ+1, 0)^2)
    penaltycorr2 = 10000000*(max(corrαβ-1, 0)^2)
    penaltyVAR = 10000000*(min(0.0, varβ)^2)

    penalty = penalty1 + penalty2 + penalty3 + penalty8 + penalty9
             + penalty10 + penalty11 + penaltyVAR + penaltycorr1 + penaltycorr2

    temp = [empir-theor][:]
    obj = [(temp'*temp)*(1+penalty)][1]
  end

  optimum = optimize(objective, x_0)
  return optimum.minimum
end

@time hip_estimates = estimARMA(CovN, Num, 1)
@time rip_estimates = estimARMA(CovN, Num, 0)
