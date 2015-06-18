using HDF5

CovN = h5read("C:/Users/tew207/My Documents/GitHub/psidJulia/output.h5", "Covariances")
Num = h5read("C:/Users/tew207/My Documents/GitHub/psidJulia/output.h5", "Observations")

# Invert each cohort's matrix, as HDF5 stores in row-major order
for i = 1:size(CovN, 3)
  CovN[:, :, i] = CovN[:, :, i]'
  Num[:, :, i] = Num[:, :, i]'
end

function estimARMA(CovN::Array{Float64, 3}, Num::Array{Int64, 3}, PH::Int64)
  # INPUTS: COVN (COVARIANCE MATRIX: AGE X AGE X TIME)
  #         NUM (# OF OBS USED TO COMPUTE EACH COV, SAME SIZE AS COVN)
  # PH=1 ESTIMATES THE HIP PROCESS
  # PH=0 ESTIMATES THE RIP PROCESS

  # Parameters
  lastcoh=1977; agecell=4; tmax=1992; Ageyes=0; choW=1; Ve=1; Vn=1;

  V = CovN
  Ivarcov = Num .> 10
  nlag = 29
  cmax = size(V,3)
  tik = lastcoh - 1966  # NOTE: this was "Ylastcoh" in original code
  hmax = size(V,1) + int(agecell/2)
  AA = Ageyes

  # start with some initial guess X0
  mypi1 = linspace(1,1.2,13)
  mypi1 = [mypi1,linspace(1.2,1.7,3)]
  if tmax > 20
    mypi1 = [mypi1,linspace(1.7,1.7,tmax-17)]
  end
  myeps1 = linspace(1,1,tmax-1)

  X0 = zeros(7 + tmax-2 + length(myeps1))

  if Ve==1 & Vn==1
      X0[7+2*AA:7+2*AA+tmax-2] = mypi1
      X0[7+2*AA+tmax-1:7+2*AA+2*tmax-3] = myeps1
  elseif Ve==1 & Vn==0
      X0[7+2*AA:7+2*AA+tmax-3] = myeps1
  elseif Ve==0 & Vn==1
      X0[7+2*AA:7+2*AA+tmax-2] = mypi1
  end

  # for m=1:5
  X0[2]=0.04; X0[3]=0.02; X0[5]=0.0005; X0[6]=-0.5

  for k=1:2, i=1:2, j=1:2
    X0[1] = 0.95 - (j-1)*0.20 # ρ (result: 0.821)
    X0[2] = 0.05 - (i-1)*0.02 # σ²(ɛ) (0.047)
    X0[3] = 0.02 + (i-1)*0.02 # σ²(η) (0.029)
    X0[4] = 0.05 + (k-1)*0.10 # σ²(α) (result 0.022)
    X0[5] = 0.0002          # σ²(β) (0.00038)
    X0[6] = -0.5            # corr(αβ) (-0.23)
    X0[1:6] = [0.80, 0.04, 0.02, 0.02, 0.00025, -0.23]

    (params,Fval,exitflag) =
      fminsearch(@MinDistOBJ_NEW_April_04,X0,options,V,Ivarcov,Num,tmax,
                 nlag,agecell,cmax,tik, Ageyes,choW,Ve,Vn,PH)

    #================= STAGE 2====================================#

    if Ageyes==1
      if Ve==1 & Vn==1
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,-1,-1,
            0.3*ones(1,tmax-1)',0.3*ones(1,tmax)']
        UB=[1.2, 2, 2, 0.5, 0.5, 1, 1, 1, 2.5*ones(1,tmax)', 4*ones(1,tmax)']
      elseif Ve==1 & Vn==0
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,-1,-1,0.3*ones(1,tmax)']
        UB=[1.2, 2,2,0.5,0.5,1,1,1,4*ones(1,tmax)']
      elseif Ve==0 & Vn==1
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,-1,-1,0.3*ones(1,tmax)']
        UB=[1.2, 2,2,0.5,0.5,1,1,1,4*ones(1,tmax)']
      else
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,-1,-1]
        UB=[1.2, 2,2,0.5,0.5,1,1,1]
      end
    else
      if Ve==1 & Vn==1
        LB=[-0.4,0.005,0.005,0.000001,0.000001,-1,
            0.3*ones(1,tmax)',0.3*ones(1,tmax)']
        UB=[1.2, 2,2,0.5,0.5,1,4*ones(1,tmax)',4*ones(1,tmax)']
      elseif Ve==1 & Vn==0
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,0.3*ones(1,tmax)']
        UB=[1.2, 2,2,0.5,0.5,1,4*ones(1,tmax)']
      elseif Ve==0 & Vn==1
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,0.3*ones(1,tmax)']
        UB=[1.2, 2,2,0.5,0.5,1,4*ones(1,tmax)']
      else
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1]
        UB=[1.2, 2,2,0.5,0.5,1]
      end
    end

    X0 = params
    (params,Fval,exitflag) =
      optimize(MinDistOBJ_DV_NEW_April_04,X0,LB,UB,
              V,Ivarcov,Num,tmax,nlag,agecell, cmax, tik, Ageyes,choW,Ve,Vn,PH)
    # @... is the function to minimize, X0 the initial condition,
    # [] empty arguments and LB and UB are lower and upper bounds
    parTEMP[:,i,j,k] = params
    FvalTEMP[i,j,k] = Fval
    flagTEMP[i,j,k] = exitflag
    X0 = params
  end
end

################################################################################

function MinDistOBJ_DV_NEW_April_04(params::Vector{Float64}, V::Array{Float64,3},
  Ivarcov, Num::Array{Float64,3}, tmax::Int64, nlag::Int64, agecell::Int64,
  cmax::Int64, tik::Int64, Ageyes::Int64, choW::Int64, Ve::Int64, Vn::Int64,
  PH::Int64)
  #V: Input, age x age matrix of empirical var-cov matrix
  #varcov is a age x age x time dimensional matrix. It has the variance
  #covariance matrix for individuals for each age for a given year, where age
  #corresponds to age-cells. last dimension controls for year and weare
  #planning to aggregare over that dimension in the estimation (not sure
  #though)
  #Ivarcov = let this be an indicator function the same size as varcov
  #indicating if we have observations in a certain age-cell time combination
  #and hence whether that cell contributes to the varcov matrix.

  AA = Ageyes

  if choW==0
    Myweight = Ivarcov
  elseif choW==1
    Myweight = Num
  end

  # Estimate with or without profile heterogeneity
  PH==1 ? qq=1 : qq=0

  # maximum experience: T + 2
  hmax = int(size(V,1)+(agecell/2))

  # varcov = zeros(T,T,size(V,3))
  varcov = similar(V)

  ρ = params[1]; varϵ = params[2]; varη = params[3]
  varα = params[4]; varβ = params[5]; corrαβ = params[6]

  if varα*varβ != 0
    covαβ = corrαβ*sqrt(varα*varβ)
  else
    covαβ = 0
  end

  hvec = Array(Float64, (3, hmax))
  myfi = Array(Float64, hmax)

  if Ageyes==1
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

  # for every t from 1 to 30, we are going to construct the A x A matrix
  # varcov. Note that for some combination of t and h, the age-cell may
  # contain no observations. be careful about that.

  function theoretical_varcov(varα::Float64, varβ::Float64, covαβ::Float64,
                          varη::Float64, varϵ::Float64, mypi::Array{Float64,1},
                          myeps::Array{Float64,1})
    # Calculate Variances
    ip_part_var = zeros(hmax, tmax)
    varz = similar(ip_part_var)
    vary = similar(ip_part_var)
    covary = Array(Float64, (hmax, nlag, tmax))

    for t = 1:tmax
      varu[1, t] = mypi[t]*varη
      for h = 3:hmax
        ip_part_var[h, t] = varα + 2*covαβ*h + varβ*h^2
        varz[h, 1] = [[mypi[1]*varη for j = 0:(h-1)]'*[ρ^(2*j) for j = 0:(h-1)]][1]
        if (t > 1) & (h>1)
          varz[h, t] = ρ^2*varz[h-1,t-1] + mypi[t]*varη
        end
        vary[h,t] = ip_part_var[h,t] + varz[h,t] + myeps[t]*varϵ
      end
    end

    # Calculate Autocovariances
    for h = int((agecell/2)+1):hmax, t = 1:tmax, n = 1:min(h-agecell/2+1, t-1, nlag)
      covary[h, n, t] = varα + 2*covαβ*(h+n) + varβ*h*(h+n)
                        + (ρ^n)*(varz[h-n, t-n])
    end
    covary

    varcov = Array(Float64, size(V))

    for coh = 1:cmax, h1 = 1:hmax-(agecell/2), h2 = max(1,h1-nlag+1):h1
      if ((h1+tik-coh >= 1) & ((h1+tik-coh) <= cmax))
        if (h1!=h2)
          varcov[h1,h2,coh] = covary[h1, h1-h2, h1+tik-coh]
        else
          varcov[h1,h2,coh] = vary[h1, h1+tik-coh]
        end
      else
        Ivarcov[h1,h2,coh]==0
      end
    end
    return varcov
  end

  countC = Array(Float64, (Tmax,Tmax))
  countD = similar(countC)
  # covTime holds the weighted theoretical covariances
  covTime = similar(countC)
  # CT holds the weighted empirical covariances (from V/CovN)
  CT = similar(countC)

  for t1=1:Tmax, t2=1:t1, coh=1:cmax
    if (t2-tik+coh>=1) & (t1-tik+coh<=hmax-(agecell/2+1))
      if (Myweight[t1-tik+coh,t2-tik+coh,coh] > 0)
        countC[t1,t2] = countC[t1,t2] + Myweight[t1-tik+coh,t2-tik+coh,coh]
        countD[t1,t2] = countD[t1,t2] + Myweight[t1-tik+coh,t2-tik+coh,coh]
        covTime[t1,t2] = covTime[t1,t2]
                        + Myweight[t1-tik+coh,t2-tik+coh,coh].*varcov[t1-tik+coh,t2-tik+coh,coh]
        CT[t1,t2] = CT[t1,t2] + Myweight[t1-tik+coh,t2-tik+coh,coh].*V[t1-tik+coh,t2-tik+coh,coh]
      end
    end

    if countC[t1,t2] > 0
      covTime[t1,t2] = covTime[t1,t2]/countD[t1,t2]
      CT[t1,t2] = CT[t1,t2]/countC[t1,t2]
    else
      covTime[t1,t2] = 0
      CT[t1,t2] = 0
    end
  end

  penalty1 = 10000000*(min(0.0, vare-0.0005))^2
  penalty2 = 100000*(min(0.0,varn-0.01))^2
  penalty8 = 100000*(max(0,ro-1.2))^2
  penalty3 = 100000*((min(0,ro+0.4))^2)

  if (max(mypi)*max(myfi)>=5)
    penalty9 = 10000000*(max(mypi)*max(myfi)-3)^2
  else
    penalty9 = 0
  end

  if (min(mypi)<0.5)
    penalty10 = 100000*(((min(mypi))-0.5)^2)
  else
    penalty10 = 0
  end

  if (min(myfi)<0.05)
    penalty11 = 100000*(((min(mypi))-0.05)^2)
  else
    penalty11 = 0
  end

  penaltycorr1 = 10000000*(min(corr01+1,0)^2)
  penaltycorr2 = 10000000*(max(corr01-1,0)^2)

  penaltyVAR = 10000000*(min(0.0,var1)^2)

  penalty = penalty1 + penalty2 + penalty3 + penalty8 + penalty9
           + penalty10 + penalty11 + penaltyVAR + penaltycorr1 + penaltycorr2

  myobjV = reshape(CT-covTime,Tmax*Tmax,1)
  size(myobjV)
  myobj = (myobjV'*myobjV)*(1+penalty)
end


beginning = 0
for i = 11:33
  if ~isnan(CovN[1,1,i])
    for j = 1:41
      if isnan(CovN[j,j,i])
        @printf "CovN[:,:,%d] from %d to %d\n" i 1 j-1
        break
      end
    end
  else
    for j = 1:41
      if ~isnan(CovN[j,j,i])
        beginning = j
        break
      end
    end
    if ~isnan(CovN[end,end,i])
      @printf "CovN[:,:,%d] from %d to %d\n" i beginning 41
    else
      for k = beginning:41
        if isnan(CovN[k,k,i])
          @printf "CovN[:,:,%d] from %d to %d\n" i beginning k-1
          break
        end
      end
    end
  end
end
