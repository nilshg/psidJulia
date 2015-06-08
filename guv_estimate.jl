CovN = randn(26,26,12)
Num = round(20*randn(26,26,12))


function estimARMA(CovN::Array{Float64, 3}, Num::Array{Int64, 3}, PH::Int64)
  # INPUTS: COVN (COVARIANCE MATRIX: AGE X AGE X TIME)
  #         NUM (# OF OBS USED TO COMPUTE EACH COV, SAME SIZE AS COVN)
  # PH=1 ESTIMATES THE HIP PROCESS
  # PH=0 ESTIMATES THE RIP PROCESS

  # Parameters
  lastcoh=1977; agecell=4; Tmax=1992; Ageyes=0; choW=1; Ve=1; Vn=1;

  V = CovN
  Ivarcov = Num .> 10
  nlag = 25
  Cmax = size(V,3)
  tik = lastcoh - 1966  # NOTE: this was "Ylastcoh" in original code
  hmax = size(V,1) + int(agecell/2)
  AA = Ageyes

  # start with some initial guess X0
  mypi1 = linspace(1,1.2,13)
  mypi1 = [mypi1,linspace(1.2,1.7,3)]
  if Tmax > 20
    mypi1 = [mypi1,linspace(1.7,1.7,Tmax-17)]
  end
  myeps1 = linspace(1,1,Tmax-1)

  X0 = zeros(7 + Tmax-2 + length(myeps1))

  if Ve==1 & Vn==1
      X0[7+2*AA:7+2*AA+Tmax-2] = mypi1
      X0[7+2*AA+Tmax-1:7+2*AA+2*Tmax-3] = myeps1
  elseif Ve==1 & Vn==0
      X0[7+2*AA:7+2*AA+Tmax-3] = myeps1
  elseif Ve==0 & Vn==1
      X0[7+2*AA:7+2*AA+Tmax-2] = mypi1
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
      fminsearch(@MinDistOBJ_NEW_April_04,X0,options,V,Ivarcov,Num,Tmax,
                 nlag,agecell,Cmax,tik, Ageyes,choW,Ve,Vn,PH)

    #================= STAGE 2====================================#

    if Ageyes==1
      if Ve==1 & Vn==1
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,-1,-1,
            0.3*ones(1,Tmax-1)',0.3*ones(1,Tmax)']
        UB=[1.2, 2, 2, 0.5, 0.5, 1, 1, 1, 2.5*ones(1,Tmax)', 4*ones(1,Tmax)']
      elseif Ve==1 & Vn==0
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,-1,-1,0.3*ones(1,Tmax)']
        UB=[1.2, 2,2,0.5,0.5,1,1,1,4*ones(1,Tmax)']
      elseif Ve==0 & Vn==1
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,-1,-1,0.3*ones(1,Tmax)']
        UB=[1.2, 2,2,0.5,0.5,1,1,1,4*ones(1,Tmax)']
      else
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,-1,-1]
        UB=[1.2, 2,2,0.5,0.5,1,1,1]
      end
    else
      if Ve==1 & Vn==1
        LB=[-0.4,0.005,0.005,0.000001,0.000001,-1,
            0.3*ones(1,Tmax)',0.3*ones(1,Tmax)']
        UB=[1.2, 2,2,0.5,0.5,1,4*ones(1,Tmax)',4*ones(1,Tmax)']
      elseif Ve==1 & Vn==0
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,0.3*ones(1,Tmax)']
        UB=[1.2, 2,2,0.5,0.5,1,4*ones(1,Tmax)']
      elseif Ve==0 & Vn==1
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1,0.3*ones(1,Tmax)']
        UB=[1.2, 2,2,0.5,0.5,1,4*ones(1,Tmax)']
      else
        LB=[-0.4,0.0005,0.0005,0.000001,0.000001,-1]
        UB=[1.2, 2,2,0.5,0.5,1]
      end
    end

    X0 = params
    (params,Fval,exitflag) =
      fmincon(@MinDistOBJ_DV_NEW_April_04,X0,[],[],[],[],LB,UB,[],options,
              V,Ivarcov,Num,Tmax,nlag,agecell, Cmax, tik, Ageyes,choW,Ve,Vn,PH)
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
  Ivarcov, Num::Array{Float64,3}, Tmax::Int64, nlag::Int64, agecell::Int64,
  Cmax::Int64, tik::Int64, Ageyes::Int64, choW::Int64, Ve::Int64, Vn::Int64,
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
  else
    error("ERROR: choW incorrectly defined")
  end

  if PH==1 # WITH profile heterogeneity (HIP process)
    qq=1
  else     # WITHOUT profile heterogeneity (RIP)
    qq=0
  end

  # maximum experience: T + 2
  hmax = size(V,1)+(agecell/2)
  # varcov = zeros(T,T,size(V,3))
  varcov = Array(Float64, (hmax-int(agecell/2),hmax-int(agecell/2),Cmax))

  ro = params[1]; vare = params[2]; varn = params[3]
  var0 = params[4]; var1 = params[5]; corr01 = params[6]

  if (var0>=0) & (var1>=0)
    cov01 = corr01*sqrt(var0*var1)
  else
    cov01 = 0
  end

  if Ageyes==1
    fi1=params[7]
    fi2=params[8]
    fivec = [1,fi1,fi2]/1
    for h = 1:hmax
      hvec[:,h] = [1,(h-1),(h-1)^2]'
      myfi[h] = fivec*hvec[:,h]
    end
  else
    myfi=ones(1,hmax)
  end

  mypi = Array(Float64, 1998)
  myeps = similar(mypi)
  if Ve==1 & Vn==1
    mypi[2:Tmax] = params[7+2*AA:7+2*AA+Tmax-2]
    mypi[1] = 1
    myeps[1] = 1
    myeps[Tmax] = 1
    myeps[2:Tmax-1] = params[7+2*AA+Tmax-1:7+2*AA+2*Tmax-4]
  elseif Ve==1 & Vn==0
    mypi=ones(1,Tmax)
    myeps[1] = 1
    myeps[Tmax] = 1
    myeps[2:Tmax-1] = params[7+2*AA:7+2*AA+Tmax-3]
  elseif Ve==0 & Vn==1
    mypi[1] = 1
    mypi[2:Tmax] = params[7+2*AA:7+2*AA+Tmax-2]
    myeps = ones[1,Tmax]
  else
    mypi = ones[1,Tmax]
    myeps = ones[1,Tmax]
  end

  myeps = myeps.^2
  mypi = mypi.^2

  # for every t from 1 to 30, we are going to construct the A x A matrix
  # varcov. Note that for some combination of t and a, the age-cell may
  # contain no observations. be careful about that.

  # find var(y_{h,t})
  mmm = Array(Float64, (hmax, Tmax))
  varu = similar(mmm)
  for t = 1:Tmax, h = (agecell/2+1):hmax
    if h <= t
      k = h-1
      mmm[h,t] = var0 + qq*var1*h^2 +  qq*2*cov01*h
      varu[h,t] = mmm[h,t] + myeps[t]*vare
                  + varn*((((ro^2).^[k:-1:0]).*mypi[t-h+1:t])'*(myfi[1:h]))
    elseif h > t
      k = t-1
      mmm[h,t] = var0 + qq*var1*h^2 + qq*2*cov01*h
      varu[h,t]
        = mmm[h,t] + myeps[t]*vare
          + varn*((((ro^2).^[k:-1:0]).*mypi[t-h+1:t])'*(myfi[1:h]))
          + varn*(mypi[1])*(((ro^2).^[h-1:-1:t])*(myfi[1:h-t]'))
    end
  end

  # find cov(y_{h,t}, y_{h+n,t+n})
  nnn = Array(Float64, (hmax, Tmax, Cmax))
  covu = similar(nnn)
  for t = 1:Tmax, h = (agecell/2+1):hmax
    for n = 1:minimum([h-(agecell/2+1),t-1,nlag])
      # This is equation (4) in the paper
      nnn[h,n,t] = var0 + qq*var1*h*(h+n) + qq*cov01*(2h+n)
      covu[h,n,t]
        = nnn[h,n,t]
          + (ro^n)*(varu[h-n,t-n] - mmm[h-n,t-n] - myeps[t-n]*vare)
    end
  end

  for coh = 1:Cmax, h1 = 1:hmax-(agecell/2), h2 = max(1,h1-nlag+1):h1
    if ((h1+tik-coh >= 1) & ((h1+tik-coh) <= Tmax))
      if (h1!=h2)
        varcov[h1,h2,coh]=covu[h1+(agecell/2),h1-h2,h1+tik-coh]
      else
        varcov[h1,h2,coh]=varu[h1+(agecell/2),h1+tik-coh]
      end
    else
      Ivarcov[h1,h2,coh]==0
    end
  end

  countC = Array(Float64, (Tmax,Tmax))
  countD = similar(countC)
  # covTime holds the weighted theoretical covariances
  covTime = similar(countC)
  # CT holds the weighted empirical covariances (from V/CovN)
  CT = similar(countC)

  for t1=1:Tmax, t2=1:t1, coh=1:Cmax
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
