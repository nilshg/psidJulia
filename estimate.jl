using HDF5, Optim, NLopt

path = "C:/Users/tew207/My Documents/GitHub/psidJulia/"

CovEmp = h5read(path*"output.h5", "Covariances")
Num = h5read(path*"output.h5", "Observations")

# Invert each cohort's matrix, as HDF5 stores in row-major order
for i = 1:size(CovEmp, 3)
  CovEmp[:, :, i] = CovEmp[:, :, i]'
  Num[:, :, i] = Num[:, :, i]'
end

function estimARMA(CovEmp::Array{Float64, 3}, Num::Array{Float64, 3}, hip::Int64,
                   use_NLopt::Int64)
  # INPUTS: CovEmp (empirical covariances, age*age*cohorts)
  #         Num (# of observations used to calculate covariances in CovEmp)
  # hip = 1 estimates HIP process
  # hip = 0 estimates RIP process

  # Parameters
  lastcoh = 1977; agecell = 4; tmax = 29; nlag = 29;
  obs_indicator = convert(Array{Int64,3}, Num .> 10)
  cmax = size(CovEmp,3)
  tik = lastcoh - 1966
  hmax = size(CovEmp,1)

  # start with some initial guess x_0
  π_1 = linspace(1, 1.2, 13)
  π_1 = [π_1, linspace(1.2, 1.7, 3)]
  if tmax > 20
    π_1 = [π_1, linspace(1.7, 1.7, tmax-16)]
  end

  ϕ_1 = linspace(1,1,tmax)

  x_0 = zeros(6 + 2*tmax)
  x_0[7:7+tmax-1] = π_1
  x_0[7+tmax:end] = ϕ_1

  x_0[2] = 0.04; x_0[3] = 0.02; x_0[5] = 0.0005; x_0[6] = -0.5
  x_0[1:6] = [0.80, 0.04, 0.02, 0.02, 0.00025, -0.23]

  ##############################################################################
  ## Function to construct theoretical var-cov matrix given parameters of income
  ## process
  ##############################################################################

  function theoretical_varcov(varα::Float64, varβ::Float64, covαβ::Float64,
    ρ::Float64, varη::Float64, varϵ::Float64, π::Array{Float64,1},
    ϕ::Array{Float64,1}, hmax::Int64, tmax::Int64, nlag::Int64, hip::Int64)

    # Calculate CovEmpariances
    ip_part_var = zeros(hmax); varz = zeros(hmax, tmax)
    vary = zeros(hmax, tmax); covary = zeros(hmax, nlag, tmax)

    for t = 1:tmax
      varz[1, t] = π[t]*varη
      for h = 1:hmax-2
        ip_part_var[h] = varα + hip*2*covαβ*h + hip*varβ*h^2
        varz[h, 1] = [[π[1]*varη for j = 0:(h-1)]'*[ρ^(2*j) for j = 0:(h-1)]][1]
        if (t > 1) && (h>1)
          varz[h, t] = ρ^2*varz[h-1,t-1] + π[t]*varη
        end
        vary[h,t] = ip_part_var[h] + varz[h,t] + ϕ[t]*varϵ
      end
    end

    # Calculate Autocovariances
    for h = 1:hmax, t = 1:tmax
      for n = 1:min(hmax-h, tmax-t, nlag)
        covary[h, n, t] = varα + hip*2*covαβ*(h+n) + hip*varβ*h*(h+n)
                          + (ρ^n)*(varz[h, t])
      end
    end

    varcov = zeros(size(CovEmp))

    for coh = 1:cmax, h1 = 1:hmax-(agecell/2), h2 = max(1,h1-nlag+1):h1
      if ((h1+tik-coh >= 1) && ((h1+tik-coh) <= tmax))
        if (h1!=h2)
          varcov[h1,h2,coh] = covary[h1, h1-h2, h1+tik-coh]
        else
          varcov[h1,h2,coh] = vary[h1, h1+tik-coh]
        end
      end
    end

    return varcov
  end

  ##############################################################################

  function theoretical_varcov_guv(varα::Float64, varβ::Float64, covαβ::Float64,
    varϵ::Float64, varη::Float64, ρ::Float64, ϕ::Array{Float64,1},
    π::Array{Float64,1}, hmax::Int64, tmax::Int64, nlag::Int64, hip::Int64)

    ip_part = zeros(hmax)
    vary = zeros(hmax, tmax)

    for t = 1:tmax, h = 3:hmax
      ip_part[h] = varα + hip*(varβ*h^2. + 2*covαβ*h)
      if h <= t
        vary[h, t] =
          (ip_part[h] + ϕ[t]*varϵ
           + sum(varη*(((ρ^2).^[(h-1):-1:0]).*π[t-h+1:t])) )
      else
        vary[h, t] =
          ( ip_part[h] + ϕ[t]*varϵ
           + sum(varη*((((ρ^2).^[(t-1):-1:0]).*π[1:t])))
           + sum(varη*π[1]*(((ρ^2).^[(h-1):-1:t]))) )
      end
    end

    covy = zeros(hmax, nlag, tmax)

    for t= 1:tmax, h = 3:hmax, n = 1:min(h-3,t-1,nlag)
      covy[h, n, t] =
        ( varα + hip*(varβ*h*(h-n) + covαβ*(2*h-n))
         + (ρ^n)*(vary[h-n,t-n]-ip_part[h-n] - ϕ[t-n]*varϵ) )
    end

    varcov = zeros(hmax-2, hmax-2, cmax)

    for coh = 1:cmax, h1 = 1:hmax-2, h2 = max(1,h1-nlag+1):h1
      if (h1+tik-coh >= 1) && (h1+tik-coh <= tmax)
        if (h1!=h2)
          varcov[h1, h2, coh] = covy[h1+2, h1-h2, h1+tik-coh]
        else
          varcov[h1, h2, coh] = vary[h1+2, h1+tik-coh]
        end
      end
    end

    return varcov
  end

  ##############################################################################
  ## Objective function, depending on observed and theoretical covariance
  ## structure
  ##############################################################################

  function objective(params::Vector{Float64}, grad, CovEmp=CovEmp,
    obs_indicator=obs_indicator, weight=Num, tmax=tmax, nlag=nlag, agecell=agecell,
    cmax=cmax, tik=tik, hip=hip)

    # maximum experience: T + 2
    hmax = size(CovEmp,1) + 2

    # varcov = zeros(T,T,size(CovEmp,3))
    varcov = zeros(size(CovEmp))

    ρ = params[1]; varϵ = params[2]; varη = params[3]
    varα = params[4]; varβ = params[5]; corrαβ = params[6]

    varα*varβ != 0 ? covαβ = corrαβ*sqrt(abs(varα*varβ)) : covαβ = 0.0

    π = zeros(tmax)
    ϕ = zeros(tmax)

    π[1] = 1
    π[2:tmax] = params[7:7+tmax-2]
    ϕ[1] = 1
    ϕ[2:tmax-1] = params[7+tmax:7+2*tmax-3]
    ϕ[tmax] = 1

    ϕ .^= 2
    π .^= 2

    # Construct theoretical var-cov matrix
    varcov =
      theoretical_varcov_guv(varα, varβ, covαβ, varϵ, varη, ρ, ϕ, π,
                             hmax, tmax, nlag, hip)

    # matrices to hold observations for
    ∑n = zeros(tmax,tmax)
    # weighted theoretical covariances (from varcov)
    theor = zeros(tmax,tmax)
    # weighted empirical covariances (from CovEmp)
    empir = zeros(tmax,tmax)

    for t1=1:tmax, t2=1:t1
      for coh=1:cmax
        if (t2-tik+coh > 0) && (t1-tik+coh+1 <= hmax-2)
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

    temp = [empir-theor][:]
    obj = [(temp'*temp)][1]
  end

  lb = [[-0.4, 5e-4, 5e-4, 1e-6, 1e-6, -1], 0.3*ones(2*(tmax))]
  ub = [[ 1.2,  2.0,  2.0,  0.5,  0.5,  1],   4*ones(2*(tmax))]

  if use_NLopt == 1
    opt = Opt(:G_MLSL, length(x_0))
    localopt = Opt(:LN_NELDERMEAD, length(x_0))
    local_optimizer!(opt, localopt)
    lower_bounds!(opt, lb)
    upper_bounds!(opt, ub)
    min_objective!(opt, objective)
    ftol_abs!(opt, 1e-12)
    maxeval!(opt, 50000)
    maxtime!(opt, 500)
    (optf, optx, flag) = optimize(opt, x_0)
  else
    optimum = Optim.optimize(objective, x_0, iterations = 1500)
  end
end

function print_results(hip, rip)
  hipf = hip.f_minimum
  ripf = rip.f_minimum
  hipx = hip.minimum
  ripx = rip.minimum
  hipflag = hip.f_converged
  ripflag = rip.f_converged
  @printf "The HIP estimates are:\n Convergence: %s\n Function minimum %.10f\n ρ = %.3f, var(ɛ)= %.3f, var(η) = %.3f,\n var(α) = %.3f, var(β) = %.4f, corr(α,β) = %.3f\n" hipflag hipf hipx[1] hipx[2] hipx[3] hipx[4] hipx[5] hipx[6]
  @printf "The RIP estimates are:\n Convergence: %s\n Function minimum %.10f\n ρ = %.3f, var(ɛ)= %.3f, var(η) = %.3f,\n var(α) = %.3f, var(β) = %.4f, corr(α,β) = %.3f\n" ripflag ripf ripx[1] ripx[2] ripx[3] ripx[4] ripx[5] ripx[6]
end

function print_results(hipflag, hipf, hipx, ripflag, ripf, ripx)
  @printf "The HIP estimates are:\n Exit flag: %s\n Function minimum %.10f\n ρ = %.3f, var(ɛ)= %.3f, var(η) = %.3f,\n var(α) = %.3f, var(β) = %.5f, corr(α,β) = %.2f\n" hipflag hipf hipx[1] hipx[2] hipx[3] hipx[4] hipx[5] hipx[6]
  @printf "Paper Results HIP:\n ρ = 0.821, var(ɛ) = 0.047, var(η) = 0.029\n var(α) = 0.022, var(β) = 0.00038, corr(α,β) = -0.23\n"
  @printf "The RIP estimates are:\n Exit flag: %s\n Function minimum %.10f\n ρ = %.3f, var(ɛ)= %.3f, var(η) = %.3f,\n var(α) = %.3f, var(β) = %.5f, corr(α,β) = %.2f\n" ripflag ripf ripx[1] ripx[2] ripx[3] ripx[4] ripx[5] ripx[6]
  @printf "Paper Results RIP:\n ρ = 0.988, var(ɛ) = 0.061, var(η) = 0.015\n var(α) = 0.058\n"
end
#@time hip =  estimARMA(CovEmp, Num, 1, 0)
#@time rip =  estimARMA(CovEmp, Num, 0, 0)
#print_results(hip, rip)

@time (hipf, hipx, hipflag) = estimARMA(CovEmp, Num, 1, 1)
@time (ripf, ripx, ripflag) = estimARMA(CovEmp, Num, 0, 1)
print_results(hipflag, hipf, hipx, ripflag, ripf, ripx)
