##############################################################################
## Function to construct theoretical var-cov matrix given parameters of income
## process
##############################################################################

function theoretical_varcov(varα::Float64, varβ::Float64, covαβ::Float64,
                            ρ::Float64, varη::Float64, varϵ::Float64, π::Array{Float64,1},
                            ϕ::Array{Float64,1}, hmax::Int64, tmax::Int64, nlag::Int64, hip::Int64)

  # Calculate CovEmpariances
  ip_part_var = zeros(hmax)
  varz = zeros(hmax, tmax)
  vary = zeros(hmax, tmax)
  covary = zeros(hmax, nlag, tmax)

  for t = 1:tmax
    varz[1, t] = π[t]*varη
    for h = 1:hmax-2
      ip_part_var[h] = varα + hip*2*covαβ*h + hip*varβ*h^2
      varz[h, 1] = [[π[1]*varη for k = 0:(h-1)]'*[ρ^(2*j) for j = 0:(h-1)]][1]
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