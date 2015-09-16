##############################################################################
## Function to construct theoretical var-cov matrix given parameters of income
## process
##############################################################################

function theoretical_varcov{T<:FloatingPoint}(ρ::T, varϵ::T, varη::T,
      varα::T, varβ::T, covαβ::T, ϕ::Array{T,1},  π::Array{T,1},
      tmax::Int64, nlag::Int64, CovEmp::Array, hip::Int64)

  agemax = size(CovEmp,1); cmax = size(CovEmp,3)
  ip_var = zeros(agemax); varz = zeros(agemax, tmax)
  vary = zeros(agemax, tmax); covy = zeros(agemax, nlag, tmax)

  # Theoretical variances age/year
  for t = 1:tmax
    varz[1, t] = π[t]*varη
    for age = 2:agemax
      expr = age + 2.
      ip_var[age] = varα + hip*( 2*covαβ*expr + varβ*expr^2 )
      varz[age, 1] = ([π[1]*varη for k = 0:(expr-1)]'*[ρ^(2*j) for j = 0:(expr-1)])[1]
      (t < 2) || (varz[age, t] = ρ^2*varz[age-1,t-1] + π[t]*varη)
      vary[age, t] = ip_var[age] + varz[age,t] + ϕ[t]*varϵ
    end
  end

  # Calculate Autocovariances
  for age = 1:agemax, t = 1:tmax
    for n = 1:nlag
      expr = age + 2.
      covy[age, n, t] = varα + hip*( 2*covαβ*(expr+n) + varβ*expr*(expr+n) )
      + (ρ^n)*(varz[age, t])
    end
  end

  varcov = zeros(size(CovEmp))

  for t = 1:tmax, c = 1:cmax
    age = c - tik + t  # How old is this cohort now?
    if (age > 0) && (age < agemax)
      varcov[age, age, c] = vary[age, t]
      for lag = 1:min(tmax-t,nlag)
        if (age+lag)<agemax
          varcov[age+lag,age,c] = covy[age,lag,t]
        end
      end
    end
  end

  return varcov
end

##############################################################################

function theoretical_varcov_guv{T<:FloatingPoint}(ρ::T, varϵ::T, varη::T,
      varα::T, varβ::T, covαβ::T, ϕ::Array{T,1},  π::Array{T,1},
      tmax::Int64, nlag::Int64, CovEmp::Array, hip::Int64)

  agemax = size(CovEmp,1)+2; cmax = size(CovEmp,3)

  ip_part = zeros(agemax)
  vary = zeros(agemax, tmax)

  for t = 1:tmax, age = 3:agemax
    ip_part[age] = varα + hip*(varβ*age^2 + 2*covαβ*age)
    if age <= t
      vary[age, t] =
        (ip_part[age] + ϕ[t]*varϵ
         + sum(varη*(((ρ^2).^collect((age-1):-1:0)).*π[t-age+1:t])) )
    else
      vary[age, t] =
        ( ip_part[age] + ϕ[t]*varϵ
         + sum(varη*((((ρ^2).^collect((t-1):-1:0)).*π[1:t])))
         + sum(varη*π[1]*(((ρ^2).^collect((age-1):-1:t)))) )
    end
  end

  covy = zeros(agemax, nlag, tmax)

  for t= 1:tmax, age = 3:agemax, n = 1:min(age-3,t-1,nlag)
    covy[age, n, t] =
      ( varα + hip*(varβ*age*(age-n) + covαβ*(2*age-n))
       + (ρ^n)*(vary[age-n,t-n]-ip_part[age-n] - ϕ[t-n]*varϵ) )
  end

  varcov = zeros(agemax-2, agemax-2, cmax)

  for coh = 1:cmax, age1 = 1:agemax-2, age2 = max(1,age1-nlag+1):age1
    if (age1+tik-coh >= 1) && (age1+tik-coh <= tmax)
      if (age1!=age2)
        varcov[age1, age2, coh] = covy[age1+2, age1-age2, age1+tik-coh]
      else
        varcov[age1, age2, coh] = vary[age1+2, age1+tik-coh]
      end
    end
  end

  return varcov
end

##############################################################################

function tvg{T<:AbstractFloat}(varα::T, varβ::T, covαβ::T, varϵ::T, varη::T,
  ρ::T, ϕ::Array{T,1}, π::Array{T,1}, hmax::Int64, tmax::Int64, nlag::Int64,
  tik::Int64, hip::Int64)

  ip_part = zeros(hmax)
  vary = zeros(hmax, tmax)

  for t = 1:tmax, h = 3:hmax
    ip_part[h] = varα + hip*(varβ*h^2. + 2*covαβ*h)
    if h <= t
      vary[h, t] =
        (ip_part[h] + ϕ[t]*varϵ
         + sum(varη*(((ρ^2).^collect((h-1):-1:0)).*π[t-h+1:t])) )
    else
      vary[h, t] =
        ( ip_part[h] + ϕ[t]*varϵ
         + sum(varη*((((ρ^2).^collect((t-1):-1:0)).*π[1:t])))
         + sum(varη*π[1]*(((ρ^2).^collect((h-1):-1:t)))) )
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
