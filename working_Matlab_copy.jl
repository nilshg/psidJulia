tmax = 29; hmax = 43; nlag = 29; hip = 1
varα = 0.02; varβ = 0.00025; corrαβ = -0.23; varϵ = 0.04; varη = 0.02
covαβ = corrαβ*sqrt(varα*varβ); ρ = 0.8
cmax = 31; tik = 11
ϕ = ones(tmax)
π_1 = linspace(1, 1.2, 13)
π_1 = [π_1, linspace(1.2, 1.7, 3)]
π = [1, π_1, linspace(1.7, 1.7, tmax-15)]
π = π.^2

function tvg(varα::Float64, varβ::Float64, covαβ::Float64,
  varϵ::Float64, varη::Float64, ρ::Float64, ϕ::Array{Float64,1},
  π::Array{Float64,1}, hmax::Int64, tmax::Int64, nlag::Int64, hip::Int64)

  ip_part = zeros(hmax)
  vary = zeros(hmax, tmax)

  for t = 1:tmax, h = 3:hmax
    ip_part[h] = varα + hip*(varβ*h^2. + 2*covαβ*h)
    if h <= t
      vary[h,t] =
        (ip_part[h] + ϕ[t]*varϵ
         + sum(varη*(((ρ^2).^[(h-1):-1:0]).*π[t-h+1:t])) )
    else
      vary[h,t] =
        ( ip_part[h] + ϕ[t]*varϵ
         + sum(varη*((((ρ^2).^[(t-1):-1:0]).*π[1:t])))
         + sum(varη*π[1]*(((ρ^2).^[h-1:-1:t]))) )
    end
  end

  covy = zeros(hmax, nlag, tmax)

  for t= 1:tmax, h = 3:hmax, n = 1:minimum([h-3,t-1,nlag])
      covy[h,n,t] =
        ( varα + hip*(varβ*h*(h-n) + covαβ*(2*h-n))
         + (ρ^n)*(vary[h-n,t-n]-ip_part[h-n] - ϕ[t-n]*varϵ) )
  end

  varcov = zeros(hmax-2, hmax-2, cmax)

  for coh = 1:cmax, h1 = 1:hmax-2, h2 = max(1,h1-nlag+1):h1
    if (h1+tik-coh >= 1) && (h1+tik-coh <= tmax)
      if (h1!=h2)
        varcov[h1,h2,coh] = covy[h1+2,h1-h2,h1+tik-coh]
      else
        varcov[h1,h2,coh] = vary[h1+2,h1+tik-coh]
      end
    end
  end

  return varcov
end

varcov = tvg(varα, varβ, covαβ, varϵ, varη, ρ, ϕ, π, hmax, tmax, nlag, hip)