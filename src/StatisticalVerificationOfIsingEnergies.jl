module StatisticalVerificationOfIsingEnergies

using Distributions
using Cumulants
using StatsBase

export estimate_ground_state_energy, squared_error, bootstrap_get_pvalue,
       bootstrap_hists_of_mins, estiamte_temperature


"""
    estimate_ground_state_energy(energies::Vector{Float64}, α::Float64)

Returns Float64, estimated ground state energy.

```jldoctest

julia> energies = [0.1, 1., 1.2, 2., 2.2, 3., 3.1, 3.5];

julia> estimate_ground_state_energy(energies, 0.19)
9.745023110079497

```
"""
function estimate_ground_state_energy(energies::Vector{Float64}, α::Float64)
    μ = mean(energies)
    σ = std(energies)
    η = skewness(energies)

    return μ-(α+2)/(α+1)*σ/η
end

"""
    squared_error(α::Float64, energies::Vector{Float64})

Returns Float64, estimated standard deviation of the estimation of the ground state energy

```jldoctest

julia> energies = [1. ,1.5, 2., 3., 4.];

julia> squared_error(0.19, energies)
8.633622351519547
```
"""
function squared_error(α::Float64, energies::Vector{Float64})
    N = length(energies)
    B = (α+2)/(α+1)
    c = [Array(x)[1] for x in cumulants(reshape(energies, :, 1), 6, 1)]
    δc2 = sqrt((c[4] + 2*c[2]^2)/N)
    δc3 = sqrt((c[6] + 9*c[2]*c[4]+ 9*c[3]^2 + 6*c[2]^3)/N)
    e2 = B*2*(c[2]/c[3])*δc2
    e3 = B*(c[2]/c[3])^2*δc3
    (e2^2+e3^2)^0.5
end

"""
    bootstrap_hists_of_mins(energies::Vector{Float64}, α::Float64, s::Int, l::Int=length(x))

Returns s-long vector of Float64, estimated ground state energies from l-long samples of dootstrap resampling


julia> energies = [1. ,1.5, 2., 3., 4.];

julia> Random.seed!(1234);

```jldoctest
julia> h = bootstrap_hists_of_mins(energies, 0.19, 3)
3-element Vector{Float64}:
 10.327274688631618
 -4.407048476374387
 -4.721357273543485
```
"""
function bootstrap_hists_of_mins(energies::Vector{Float64}, α::Float64, s::Int, l::Int=length(energies))
    αs = fill(α, s)
    ret = zeros(length(αs))
    i = 1
    Threads.@threads for α ∈ αs
        ret[i] = estimate_ground_state_energy(rand(energies, l), α)
        i = i+1
    end
    ret
end

"""
    bootstrap_get_pvalue(energies::Vector{Float64}, α::Float64, s::Int=1000)

returns Float64, p-value that energies contains the ground state if the Ising sampler,
s is the number ob bootstrap samples.

```jldoctest
julia> Random.seed!(1234);

julia> energies = [0.1, 1., 1.2, 2., 2.2, 3., 3.1, 3.5];

julia> bootstrap_get_pvalue(energies, 0.19)
0.655
```
"""
function bootstrap_get_pvalue(energies::Vector{Float64}, α::Float64, s::Int=1000)
    y = bootstrap_hists_of_mins(energies, α, s)
    f = ecdf(y)
    Hmin = minimum(energies)
    1.0- f(Hmin)
end

"""
    estiamte_temperature(ground_energy::Float64, energies::Vector{Float64})

returns Float64, estimated temperature.


```jldoctest
julia> energies = [1. ,1.5, 2., 3., 4.];

julia> estiamte_temperature(0., energies)
0.8958957699200034
```
"""
function estiamte_temperature(ground_energy::Float64, energies::Vector{Float64})
    EH = ground_energy - mean(energies)
    EH/(std(energies)*skewness(energies)*EH-var(energies))
end

end
