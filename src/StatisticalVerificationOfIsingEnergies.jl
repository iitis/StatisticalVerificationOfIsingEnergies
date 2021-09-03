module StatisticalVerificationOfIsingEnergies

using Distributions
using Cumulants

export estimate_ground_state_energy, squared_error, bootstrap_get_pvalue,
       bootstrap_hists_of_mins, estiamte_temperature

function estimate_ground_state_energy(energies::Vector{Float64}, α::Float64)
    μ = mean(energies)
    σ = std(energies)
    η = skewness(energies)

    return μ-(α+2)/(α+1)*σ/η
end

"""
    squared_error(α::Float64, x::Vector{Float64})

Returns Float64, estimated standard deviation of the estimation of the ground state energy
"""
function squared_error(α::Float64, x::Vector{Float64})
    N = length(x)
    B = (α+2)/(α+1)
    c = [Array(x)[1] for x in cumulants(reshape(x, :, 1), 6, 1)]
    δc2 = sqrt((c[4] + 2*c[2]^2)/N)
    δc3 = sqrt((c[6] + 9*c[2]*c[4]+ 9*c[3]^2 + 6*c[2]^3)/N)
    e2 = B*2*(c[2]/c[3])*δc2
    e3 = B*(c[2]/c[3])^2*δc3
    (e2^2+e3^2)^0.5
end


function bootstrap_hists_of_mins(x::Vector{Float64}, α::Float64, n::Int, l::Int=length(x))
    αs = fill(α, n)
    ret = zeros(length(αs))
    Threads.@threads for α ∈ αs
        estimate_ground_state_energy(rand(x, l), α)
    end
    ret
end

function bootstrap_get_pvalue(energies::Vector{Float64}, α::Float64, s::Int=1000)
    y = bootstrap_hists_of_mins(energies, α, s)
    f = ecdf(y)
    Hmin = minimum(energies)
    1.0- f(Hmin)
end

function estiamte_temperature(ground_energy::Float64, energies::Vector{Float64})
    EH = ground_energy - mean(energies)
    EH/(std(energies)*skewness(energies)*EH-var(energies))
end

end
