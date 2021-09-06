using CSV
using StatsBase
using Distributions
using NPZ
using ArgParse

using StatisticalVerificationOfIsingEnergies
import StatisticalVerificationOfIsingEnergies: bootstrap_hists_of_mins, squared_error



function read_Pegasus(f::String)
    csv_reader = CSV.File(f)
    energies = Float64[]

    for i in 1:length(csv_reader)
        if csv_reader[i][:chain_break_fraction] < 0.01
            for _ in 1:csv_reader[i][:num_occurrences]
                push!(energies, csv_reader[i][:energy])
            end
        end
    end
    energies
end


function estimate_min_pegasusu_data(ats::Vector{Int}, css::String, α::Float64, folder::String)

    D = Dict{String, Any}()
    push!(D, "annealing_times" => ats)
    push!(D, "alphas" => α)
    push!(D, "true_ground" => -92.43)

    energies = []

    for i in 1:length(ats)
        f = folder*"/Qfile_case7_wisla_20_10_1.75_1.75_at_$(ats[i])_css_"*css*".csv"
        push!(energies, read_Pegasus(f))

    end

    S = 1_000

    min_data = [minimum(e) for e in energies]
    min_estimated = [estimate_ground_state_energy(e, α) for e in energies]
    p_values = [bootstrap_get_pvalue(e, α, S) for e in energies]

    ys = [bootstrap_hists_of_mins(e, α, S) for e in energies]
    bootstrap_std = [std(y) for y in ys]
    squared_error_from_cums = [squared_error(α, e) for e in energies]

    e_min = D["true_ground"]
    betas = [estiamte_temperature(e_min, e) for e in energies]

    push!(D, "minimum_from_data" => min_data)
    push!(D, "minimum_estimated" => min_estimated)
    push!(D, "p_values" => p_values)
    push!(D, "estimated_betas" => betas)
    push!(D, "bootstrap_std" => bootstrap_std)
    push!(D, "squared_error" => squared_error_from_cums)

    npzwrite("output/case7_wisla.npz", D)
end


ats = [20, 200, 1400]
α = 0.19
estimate_min_pegasusu_data(ats, "2", α, "input_data/trains_data/")
