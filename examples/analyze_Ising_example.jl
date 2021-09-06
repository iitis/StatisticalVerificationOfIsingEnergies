using CSV
using StatsBase
using Distributions
using NPZ
using ArgParse

using StatisticalVerificationOfIsingEnergies
import StatisticalVerificationOfIsingEnergies: bootstrap_hists_of_mins, squared_error


function estimate_min(α::Float64, folder::String)

    D = npzread(folder)

    ens = 0.
    try
        ens = transpose(D["energies01"])
    catch
        ens = D["energies"]
    end

    S = 1_000

    l = size(ens, 1)

    min_data = [minimum(ens[i,:]) for i in 1:l]
    min_estimated = [estimate_ground_state_energy(ens[i,:], α) for i in 1:l]
    p_values = [bootstrap_get_pvalue(ens[i,:], α, S) for i in 1:l]

    ys = [bootstrap_hists_of_mins(ens[i,:], α, S) for i in 1:l]
    bootstrap_std = [std(y) for y in ys]
    squared_error_from_cums = [squared_error(α, ens[i,:]) for i in 1:l]

    e_min = 0.
    try
        e_min = D["true_ground"]
    catch
        e_min = D["ground"]
    end

    betas = [estiamte_temperature(e_min, ens[i,:]) for i in 1:l]

    push!(D, "alpha" => α)
    push!(D, "minimum_from_data" => min_data)
    push!(D, "minimum_estimated" => min_estimated)
    push!(D, "p_values" => p_values)
    push!(D, "estimated_betas" => betas)
    push!(D, "bootstrap_std" => bootstrap_std)
    push!(D, "squared_error" => squared_error_from_cums)

    if contains(file, "artificial")
        outfile = "artificial_case1_wisla.npz"
    else
        outfile = "case7_wisla.npz"
    end

    npzwrite("output/"*outfile, D)
end


s = ArgParseSettings("description")
  @add_arg_table! s begin
  "file"
  arg_type = String
  help = "file with data to be plotted"
end

file = parse_args(s)["file"]
estimate_min(0.19, file)
