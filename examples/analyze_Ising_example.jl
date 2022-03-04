using CSV
using StatsBase
using Distributions
using NPZ
using ArgParse

using StatisticalVerificationOfIsingEnergies
import StatisticalVerificationOfIsingEnergies: bootstrap_hists_of_mins, squared_error


function estimate_min(file::String)

    D = npzread(file)
    α = 0.19

    ens = D["energies01"]

    S = 1_000

    l = size(ens, 2)

    min_data = [minimum(ens[:,i]) for i in 1:l]
    eta = [skewness(ens[:,i]) for i in 1:l]
    min_estimated = [estimate_ground_state_energy(ens[:,i], α) for i in 1:l]
    p_values = [bootstrap_get_pvalue(ens[:,i], α, S) for i in 1:l]
    p_values_10 = [bootstrap_get_pvalue(ens[:,i], 0.10, S) for i in 1:l]
    p_values_39 = [bootstrap_get_pvalue(ens[:,i], 0.39, S) for i in 1:l]

    ys = [bootstrap_hists_of_mins(ens[:,i], α, S) for i in 1:l]
    q95 = [quantile(y, 0.95) for y in ys]
    bootstrap_std = [std(y) for y in ys]
    bootstrap_mean = [mean(y) for y in ys]
    squared_error_from_cums = [squared_error(ens[:,i], α) for i in 1:l]

    e_min = D["ground"]

    betas = [estiamte_temperature(ens[:,i], e_min) for i in 1:l]

    push!(D, "eta" => eta)
    push!(D, "alpha" => α)
    push!(D, "minimum_from_data" => min_data)
    push!(D, "minimum_estimated" => min_estimated)
    push!(D, "p_values" => p_values)
    push!(D, "p_values_10" => p_values_10)
    push!(D, "p_values_39" => p_values_39)
    push!(D, "estimated_betas" => betas)
    push!(D, "bootstrap_std" => bootstrap_std)
    push!(D, "bootstrap_mean" => bootstrap_mean)
    push!(D, "bootstrap_q95" => q95)
    push!(D, "squared_error" => squared_error_from_cums)

    outfile = replace(file, "input_data" => "output")
    npzwrite(outfile, D)
end


s = ArgParseSettings("description")
  @add_arg_table! s begin
  "file"
  arg_type = String
  help = "file with data to be plotted"
end

file = parse_args(s)["file"]
estimate_min(file)
