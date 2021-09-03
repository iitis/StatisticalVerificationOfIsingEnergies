using DataFrames
using CSV
using Plots
using Random
using StatsBase
using Statistics
using Distributions
using HypothesisTests
using NPZ
using ArgParse
using Distributed
using DataFrames
using DelimitedFiles
using LightGraphs
using StatisticalVerificationOfIsingEnergies
import StatisticalVerificationOfIsingEnergies: bootstrap_hists_of_mins, squared_error
#using SpinGlassMetropolisHastings
#import SpinGlassMetropolisHastings: bootstrap_hists_of_mins, squared_error


function read_Pegasus(f::String)
    println(f)
    csv_reader = CSV.File(f)
    system_size = length(keys(csv_reader[1]))-3
    energies = Float64[]
    sols = Vector{Int}[]
    ke = keys(csv_reader[1])

    for i in 1:length(csv_reader)
        if csv_reader[i][:chain_break_fraction] < 0.01
            for _ in csv_reader[i][:num_occurrences]
                sol =  [csv_reader[i][k] for k in ke[1:end-3]]
                push!(energies, csv_reader[i][:energy])
                push!(sols, sol)
            end
        end
    end
    energies, system_size, sols
end



function read_chimera(f::String)
    csv_reader = CSV.File(f)

    ke = keys(csv_reader[1])

    e = zeros(length(csv_reader))
    sols = Vector{Int}[]
    for i in 1:length(csv_reader)

        sol =  [csv_reader[i][k] for k in ke[1:end-2]]
        push!(sols, sol)
        e[i] = csv_reader[i][:energy]
    end
    e, length(sols[1]), sols
end



function file_names(folder::String, at::Int, css::String, p_size::String)
    if p_size == ""
        f = folder*"/Qfile_case7_wisla_20_10_1.75_1.75_at_$(at)_css_"*css*".csv"
    end
    f
end


function make_dict(ats::Vector{Int}, all_css::Vector{String}, αs::Vector{Float64}, problem_type::String, p_size::String)
    D = Dict{String, Any}()
    push!(D, "annealing_times" => ats)
    push!(D, "css" => [parse(Float64, css) for css in all_css])
    push!(D, "alphas" => αs)

    if occursin("trains", problem_type)
        push!(D, "true_ground" => -92.43)
    end
    D
end



function estimate_min_pegasusu_data(ats::Vector{Int}, all_css::Vector{String}, αs::Vector{Float64}, problem_type::String, p_size::String)

    D = make_dict(ats, all_css, αs, problem_type, p_size)


    min_data = zeros(length(ats), length(all_css))
    min_estimated = zeros(length(ats), length(all_css), length(αs))
    p_values = zeros(length(ats), length(all_css), length(αs))

    bootstrap_std = zeros(length(ats), length(all_css), length(αs))
    squared_error_from_cums = zeros(length(ats), length(all_css), length(αs))


    q99 = zeros(length(ats), length(all_css), length(αs))
    q95 = zeros(length(ats), length(all_css), length(αs))
    q75 = zeros(length(ats), length(all_css), length(αs))
    betas = zeros(length(ats), length(all_css))

    system_size = 0


    for i in 1:length(ats)
        at = ats[i]

        for i1 in 1:length(all_css)

            css = all_css[i1]

            f = file_names(problem_type, at, all_css[i1], p_size)

            energies, system_size, sols = read_Pegasus(f)

            min_data[i, i1] = minimum(energies)

            for i2 in 1:length(αs)
                α = αs[i2]

                mini = estimate_ground_state_energy(energies, α)
                min_estimated[i, i1, i2] = mini

                y = bootstrap_hists_of_mins(energies, α, 1_000)

                q75[i, i1, i2] = quantile(y, 0.75)
                q95[i, i1, i2] = quantile(y, 0.95)
                q99[i, i1, i2] = quantile(y, 0.99)
                p_values[i, i1, i2] = bootstrap_get_pvalue(energies, α, 1_000)

                bootstrap_std[i, i1, i2] = std(y)
                squared_error_from_cums[i, i1, i2] = squared_error(α, energies)


            end

            e_min = D["true_ground"]
            betas[i, i1] = βfromcumulants(e_min, energies)


        end
    end

    if p_size == ""
        f2 = "output/case7_wisla.npz"
    end

    push!(D, "size" => system_size)
    push!(D, "estimated_betas" => betas)
    push!(D, "minimum_from_data" => min_data)
    push!(D, "minimum_estimated" => min_estimated)
    push!(D, "p_values" => p_values)

    push!(D, "bootstrap_std" => bootstrap_std)
    push!(D, "squared_error" => squared_error_from_cums)

    push!(D, "q75" => q75)
    push!(D, "q95" => q95)
    push!(D, "q99" => q99)

    npzwrite(f2, D)
end


function main(folder::String, αs::Vector{Float64})


    if occursin("trains", folder)
        println("trains")
        all_css = ["2"]
        ats = [20, 200, 1400]
        estimate_min_pegasusu_data(ats, all_css, αs, folder, "")
    end
end


s = ArgParseSettings("description")
  @add_arg_table! s begin
  "folder"
  default = "input_data/trains_data/"
  arg_type = String
  help = "folder with pegasus data"
  "--alphas", "-a"
  nargs = '*'
  default = [0.19]
  arg_type = Float64
end

folder = parse_args(s)["folder"]

as = parse_args(s)["alphas"]


main(folder, as)
