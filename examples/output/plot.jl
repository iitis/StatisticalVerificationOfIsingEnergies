using NPZ
using Plots
using ArgParse
using StatsBase
using HypothesisTests


function plot_minimal_energies(file::String)
        D = npzread(file)

        p = plot(size = (300, 200))

        x = 0.
        try
            x = D["annealing_times"]
            xlabel!("annelaing time μs")
        catch
            x = D["betas"]
            xlabel!("Metropolis Hastings β")
        end

        y =  [D["ground"] for _ in 1:length(x)]
        α = D["alpha"]
        Z = D["minimum_estimated"]

        plot!(p, x, y, label = "true ground", line = (:green, 1.5), legend=false)

        plot!(p, x, D["minimum_from_data"], line = (:red, 1.), marker = (:dot, :red), label = "minimum from D-Wave samples")
        plot!(p, x, Z, label = "minimum estimated", marker = (:dot, :black),  color = "black")

        plot!(p, x, D["bootstrap_q95"], color = "gray", label = "95 percentile of estimated min.", line = (:dash))
        ylabel!("minimal energy")

        str = "minimal_energies"
        file1 = replace(file, ".npz" => str*".pdf")
        savefig(p, file1)

end


function plot_bootstrad_std(file::String)
        D = npzread(file)

        p = plot(size = (300, 200))
        x = 0.

        try
            x = D["annealing_times"]
            xlabel!("annelaing time μs")
        catch
            x = D["betas"]
            xlabel!("Metropolis Hastings β")
        end


        y = D["bootstrap_std"]
        y1 = D["squared_error"]

        α = D["alpha"]


        plot!(p, x, y, label = "std from bootsrap resampling", line = (:green, 1.5), marker = (:dot, :green), legend=(:topright))
        plot!(p, x, y1, label = "std from error calculus", line = (:red, 1.5), marker = (:dot, :red) , legend=(:topright))

        ylabel!("energy")


        str = "_bootstrad_std"
        file1 = replace(file, ".npz" => str*".pdf")
        savefig(p, file1)

end

function plot_minenergy_vs_ground(file::String)

    D = npzread(file)

    p = plot(size = (300, 200))

    x = 0.

    try
        x = D["annealing_times"]
        xlabel!("annelaing time μs")
    catch
        x = D["betas"]
        xlabel!("Metropolis Hastings β")
    end

    println(D["minimum_from_data"])
    println(D["ground"])

    Z = ( D["minimum_from_data"] .- D["ground"])/abs(D["ground"])

    plot!(p, x, Z, markershape = :circle, legend=(:topright), color = "red", label = false)

    ylabel!("(Hₘᵢₙ - H₀)/|H₀|")

    str = "_energies_vsground"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(p, file1)

end

function plot_betas(file::String)

    D = npzread(file)

    x = 0.

    try
        x = D["annealing_times"]
    catch
        x = D["betas"]
    end

    Z = D["estimated_betas"]


    p = plot(x, Z, legend=(:topright), size = (300, 200), color = "red", label = "β")

    plot!(p, x, Z, label = false, marker = (:dot, :red),  color = "red")

    try
        x = D["annealing_times"]
        xlabel!("annelaing time μs")
    catch
        xlabel!("Metropolis Hastings β")
    end
    ylabel!("β estimated")


    str = "_betas"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(p, file1)

end

function plot_p_values(file::String)

    D = npzread(file)

    α = D["alpha"]

    x = 0.
    try
        x = D["annealing_times"]
    catch
        x = D["betas"]
    end

    Z = D["p_values"]

    println(Z)
    Z = Z .+ 1*10^-8


    p = plot(x, Z, markershape = :circle, size = (300, 200), legend=(:topright), label = "α = $α", color = "red", ylims = (10^-5,1.), yaxis=:log)

    try
        x = D["annealing_times"]
        xlabel!("annelaing time μs")
    catch
        xlabel!("Metropolis Hastings β")
    end
    ylabel!("p - value")


    str = "_p_values"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(p, file1)


end



function main(file::String)

    plot_minimal_energies(file)
    plot_bootstrad_std(file)

    plot_betas(file)
    plot_p_values(file)
    plot_minenergy_vs_ground(file)

end


s = ArgParseSettings("description")
  @add_arg_table! s begin
  "file"
  arg_type = String
  help = "file with data to be plotted"
end

file = parse_args(s)["file"]
main(file)
