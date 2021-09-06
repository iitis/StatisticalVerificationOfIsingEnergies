using NPZ
using Plots
using ArgParse
using StatsBase
using HypothesisTests


function plot_vs_at(file::String)
        D = npzread(file)

        x = 0.
        y = 0.

        try
            x = D["annealing_times"]
            y =  [D["true_ground"] for _ in 1:length(x)]
        catch

            x = D["betas"]
            y =  [D["ground"] for _ in 1:length(x)]
        end


        α = D["alpha"]
        Z = D["minimum_estimated"]

        #p = plot(x, y, label = "true ground", line = (:green, 1.5), title = " α = $α, Pegasus", legend=(:topright), ylims = (-150, 200))
        p = plot(x, y, label = "true ground", line = (:green, 1.5), title = " α = $α", legend=(:topright))

        plot!(p, x, D["minimum_from_data"], line = (:red, 1.), marker = (:dot, :red), label = "minimum from D-Wave samples")
        plot!(p, x, Z, label = "minimum estimated", marker = (:dot, :black),  color = "black")


        try
            x = D["annealing_times"]
            xlabel!("annelaing time μs")
        catch
            xlabel!("Metropolis Hastings β")
        end
        xlabel!("Metropolis Hastings β")
        ylabel!("minimal energy")


        str = "_css_$(α)"
        file1 = replace(file, ".npz" => str*".pdf")
        savefig(p, file1)

end


function plot_vs_at(file::String)
        D = npzread(file)


        x = 0.

        try
            x = D["annealing_times"]
        catch
            x = D["betas"]
        end


        y = D["bootstrap_std"]
        y1 = D["squared_error"]

        α = D["alpha"]


        p = plot(x, y, label = "std from bootsrap resampling", line = (:green, 1.5), title = " α = $α", legend=(:topright))
        plot!(p, x, y1, label = "std from error calculus", line = (:red, 1.5), title = " α = $α" , legend=(:topright))


        try
            x = D["annealing_times"]
            xlabel!("annelaing time μs")
        catch
            xlabel!("Metropolis Hastings β")
        end
        ylabel!("energy")


        str = "_bootstrad_std_$(α)"
        file1 = replace(file, ".npz" => str*".pdf")
        savefig(p, file1)

end

function plot_en_vs_ground(file::String)

    D = npzread(file)

    x = 0.
    Z = 0.

    try
        x = D["annealing_times"]
        Z =  D["minimum_from_data"] .- D["true_ground"]
    catch
        x = D["betas"]
        Z =  D["minimum_from_data"] .- D["ground"]
        println(D["ground"])
    end


    p = plot(x, Z, markershape = :circle, legend=(:topright), color = "red", label = "energy")

    try
        x = D["annealing_times"]
        xlabel!("annelaing time μs")
    catch
        xlabel!("Metropolis Hastings β")
    end
    ylabel!("minimal energy - true ground")

    str = "_energies"
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


    p = plot(x, Z, legend=(:topright), color = "red", label = "β")

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

    p = plot(x, Z, title = "α = $α", markershape = :circle, legend=(:topright), label = "p-value", color = "red", ylims = (-0.1,1.1))

    try
        x = D["annealing_times"]
        xlabel!("annelaing time μs")
    catch
        xlabel!("Metropolis Hastings β")
    end
    ylabel!("p - value")


    str = "p_values$(α)"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(p, file1)


end



function main(file::String)

    plot_vs_at(file)

    plot_betas(file)
    plot_p_values(file)
    plot_en_vs_ground(file)

end


s = ArgParseSettings("description")
  @add_arg_table! s begin
  "file"
  arg_type = String
  help = "file with data to be plotted"
end

file = parse_args(s)["file"]
main(file)
