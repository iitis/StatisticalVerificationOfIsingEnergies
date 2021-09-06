using NPZ
using Plots
using ArgParse
using StatsBase
using HypothesisTests


function plot_vs_at(file::String)
    D = npzread(file)


        x = D["annealing_times"]

        α = D["alpha"]

        y = [D["true_ground"] for _ in 1:length(x)]
        Z = D["minimum_estimated"]


        p = plot(x, y, label = "true ground", line = (:green, 1.5), title = " α = $α, Pegasus" legend=(:topright), ylims = (-150, 200))

        plot!(p, x, D["minimum_from_data"][:,i1], line = (:red, 1.), marker = (:dot, :red), label = "minimum from D-Wave samples")
        plot!(p, x, Z, label = "minimum estimated", marker = (:dot, :black),  color = "black")


        xlabel!("annelaing time μs")
        ylabel!("minimal energy")


        str = "_css_$(css)_$(α)"
        file1 = replace(file, ".npz" => str*".pdf")
        savefig(p, file1)

end


function plot_en_vs_ground(file::String)

    D = npzread(file)
    no_qbits = D["size"]


    x = D["annealing_times"]
    Z =  D["minimum_from_data"] .- D["true_ground"]

    p = plot(x, Z, markershape = :circle, title = " Pegasus", legend=(:topright), color = "red")


    xlabel!("annelaing time μs")
    ylabel!("minimal energy - true ground")


    str = "_energies"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(p, file1)

end

function plot_betas(file::String)

    D = npzread(file)
    no_qbits = D["size"]

    css = D["css"]

    x = D["annealing_times"]
    Z = D["estimated_betas"]


    p = plot(x, Z, label = "css = $(css)", title = " Pegasus", legend=(:topright), color = "red")

    plot!(p, x, Z, label = false, marker = (:dot, :red),  color = "red")


    xlabel!("annelaing time μs")
    ylabel!("β estimated")


    str = "_betas"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(p, file1)

end

function plot_p_values(file::String)

    D = npzread(file)
    no_qbits = D["size"]


    α = D["alpha"]


    x = D["annealing_times"]
    Z = D["p_values"]

    p = plot(x, Z, title = "α = $α, Pegasus", markershape = :circle, legend=(:topright), color = "red", ylims = (-0.1,1.1))


    xlabel!("annelaing time μs")
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
