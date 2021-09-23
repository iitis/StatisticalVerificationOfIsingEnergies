using NPZ
using Plots
using Plots.PlotMeasures
using ArgParse
using StatsBase
using HypothesisTests
using LsqFit
using Compose


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
            plot!(xscale = (:log))
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


function plot_bootstrad_std(file::String, l::Int)
        D = npzread(file)

        p = plot(size = (400, 250))
        x = 0.

        y = D["bootstrap_std"]
        y1 = D["squared_error"]

        try
            x = D["annealing_times"]
            xlabel!("annelaing time μs")
        catch
            x = D["betas"]
            xlabel!("Metropolis Hastings β")
            #plot!(p, ylims = (0, 1.1*maximum(y1[1:10])))
            plot!(p, xaxis = (:log))
            #vline!([x[l]], style = :dot, linewidth = 2., color = "green", label = "model limit")
        end


        α = D["alpha"]


        plot!(p, x, y, label = "bootsrap", line = (:green, 1.5), marker = (:dot, :green), legend=(:topright))
        plot!(p, x, y1, label = "Eq.(12)", line = (:red, 1.5), marker = (:dot, :red) , legend=(:topright))

        ylabel!("std(E₀)")


        str = "_bootstrad_std"
        file1 = replace(file, ".npz" => str*".pdf")
        savefig(p, file1)

end


function plot_skewness(file::String, l::Int )

    D = npzread(file)

    try
        p = plot(size = (300, 200))
        x = D["betas"]
        xlabel!("Metropolis Hastings β")

        z = D["eta"]

        @. model(x, p) = p[1]*x^(p[2]/2)
        p0 = [0.5, 0.5]

        fit = curve_fit(model, x[1:l], z[1:l], p0)
        f(x, a = fit.param[1], b = fit.param[2]) = a*x^(b/2)

        plot!(p, x, z, markershape = :circle, legend=(:topleft), color = "red", label = "epmpirical data", xaxis = (:log))
        #plot!(p, ylims = (0., 2.5))
        plot!(p, x, [f(i) for i in x], linewidth = 2, style = :dot, color = "black", label = "fit with α = $(round(fit.param[2], digits=2))")

        vline!([x[l]], style = :dot, linewidth = 2., color = "green", label = "fitting limit β = $(x[l])")

        ylabel!("η(H)")

        str = "_etas"
        file1 = replace(file, ".npz" => str*".pdf")
        savefig(p, file1)


    catch
        0
    end

end



function plot_betas(file::String, l::Int)

    D = npzread(file)
    x = 0.
    p = plot(size = (300, 200))

    Z = D["estimated_betas"]

    try
        x = D["annealing_times"]
    catch
        x = D["betas"]
        #plot!(p, x, x, style = :dot, linewidth = 2., color = "black", label = "1 to 1")
        #plot!(xlims = (0.04, 2.), ylims = (0, 1.))
        vline!([x[l]], style = :dot, linewidth = 2., color = "green", label = "model limit")

        @. model(x, p) = p[1]+p[2]*x
        p0 = [0.5, 0.5]

        fit = curve_fit(model, x[1:l], Z[1:l], p0)
        f(x, a = fit.param[1], b = fit.param[2]) = a+x*b

        plot!(p, x, [f(i) for i in x], linewidth = 2, style = :dot, color = "blue", label = "linear model, slope $(round(fit.param[2], digits = 2))")
    end


    plot!(p, x, Z, legend=(:topright), color = "red", label = "β")
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

function plot_p_values(file::String, l::Int)

    D = npzread(file)
    ZZ = ( D["minimum_from_data"] .- D["ground"])/abs(D["ground"])
    li = maximum(ZZ)
    p = plot(size = (400, 250))
    p1 = twinx()
    α = D["alpha"]

    x = 0.
    try
        x = D["annealing_times"]
        xlabel!("annelaing time μs")
    catch
        x = D["betas"]
        plot!(p, xaxis = (:log))
        vline!([x[l]], style = :dot, linewidth = 2., color = "green", label = "model limit")

        xlabel!("Metropolis Hastings β")
    end

    plot!(p, legend=(0.15, 0.45), ylims = (-0.01, 1.05*li))
    plot!(p1, ylims = (-0.01, 1.01), legend=(0.67, 0.45))
    plot!(p, x, ZZ, markershape = :square, linewidth = 2, markersise = 10., color = "blue", label = "Hs", ylabel = "(Hₘᵢₙ - H₀)/|H₀|", right_margin=12mm)

    Z1 = D["p_values_10"]
    plot!(p1, x, Z1, markershape = :circ, label = "p-val., α = 0.10", color = "orange", ylabel = "p - value", right_margin=12mm)

    Z = D["p_values"]
    plot!(p1, x, Z, markershape = :diamond, label = "p-val., α = $α", color = "red")

    Z2 = D["p_values_39"]
    plot!(p1, x, Z2, markershape = :star, label = "p-val., α = 0.39", color = "brown")

    str = "_p_values"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(file1)
end


function main(file::String)
    β_lim_ind = 0
    if occursin("artificial_trains_case1", file)
        β_lim_ind = 8
    elseif occursin("artificial_trains_short", file)
        β_lim_ind = 5
    end
    plot_minimal_energies(file)
    plot_bootstrad_std(file, β_lim_ind)
    plot_betas(file, β_lim_ind)
    plot_p_values(file, β_lim_ind)
    plot_skewness(file, β_lim_ind)

end


s = ArgParseSettings("description")
  @add_arg_table! s begin
  "file"
  arg_type = String
  help = "file with data to be plotted"
end

file = parse_args(s)["file"]
main(file)
