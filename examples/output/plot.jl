using NPZ
using Plots
using Plots.PlotMeasures
using ArgParse
using StatsBase
using HypothesisTests
using LsqFit


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


function plot_bootstrad_std(file::String)
        D = npzread(file)

        p = plot(size = (300, 200))
        x = 0.

        y = D["bootstrap_std"]
        y1 = D["squared_error"]

        try
            x = D["annealing_times"]
            xlabel!("annelaing time μs")
        catch
            x = D["betas"]
            xlabel!("Metropolis Hastings β")
            plot!(p, ylims = (0, 1.1*maximum(y1[1:10])))
            plot!(p, xaxis = (:log), xlims = (0.04, 1.))
            l = 5
            vline!([x[l]], style = :dot, linewidth = 2., color = "green", label = "model limit")
        end


        α = D["alpha"]


        plot!(p, x, y, label = "bootsrap", line = (:green, 1.5), marker = (:dot, :green), legend=(:topright))
        plot!(p, x, y1, label = "Eq.(12)", line = (:red, 1.5), marker = (:dot, :red) , legend=(:topright))

        ylabel!("std(E₀)")


        str = "_bootstrad_std"
        file1 = replace(file, ".npz" => str*".pdf")
        savefig(p, file1)

end


function plot_skewness(file::String, l::Int = 0)

    D = npzread(file)

    try

        p = plot(size = (300, 200))

        x = D["betas"]
        xlabel!("Metropolis Hastings β")

        l = 5

        z = D["eta"]

        @. model(x, p) = p[1]*x^(p[2]/2)
        p0 = [0.5, 0.5]

        fit = curve_fit(model, x[1:l], z[1:l], p0)

        f(x, a = fit.param[1], b = fit.param[2]) = a*x^(b/2)


        plot!(p, x, z, markershape = :circle, legend=(:topleft), color = "red", label = "epmpirical data", xaxis = (:log), ylims = (0., 2.5))
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
        plot!(p, xaxis = (:log), ylims = (-0.02, 1.1))
        l = 5
        vline!([x[l]], style = :dot, linewidth = 2., color = "green", label = "model limit")
    end


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
    p = plot(size = (300, 200))

    try
        x = D["annealing_times"]
    catch
        x = D["betas"]
        plot!(p, x, x, style = :dot, linewidth = 2., color = "black", label = "expected")
        plot!(xlims = (0.04, 1.), ylims = (0, .7))
        l = 5
        vline!([x[l]], style = :dot, linewidth = 2., color = "green", label = "model limit")
    end

    Z = D["estimated_betas"]


    plot!(p, x, Z, legend=(:bottomright), color = "red", label = "β")

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

    ZZ = ( D["minimum_from_data"] .- D["ground"])/abs(D["ground"])
    li = maximum([1.3*maximum(ZZ), 0.4])

    p = plot(size = (400, 280), ylims = (-0.01, li))
    p1 = twinx()
    α = D["alpha"]

    x = 0.
    try
        x = D["annealing_times"]
        xlabel!("annelaing time μs")
        plot!(p1, ylims = (-0.1, 1.3))
    catch
        x = D["betas"]
        l = 5
        plot!(p, xaxis = (:log), ylims = (-0.01, 1.5))
        plot!(p1, xaxis = (:log), ylims = (-0.1, 1.3))

        vline!([x[l]], style = :dot, linewidth = 2., color = "green", label = "model limit")
        xlabel!("Metropolis Hastings β")
    end

    plot!(p, x, ZZ, markershape = :square, legend=(:topleft), markersise = 10., color = "black", label = "enegies", ylabel = "(Hₘᵢₙ - H₀)/|H₀|", right_margin=12mm)

    Z1 = D["p_values_14"]
    plot!(p1, x, Z1, markershape = :circ, legend=(:topright), label = "p-val., α = 0.14", color = "orange", ylabel = "p - value", right_margin=12mm)

    Z = D["p_values"]
    plot!(p1, x, Z, markershape = :diamond, label = "p-val., α = $α", color = "red")

    Z2 = D["p_values_26"]
    plot!(p1, x, Z2, markershape = :star, label = "p-val., α = 0.26", color = "brown")

    str = "_p_values"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(file1)


end



function main(file::String)

    #plot_minimal_energies(file)
    #plot_bootstrad_std(file)

    #plot_betas(file)
    plot_p_values(file)
    #plot_minenergy_vs_ground(file)
    #plot_skewness(file)

end


s = ArgParseSettings("description")
  @add_arg_table! s begin
  "file"
  arg_type = String
  help = "file with data to be plotted"
end

file = parse_args(s)["file"]
main(file)
