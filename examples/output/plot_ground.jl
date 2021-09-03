using NPZ
using Plots
using ArgParse
using StatsBase
using HypothesisTests


function plot_vs_at(file::String)
    D = npzread(file)


    for i2 in 1:length(D["alphas"])
        α = D["alphas"][i2]

        for i1 in 1:length(D["css"])

            css = D["css"][i1]

            x = D["annealing_times"]
            no_qbits = D["size"]

            y = [D["true_ground"] for _ in 1:length(x)]
            Z = D["minimum_estimated"][:,i1,i2]

            p = ()

            if contains(file, "trains")
                print(file)
                if contains(file, "short_15")
                    p = plot(x, y, label = "true ground", line = (:green, 1.5), title = " α = $α, Chimera, short trains, css = $(css)", legend=(:topright), ylims = (-14, -2))
                elseif contains(file, "case2")
                    p = plot(x, y, label = "true ground", line = (:green, 1.5), title = " α = $α, Chimera, case2 trains, css = $(css)")
                elseif contains(file, "case1")
                    p = plot(x, y, label = "true ground", line = (:green, 1.5), title = " α = $α, Chimera, case1 trains, css = $(css)", legend=false)
                else
                    p = plot(x, y, label = "true ground", line = (:green, 1.5), title = " α = $α, Pegasus, $(no_qbits) l. bits, css = $(css)", legend=(:topright), ylims = (-150, 200))
                end
            else
                p = plot(x, y, label = "true ground", line = (:green, 1.5), title = " α = $α, Pegasus, $(no_qbits) l. bits, css = $(css)", legend=(:topright))
            end

            plot!(p, x, D["minimum_from_data"][:,i1], line = (:red, 1.), marker = (:dot, :red), label = "minimum from D-Wave samples")
            plot!(p, x, Z, label = "minimum estimated", marker = (:dot, :black),  color = "black")



            q99 = D["q99"][:,i1,i2]
            q95 = D["q95"][:,i1,i2]
            q75 = D["q75"][:,i1,i2]


            plot!(p, x, q99, color = "brown", label = "percentile 99", line = (:dash))
            plot!(p, x, q95, color = "gray", label = "percentile 95", line = (:dash))
            plot!(p, x, q75, color = "silver", label = "percentile 75", line = (:dash))


            xlabel!("annelaing time μs")
            ylabel!("minimal energy")


            str = "_css_$(css)_$(α)"
            file1 = replace(file, ".npz" => str*".pdf")
            savefig(p, file1)
        end
    end
end


function plot_en_vs_ground(file::String)

    D = npzread(file)
    no_qbits = D["size"]

    p = ()

    for i1 in 1:length(D["css"])

        css = D["css"][i1]

        x = D["annealing_times"]
        Z =  D["minimum_from_data"][:,i1] .- D["true_ground"]

        if css == 2.
            if contains(file, "short_15")
                p = plot(x, Z, label = "short, css = $(css)", title = " Chimera trains, $(no_qbits) l. bits", legend=(:topright), color = "red", ylims = (0, 17))
            elseif contains(file, "case2")
                p = plot(x, Z, label = "css = $(css)", title = " Chimera case2 trains, $(no_qbits) l. bits", legend=(:topright), color = "red", ylims = (0, 17))
            elseif contains(file, "case1")

                file1 = replace(file, "case1" => "case2")
                D1 = npzread(file1)

                Z1 = D1["minimum_from_data"][:,i1] .- D1["true_ground"]

                p = plot(x, Z, label = "case 1, css = $(css)", title = " Chimera trains, $(no_qbits) l. bits", markershape = :circle, legend=(:topright), color = "red", ylims = (0, 17))

                plot!(p, x, Z1, label = "case 2, css = $(css)", markershape = :circle, color = "blue")
            elseif contains(file, "setcov")
                p = plot(x, Z, label = "css = $(css)", markershape = :circle, title = " setcover, Pegasus, $(no_qbits) l. bits", legend=(:topright), color = "red")
            else
                p = plot(x, Z, label = "css = $(css)", markershape = :circle, title = " Pegasus, $(no_qbits) l. bits", legend=(:topright), color = "red")
            end
        elseif css == 3.5
            plot!(p, x, Z, label = "css = $(css)", markershape = :circle,  color = "blue")
        end
    end

    xlabel!("annelaing time μs")
    ylabel!("minimal energy - true ground")


    str = "_energies"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(p, file1)

end

function plot_betas(file::String)

    D = npzread(file)
    no_qbits = D["size"]

    p = ()

    for i1 in 1:length(D["css"])

        css = D["css"][i1]

        x = D["annealing_times"]
        Z = D["estimated_betas"][:,i1]

        if css == 2.
            if contains(file, "short_15")
                p = plot(x, Z, label = "css = $(css)", title = " Chimera trains, $(no_qbits) l. bits", legend=(:topright), color = "red")
            elseif contains(file, "case2")
                p = plot(x, Z, label = "css = $(css)", title = " Chimera case2 trains, $(no_qbits) l. bits", legend=(:topright), color = "red")
            elseif contains(file, "case1")

                file1 = replace(file, "case1" => "case2")
                D1 = npzread(file1)

                Z1 = D1["estimated_betas"][:,i1]

                p = plot(x, Z, label = "case 1, css = $(css)", title = " Chimera trains, $(no_qbits) l. bits", markershape = :circle, legend=(:topright), color = "red")

                plot!(p, x, Z1, label = "case 2, css = $(css)", markershape = :circle, color = "blue")

            else
                p = plot(x, Z, label = "css = $(css)", title = " Pegasus, $(no_qbits) l. bits", legend=(:topright), color = "red")
            end
            plot!(p, x, Z, label = false, marker = (:dot, :red),  color = "red")
        elseif css == 3.5
            plot!(p, x, Z, label = "css = $(css)", title = " Pegasus, $(no_qbits) l. bits",  color = "blue")
            plot!(p, x, Z, label = false, marker = (:dot, :blue),  color = "blue")
        end
    end

    xlabel!("annelaing time μs")
    ylabel!("β estimated")


    str = "_betas"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(p, file1)

end

function plot_p_values(file::String)


    D = npzread(file)
    no_qbits = D["size"]

    p = ()

    for i2 in 1:length(D["alphas"])
        α = D["alphas"][i2]

        for i1 in 1:length(D["css"])

            css = D["css"][i1]

            x = D["annealing_times"]
            Z = D["p_values"][:,i1, i2]

            if css == 2.
                if contains(file, "short_15")
                    p = plot(x, Z, label = "short, css = $(css)", title = "α = $α,  Chimera trains, $(no_qbits) l. bits", markershape = :circle, legend=(:topright), color = "red", ylims = (-0.1,1.1))
                elseif contains(file, "case2")
                    p = plot(x, Z, label = "css = $(css)", title = "α = $α,  Chimera case2 trains, $(no_qbits) l. bits", markershape = :circle, legend=(:topright), color = "red", ylims = (-0.1,1.1))
                elseif contains(file, "case1")

                    file1 = replace(file, "case1" => "case2")
                    D1 = npzread(file1)

                    Z1 = D1["p_values"][:,i1, i2]

                    p = plot(x, Z, label = "case1, css = $(css)", title = "α = $α, Chimera trains, $(no_qbits) l. bits", markershape = :circle, legend=(:topright), color = "red", ylims = (-0.1,1.1))
                    plot!(p, x, Z1, label = "case2, css = $(css)", markershape = :circle, color = "blue")

                elseif contains(file, "setcov")
                    p = plot(x, Z, label = "css = $(css)", title = "α = $α, stecover Pegasus, $(no_qbits) l. bits", markershape = :circle, legend=(:topright), color = "red", ylims = (-0.1,1.1))
                else
                    p = plot(x, Z, label = "css = $(css)", title = "α = $α, Pegasus, $(no_qbits) l. bits", markershape = :circle, legend=(:topright), color = "red", ylims = (-0.1,1.1))
                end

            elseif css == 3.5
                plot!(p, x, Z, label = "css = $(css)", title = "α = $α, Pegasus, $(no_qbits) l. bits",  markershape = :circle, color = "blue")

            end
        end

        xlabel!("annelaing time μs")
        ylabel!("p - value")


        str = "p_values$(α)"
        file1 = replace(file, ".npz" => str*".pdf")
        savefig(p, file1)

    end

end


function plot_betas_vs_MH(file::String, mhβ::Float64)

    D = npzread(file)
    no_qbits = D["size"]

    p = ()

    for i1 in 1:length(D["css"])

        css = D["css"][i1]

        x = D["annealing_times"]
        Z = D["estimated_betas"][:,i1]./mhβ

        if css == 2.
            p = plot(x, Z, label = "css = $(css)", title = " Pegasus, $(no_qbits) l. bits", legend=(:topright), xaxis=:log)
        elseif css == 3.5
            plot!(p, x, Z, label = "css = $(css)", title = " Pegasus, $(no_qbits) l. bits", xaxis=:log)
        end
    end

    xlabel!("annelaing time μs")
    ylabel!("β estimated / β from MH")


    str = "_betas_vsMH"
    file1 = replace(file, ".npz" => str*".pdf")
    savefig(p, file1)

end



function main(file::String)

    plot_vs_at(file)


    if contains(file, "setcov/07_medium_small")

        plot_betas_vs_MH(file, 0.00144)
    end

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
