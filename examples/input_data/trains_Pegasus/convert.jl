using CSV
using NPZ



function read_Pegasus(f::String)
    csv_reader = CSV.File(f)
    energies = Float64[]


    for i in 1:length(csv_reader)
        if csv_reader[i][:chain_break_fraction] < 0.5
            for _ in 1:csv_reader[i][:num_occurrences]
                if typeof(csv_reader[i][:energy]) == Float64
                    push!(energies, csv_reader[i][:energy])
                end
            end
        end
    end

    energies
end


ats = [20, 200, 1400]
energies = []

for i in 1:length(ats)
    f = "Qfile_case7_wisla_20_10_1.75_1.75_at_$(ats[i])_css_2.csv"
    x = read_Pegasus(f)
    push!(energies, x)
end

l = minimum([length(e) for e in energies])
ret = zeros(l, length(ats))


for i in 1:length(ats)
    ret[:,i] = energies[i]
end

D = Dict{String, Any}()

push!(D, "energies01" => ret)
push!(D, "annealing_times"  => ats)
push!(D, "ground" => -92.43)


npzwrite("../pegasus_trains_case7.npz", D)
