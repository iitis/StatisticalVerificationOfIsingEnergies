using CSV
using NPZ



function read_Chimera(f::String)
    csv_reader = CSV.File(f)
    energies = Float64[]

    for i in 1:length(csv_reader)
        if typeof(csv_reader[i][:energy]) == Float64
            push!(energies, csv_reader[i][:energy])
        end
    end
    energies
end


ats = [5, 20, 200, 500, 1000, 2000]
energies = []

for i in 1:length(ats)
    f = "Qfile_short_15_$(ats[i])_css_2.0.csv"
    #f = "Qfile_case1_$(ats[i])_css_2.0.csv"
    x = read_Chimera(f)
    println(length(x))
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
push!(D, "ground" => -11.985714)
#push!(D, "ground" => -30.96)

npzwrite("../chimera_trains_short.npz", D)
#npzwrite("../chimera_trains_case1.npz", D)
