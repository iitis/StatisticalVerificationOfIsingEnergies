# StatisticalVerificationOfIsingEnergies
implementation of cumulant's based method to verify whether the Ising output contains a ground state energy


# Analysing Ising energy spectrum

### Estimating ground state energy


To estimate ground state energy from data

```julia

estimate_ground_state_energy(energies::Vector{Float64}, α::Float64)


julia> energies = [0.1, 1., 1.2, 2., 2.2, 3., 3.1, 3.5];

julia> estimate_ground_state_energy(energies, 0.19)
9.745023110079497

```

### Bootstrap series of samples of estimated minimal energies

Statistics and histogram can be computed from these


```julia

bootstrap_hists_of_mins(energies::Vector{Float64}, α::Float64, s::Int, l::Int=length(x))

julia> energies = [1. ,1.5, 2., 3., 4.];

julia> Random.seed!(1234);

julia> h = bootstrap_hists_of_mins(energies, 0.19, 3)
3-element Vector{Float64}:
 10.327274688631618
 -4.407048476374387
 -4.721357273543485
```

### Squared error of estimated minimal energy, for comparison with bootstrap

```julia

julia> squared_error(α::Float64, energies::Vector{Float64})

julia> energies = [1. ,1.5, 2., 3., 4.];


julia> squared_error(0.19, energies)
8.633622351519547

```


### Testing energy spectrum

To get the p-value, the probability that the minimal enegry can be the ground state energy (using the Boodstad resampling):

```julia

 bootstrap_get_pvalue(energies::Vector{Float64}, α::Float64, s::Int = 1_000)


 julia> Random.seed!(1234);

 julia> energies = [0.1, 1., 1.2, 2., 2.2, 3., 3.1, 3.5];

 julia> bootstrap_get_pvalue(energies, 0.19)
 0.655


```


### Estimate temperature, given a ground state energy

```julia

estiamte_temperature(ground_energy::Float64, energies::Vector{Float64})


julia> energies = [1. ,1.5, 2., 3., 4.];


julia> estiamte_temperature(0., energies)
0.8958957699200034


```

This project was supported by the Foundation for Polish Science (FNP) under grant number TEAM NET POIR.04.04.00-00-17C1/18-00
