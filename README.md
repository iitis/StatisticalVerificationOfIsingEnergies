# StatisticalVerificationOfIsingEnergies
implementation of cumulant's based method to verify whether the Ising output contains a ground state energy


# Analysing Ising energy spectrum

### Estimating ground state energy


To estimate ground state energy from data

```julia

estimate_ground_state_energy(energies::Vector{Float64}, α::Float64)


julia> x = [0.1, 1., 1.2, 2., 2.2, 3., 3.1, 3.5];

julia> estimate_ground_state_energy(x, 0.19)
9.745023110079497

```

### Bootstrap series of samples of estimated minimal energies

Statistics and histogram can be computed from these


```julia

bootstrap_hists_of_mins(x::Vector{Float64}, α::Float64, n::Int, l::Int=length(x))

julia> x = [1. ,1.5, 2., 3., 4.];

julia> h = bootstrap_hists_of_mins(x, 0.19, 10)
10-element Vector{Float64}:
 -30.312380622438372
  -9.205220830116549
  -0.143734959618822
   0.2968807971118199
  48.92088332046563
  -2.429249461334908
   5.837040780975512
   4.3486822264117135
   0.3245649308236207
   0.9309299811472944

julia> std(h)
19.500184919003857

```

### Squared error of estimated minimal energy, for comparison with bootstrap

```julia

julia> squared_error(α::Float64, x::Vector{Float64})

julia> x = [1. ,1.5, 2., 3., 4.];


julia> squared_error(0.19, x)
8.633622351519547


```


### Testing energy spectrum

To get the p-value, the probability that the minimal enegry can be the ground state energy (using the Boodstad resampling):

```julia

 bootstrap_get_pvalue(energies::Vector{Float64}, α::Float64, s::Int = 1_000)



 julia> Random.seed!(1234);

 julia> x = [0.1, 1., 1.2, 2., 2.2, 3., 3.1, 3.5];

 julia> bootstrap_get_pvalue(x, 0.19)
 0.6659999999999999


```


### Estimate temperature, given a ground state energy

```julia

estiamte_temperature(ground_energy::Float64, energies::Vector{Float64})


julia> x = [1. ,1.5, 2., 3., 4.]
5-element Vector{Float64}:
 1.0
 1.5
 2.0
 3.0
 4.0

julia> estiamte_temperature(0., x)
0.8958957699200034


```

This project was supported by the Foundation for Polish Science (FNP) under grant number TEAM NET POIR.04.04.00-00-17C1/18-00
