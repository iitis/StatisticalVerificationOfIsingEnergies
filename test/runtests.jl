using Test
using Random
using StatisticalVerificationOfIsingEnergies

@testset "estimation minimum" begin

    x = [0.1, 0.1, 0.1, 0.2, 0.2, 0.3]
    @test estimate_ground_state_energy(x, 0.34) ≈ -0.06106446065424451
end



@testset "error calculation" begin
    x = [0.1, 0.1, 0.1, 0.2, 0.2, 0.3]
    @test squared_error(0.34, x) ≈ 0.2614025705829984
end


@testset "bootstrap histogram of minimas" begin
    Random.seed!(1234)
    x = [1., 1., 2., 2., 3., 4., 5.]
    y = bootstrap_hists_of_mins(x, 0.19, 2)
    @test y ≈ [-6.276354333015176, -6.893114986285872]
end

@testset "p-value" begin
    Random.seed!(1234)
    x = [0.1, 1., 1.2, 2., 2.2, 3., 3.1, 3.5];
    @test bootstrap_get_pvalue(x, 0.19) ≈ 0.655
end


@testset "estimate temparature" begin
    x = [1. ,1.5, 2., 3., 4.]
    @test estiamte_temperature(0., x) ≈ 0.8958957699200034

end
