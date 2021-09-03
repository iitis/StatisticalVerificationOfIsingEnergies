using Test
using StatisticalVerificationOfIsingEnergies

@testset "estimation minimum" begin

    x = [0.1, 0.1, 0.1, 0.2, 0.2, 0.3]
    @test estimate_ground_state_energy(x, 0.34) ≈ -0.06106446065424451
end



@testset "error calculation" begin
    x = [0.1, 0.1, 0.1, 0.2, 0.2, 0.3]
    @test squared_error(0.34, x) ≈ 0.2614025705829984
end
