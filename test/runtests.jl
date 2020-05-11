using Powells_minimization_travis
using Test

f(x) = (x-3)^2+5
g(x) = x^5+3*x^2+x+9
h(x) = sin(x)

@testset "Powells_minimization_travis.jl" begin
    @test Powell_minimize(f, 8)==3
    @test Powell_minimize(f, -5)==3
    @test round(Powell_minimize(g, -1), digits=5)==-0.16732
    @test round(Powell_minimize(g, 2), digits=5)==-0.16732
    @test round(Powell_minimize(g, 5), digits=5)==-0.16732
    @test round(Powell_minimize(h, 1), digits=5)==round(-π/2, digits=5)
end
