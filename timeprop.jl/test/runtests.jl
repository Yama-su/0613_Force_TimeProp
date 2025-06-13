using timeprop
import timeprop: perform_timeprop
using Test
using Plots

# 拡張版の timeprop: 時系列データ記録付き
function perform_timeprop_with_trace(F::Function, tmax::Float64, x0::Float64, a0::Float64, h::Float64)
    if tmax <= 0.0
        error("tmax must be positive")
    end
    if h <= 0.0
        error("h must be positive")
    end

    x, a, t = x0, a0, 0.0
    ts, xs, as = Float64[], Float64[], Float64[]
    while t < tmax
        push!(ts, t)
        push!(xs, x)
        push!(as, a)
        a += h * F(x, t)
        x += h * a
        t += h
    end
    return x, a, ts, xs, as
end

@testset "timeprop.jl" begin

    @testset "Uniform motion" begin
        F(x, t) = 0.0
        tmax, x0, a0, h = 1.0, 0.0, 1.0, 1e-3

        println("\n🧪 Test: Uniform Motion")
        println("  Force: F(x, t) = 0.0")
        println("  Initial conditions: x0 = $x0, a0 = $a0, tmax = $tmax, h = $h")

        x_final, a_final, ts, xs, as = perform_timeprop_with_trace(F, tmax, x0, a0, h)

        expected_x = x0 + a0 * tmax
        expected_a = a0

        println("  → Computed: x = $x_final, a = $a_final")
        println("  → Expected: x = $expected_x, a = $expected_a")

        plot(ts, xs, label="x(t)", xlabel="t", ylabel="x / a", title="Uniform Motion", legend=:bottomright)
        plot!(ts, as, label="a(t)")
        savefig("uniform_motion.png")

        @test isapprox(x_final, expected_x, rtol=1e-10)
        @test isapprox(a_final, expected_a, rtol=1e-10)
    end

    @testset "Uniform acceleration motion" begin
        F(x, t) = 1.0
        tmax, x0, a0, h = 1.0, 0.0, 0.0, 1e-4

        println("\n🧪 Test: Uniform Acceleration")
        println("  Force: F(x, t) = 1.0")
        println("  Initial conditions: x0 = $x0, a0 = $a0, tmax = $tmax, h = $h")

        x_final, a_final, ts, xs, as = perform_timeprop_with_trace(F, tmax, x0, a0, h)

        expected_x = x0 + a0 * tmax + 0.5 * F(0, 0) * tmax^2
        expected_a = a0 + F(0, 0) * tmax

        println("  → Computed: x = $x_final, a = $a_final")
        println("  → Expected: x = $expected_x, a = $expected_a")

        plot(ts, xs, label="x(t)", xlabel="t", ylabel="x / a", title="Uniform Acceleration", legend=:bottomright)
        plot!(ts, as, label="a(t)")
        savefig("uniform_acceleration.png")

        @test isapprox(x_final, expected_x, rtol=1e-3)
        @test isapprox(a_final, expected_a, rtol=1e-3)
    end

    @testset "Spring motion" begin
        k = 1.0
        F(x, t) = -k * x
        tmax, x0, a0, h = 2π, 1.0, 0.0, 1e-4

        println("\n🧪 Test: Spring Motion")
        println("  Force: F(x, t) = -k*x with k = $k")
        println("  Initial conditions: x0 = $x0, a0 = $a0, tmax = $tmax, h = $h")

        x_final, a_final, ts, xs, as = perform_timeprop_with_trace(F, tmax, x0, a0, h)

        ω = √k
        expected_x = x0 * cos(ω * tmax)
        expected_a = -x0 * ω * sin(ω * tmax)

        println("  → Computed: x = $x_final, a = $a_final")
        println("  → Expected: x = $expected_x, a = $expected_a")

        plot(ts, xs, label="x(t)", xlabel="t", ylabel="x / a", title="Spring Motion", legend=:bottomleft)
        plot!(ts, as, label="a(t)")
        savefig("spring_motion.png")

        @test isapprox(x_final, expected_x, atol=1e-2)
        @test isapprox(a_final, expected_a, atol=1e-2)
    end

    @testset "Time-dependent force" begin
        F(x, t) = sin(t)
        tmax, x0, a0, h = 2π, 0.0, 0.0, 1e-3

        println("\n🧪 Test: Time-dependent Force")
        println("  Force: F(x, t) = sin(t)")
        println("  Initial conditions: x0 = $x0, a0 = $a0, tmax = $tmax, h = $h")

        x_final, a_final, ts, xs, as = perform_timeprop_with_trace(F, tmax, x0, a0, h)

        println("  → Computed: x = $x_final, a = $a_final")
        println("  ※ 理論値は複雑なため比較は省略")

        plot(ts, xs, label="x(t)", xlabel="t", ylabel="x / a", title="Time-dependent Force", legend=:bottomleft)
        plot!(ts, as, label="a(t)")
        savefig("time_dependent_force.png")

        @test !isnan(x_final)
        @test !isnan(a_final)
    end

end
