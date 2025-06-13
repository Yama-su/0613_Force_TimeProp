module timeprop

function perform_timeprop(F::Function, tmax::Float64, x0::Float64, a0::Float64, h::Float64)::Tuple{Float64, Float64}
    if tmax <= 0.0
        error("tmax must be positive")
    end
    if h <= 0.0
        error("h must be positive")
    end

    x::Float64 = x0
    a::Float64 = a0
    t::Float64 = 0.0
    while t < tmax
        a = a + h * F(x, t)
        x = x + h * a
        t = t + h
    end
    return (x, a)
end
println("This is a test run")

# 関数を定義しているならそれを使って出力
# println(my_function(123))

end
