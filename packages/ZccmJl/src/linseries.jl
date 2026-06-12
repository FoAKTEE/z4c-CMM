# linseries.jl — exact retarded-series algebra over Rational{BigInt} in the
# basis F^(n)(u) / r^k (port of the A1/A2 dict machinery; the workhorse for
# linear characteristic derivations).

export LinSeries, dr_series, du_series, mul_rpow, series_add, series_scale,
       integrate_r

const LinSeries = Dict{Tuple{Int,Int},Rational{BigInt}}

series_add(a::LinSeries, b::LinSeries) = begin
    out = copy(a)
    for (k, c) in b
        v = get(out, k, 0//1) + c
        v == 0 ? delete!(out, k) : (out[k] = v)
    end
    out
end

series_scale(a::LinSeries, f::Union{Rational,Integer}) = begin
    ff = Rational{BigInt}(f)
    LinSeries(k => c*ff for (k, c) in a if c*ff != 0)
end

"""d/dr at fixed u (only the r^-k factors differentiate; retarded basis)."""
dr_series(a::LinSeries) =
    LinSeries((n, k+1) => -k*c for ((n, k), c) in a if k != 0)

"""d/du: F^(n) -> F^(n+1)."""
du_series(a::LinSeries) = LinSeries((n+1, k) => c for ((n, k), c) in a)

mul_rpow(a::LinSeries, p::Int) =
    LinSeries((n, k - p) => c for ((n, k), c) in a)

"""Antiderivative in r: c r^{1-k}/(1-k); the log channel (k = 1) throws —
it must be handled as a closure condition upstream."""
function integrate_r(a::LinSeries, label::String)
    out = LinSeries()
    for ((n, k), c) in a
        k == 1 && error("log channel in $label at n=$n")
        out[(n, k-1)] = c//(1 - k)
    end
    out
end

"""Evaluate a series numerically at (u, r) with pulse params (BigFloat-safe)."""
function evaluate(a::LinSeries, u::T, r::T, X::T, rc::T, tau::T) where {T<:Real}
    out = zero(T)
    for ((n, k), c) in a
        out += T(c)*teuk_F(n, u, X, rc, tau)/r^k
    end
    out
end
export evaluate
