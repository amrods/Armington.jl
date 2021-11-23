#! basic implementation

"""
basic cobb douglas
"""
function cobbdouglas(x, w)
    z = one(eltype(x))
    for i in eachindex(x)
        z *= x[i]^w[i]
    end
    return z
end

"""
basic ces
"""
function _ces(x, w, r)
    z = zero(eltype(x))
    for i in eachindex(x)
        z += w[i] * x[i]^r
    end
    return z^inv(r)
end

"""
constant elasticity of substitution function
"""
function ces(x, w, r)
    n = length(x)
    (n == length(w)) || throw(DimensionMismatch())
    r == zero(eltype(x)) ? cobbdouglas(x, w) : _ces(x, w, r)
end
