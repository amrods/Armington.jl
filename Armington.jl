#! basic implementation

"""
constant elasticity of substitution function
"""
function ces(x, w, r)
    n = length(x)
    (n == length(w)) || throw(DimensionMismatch())
    
    if r == zero(eltype(x))
        z = one(eltype(x))
        for i in 1:n
            z *= x[i]^w[i]
        end
        return z
    else
        z = zero(eltype(x))
        for i in 1:n
            z += w[i] * x[i]^r
        end
        return z^inv(r)
    end
end
