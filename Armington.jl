#! basic implementation

function _cessum(x, shares, sub)
    n = length(x)
    (n == length(shares)) || throw(DimensionMismatch())
    s = zero(eltype(x))
    @inbounds for i in 1:n
        s += shares[i]*x[i]^sub
    end
    return s
end

"""
constant elasticity of substitution function
"""
function ces(x, shares, sub; normalize=false)
    n = length(x)
    (normalize == false) && return _cessum(x, shares, sub)^(1/sub)
    vx = @view(x[1:n-1])
    vshares = @view(shares[1:n-1])
    sharessum = zero(eltype(x))
    @inbounds for i in 1:n-1
        sharessum += shares[i]
    end
    return (_cessum(vx, vshares, sub) + (1 - sharessum) * x[n]^sub)^(1/sub)
end
