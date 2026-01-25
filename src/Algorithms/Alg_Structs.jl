struct GD_Prob{T<:Real}
    fobj::Function
    fgrad::Function
    projb::Function
    pdata::NamedTuple
    u0::Array{T,1}
    mu::T
end
