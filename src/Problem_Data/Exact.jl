function lapeval_ex1(r; p0 = 0.5)
    sigma = 2.0 * p0 / (1.0 - p0)
    mask = (r .> 1.0 / 3.0)
    T = abs(mask .* (3.0 * r .- 1.0) * 0.5)
    (T == NaN) && println(T)
    #u=(r-1.0/3.0)*T^sigma
    u = (2.0 / 3.0) * T^(sigma + 1.0)
    dudr = (sigma + 1.0) * T^sigma
    du2dr2 = sigma * (sigma + 1.0) * 1.5 * T^(sigma - 1.0)
    return lapu = du2dr2 + (dudr / r)
end

function lapeval_ex12d(x, y; p0 = 0.5)
    r = sqrt(x^2 + y^2)
    lap2d = lapeval_ex1(r; p0 = p0)
    return lap2d
end

function uefun_ex1(r; p0 = 0.5)
    sigma = 2.0 * p0 / (1.0 - p0)
    T = zeros(size(r))
    mask = (r .> 1.0 / 3.0)
    T = abs.(mask .* (3.0 * r .- 1.0)) * 0.5
    ue2 = (2.0 / 3.0) * T .^ (sigma + 1.0)
    return ue2
end

function uefun_ex12d(x, y; p0 = 0.5)
    r = sqrt(x^2 + y^2)
    u2d = uefun_ex1(r; p0 = p0)
    return u2d
end

function u2dex1(x, y; p0 = 0.5)
    r = sqrt(x^2 + y^2)
    u2d = uefun_ex1(r; p0 = p0)
    return u2d
end
