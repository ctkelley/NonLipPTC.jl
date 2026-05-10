#
# Solve the BVP and compare to the exact solution
#
function LSolve(rhs, uefun; p0 = 0.5)
    N = length(rhs)
    n = Int(sqrt(N))
    rhs = fix_rhs!(rhs, uefun; p0 = p0)
    D2 = Lap2d(n)
    Psol = D2 \ rhs
    return Psol
end

function fix_rhs!(rhs, u2d; p0 = 0.5)
    N = length(rhs)
    n = Int(sqrt(N))
    h = 1.0 / (n + 1.0)
    X = h:h:(1.0 - h)
    ulow = u2d.(X, 0.0)
    uhigh = u2d.(X, 1.0)
    uleft = u2d.(0.0, X)
    uright = u2d.(1.0, X)
    add_boundary!(rhs, -uleft, -uright, -ulow, -uhigh)
    return rhs
end

function add_boundary!(u, uleft, uright, ulow, uhigh)
    N = length(u)
    n = Int(sqrt(N))
    hm2 = (n + 1) * (n + 1)
    u2 = reshape(u, (n, n))
    for i in 1:n
        u2[i, 1] -= hm2 * ulow[i]
        u2[i, n] -= hm2 * uhigh[i]
        u2[1, i] -= hm2 * uleft[i]
        u2[n, i] -= hm2 * uright[i]
    end
    return u = reshape(u2, (N, 1))
end
