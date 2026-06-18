#
# Residual
#
function FEX1(u, pdata)
    precond = pdata.precond
    N = length(u); n = Int(sqrt(N)); h = 1.0 / (1.0 + n)
    D2 = pdata.D2
    FC1 = FGen(u, pdata) - pdata.rhs_eg1
    v1 = reshape(FC1, (n, n))
    fc21 = fish2d(v1, pdata.fdata)
    FC2 = reshape(fc21, (n * n, 1))
    #    FC2 = D2\FC1
    precond ? (FC = FC2) : (FC = FC1)
    return FC
end

#
# Easy to understand Fixed point map
#
function GF1!(G, u, pdata)
    N = length(u); n = Int(sqrt(N)); 
    h=1.0/(n+1); h2=h*h
    tau0 = pdata.tau0
    precond=pdata.precond
    precond ? (tau = tau0) : (tau = tau0*h2)
#    G .= -tau*FEX1(u, pdata)
#    G .+= u
    G .= proj0(u - tau*FEX1(u,pdata))
    return G
end

#
# Fixed point map for two stage method
#
function GF2!(G,u,pdata)
    N = length(u); n = Int(sqrt(N)); 
    h=1.0/(n+1); h2=h*h
    tau0 = pdata.tau0
    precond=pdata.precond
    precond ? (tau = tau0) : (tau = tau0*h2)
#    FC1 = FEX1(u, pdata)
    y = proj0(u - tau*FEX1(u,pdata))
    G .= proj0(u - tau*FEX1(y,pdata))
end
#
# Old fixed point map 
#
 function GFx!(G,u,pdata)
    precond = pdata.precond
    N = length(u); n = Int(sqrt(N)); 
    FC1 = FGen(u, pdata) - pdata.rhs_eg1
    v1 = reshape(FC1, (n, n))
    fc21 = fish2d(v1, pdata.fdata)
    FC2 = reshape(fc21, (n * n, 1))
    precond ? (G .= u - FC2) : (G .= u - FC1)
    return G
end
    


#
# Generic residual
#
function FGen(u, pdata)
    #
    # Problem parameters
    #
    precond = pdata.precond
    nu = pdata.nu
    p = pdata.alpha
    #
    # boundary data
    #
    bvec = pdata.bvec
    #
    # negative Laplacian
    #
    D2 = pdata.D2
    #
    up = projp(u)
    FG = D2 * u + nlterm_ex1.(up, p, nu) - bvec
    #    FG2 = D2\FG1
    #    precond ? (FG=FG2) : (FG=FG1)
    return FG
end
#
# This is the nonlip term in the operator.
# the parameter delta depends on the example
#
function nlterm_ex1(u, p, nu)
    nlout = nu * u^p
    return nlout
end
#
# Problem needs the objective and the projection
#
function fobj(u, pdata)
    #
    # objective function
    #
    chat = pdata.bvec + pdata.rhs_eg1
    p = pdata.alpha
    precond = pdata.precond
    D2 = pdata.D2
    nu = pdata.nu
    if precond == true
        chat .= D2 \ chat
    end
    fcons = dot(u, chat)
    up = projp(u)
    upp1 = (up) .^ (p + 1.0)
    fnlv = (nu / (p + 1.0)) * upp1
    if precond == true
        fnlv .= D2 \ fnlv
    end
    fnl = sum(fnlv)
    D2 = pdata.D2
    fl = 0.5 * (u' * D2 * u)
    if precond == true
        fl = 0.5 * (u' * u)
    end
    fobj = fl + fnl - fcons
    return fobj
end
function projp(y)
    return p = max.(y, 0.0)
end
function proj0(y)
    #    p=max.(y, 0.0)
    return p = y
end
