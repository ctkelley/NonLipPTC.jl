function alg1(u0, fobj, fgrad, proj, pdata, R, tau, maxit; epsilon = 0.001)
    uex = pdata.uex1d
    E0 = norm(u0 - uex)
    E = E0
    reshist = Float64[]
    errhist = Float64[]
    v = copy(u0)
    R .= fgrad(v, pdata)
    N0 = norm(R)
    push!(reshist, 1.0)
    push!(errhist, E / E0)
    RX = copy(R)
    u = copy(u0)
    ux = copy(u)
    fc = fobj(u, pdata)
    ix = 1
    rrnrm = 1.0
    while ((ix <= maxit) && (rrnrm > epsilon))
        #    for ix = 1:maxit
        v .= proj(v - tau * R)
        R .= fgrad(v, pdata)
        ft = fobj(v, pdata)
        #        if (ft > fc)
        if (norm(R) > norm(RX))
            u .= ux
            R .= RX
        else
            RX .= R
            fc = ft
            u .= v
            ux .= u
        end
        rrnrm = norm(R) / N0
        push!(reshist, rrnrm)
        E = norm(u - uex)
        push!(errhist, E / E0)
        ix += 1
    end
    return (sol = u, reshist = reshist, errhist = errhist)
end

function alg1(GP::GD_Prob, R, tau, maxit; epsilon = 0.001)
    u0 = GP.u0
    fobj = GP.fobj
    fgrad = GP.fgrad
    pdata = GP.pdata
    proj = GP.projb
    aout = alg1(u0, fobj, fgrad, proj, pdata, R, tau, maxit; epsilon = epsilon)
    return aout
end
