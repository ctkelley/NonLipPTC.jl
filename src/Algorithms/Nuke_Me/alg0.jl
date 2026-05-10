function alg0(u0, fobj, fgrad, proj, pdata, R, tau, maxit)
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
    dv = copy(u)
    fc = fobj(u, pdata)
    NC = N0
    Erest = 1.0
    for ix in 1:maxit
        dv .= u
        v .= proj(v - tau * R)
        dv .-= v
        R .= fgrad(v, pdata)
        NT = norm(R)
        ft = fobj(v, pdata)
        #        if (ft > fc)
        if (NT > NC)
            u .= ux
        else
            #            println(norm(v))
            fc = ft
            NC = NT
            u .= v
            ux .= u
            Erest = norm(dv, Inf)
        end
        rrnrm = NC / N0
        #        rrnrm=norm(R) / N0
        push!(reshist, rrnrm)
        push!(errhist, Erest)
        E = norm(u - uex)
    end
    return (sol = u, reshist = reshist, errhist = errhist)
end
