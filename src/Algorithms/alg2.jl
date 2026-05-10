function alg2(u0, fobj, fgrad, proj, pdata, R, tau, epsilon, maxit)
    uex = pdata.uex1d
    E0 = norm(u0 - uex)
    E = E0
    mu = pdata.mu
    reshist = Float64[]
    errhist = Float64[]
    push!(reshist, 1.0)
    push!(errhist, E / E0)
    u = copy(u0)
    fc = fobj(u, pdata)
    v = copy(u0)
    ux = copy(u0)
    vx = copy(u0)
    R .= fgrad(v, pdata)
    RX = copy(R)
    N0 = norm(R)
    N = length(u)
    #tau=tau0/N
    #    tau=tauev(tau0, N)
    f_best = fc
    for ix in 1:maxit
        v .= proj0(v - tau * R)
        ft = fobj(v, pdata)
        df = ft - fc
        # The sufficient decrease condition is WEAKER that simple decrease
        # because of the extra positive terms
        SD1 = dot(RX, v - vx) + 0.5 * dot(v - vx, v - vx) / tau + mu * epsilon^2 / 4
        R .= fgrad(v, pdata)
        if (df > SD1)
            v .= vx
            R .= RX
            tau *= 0.5
        else
            RX .= R
            vx .= v
            fc = ft
            if (f_best < fc)
                u .= ux
            else
                f_best = fc
                u .= v
                ux .= u
            end
        end
        rrnrm = norm(R) / N0
        push!(reshist, rrnrm)
        E = norm(u - uex)
        push!(errhist, E / E0)
    end
    return alg2ouot = (sol = u, reshist = reshist, errhist = errhist)
end

function alg2(GP::GD_Prob, R, tau, epsilon, maxit)
    u0 = GP.u0
    fobj = GP.fobj
    fgrad = GP.fgrad
    pdata = GP.pdata
    proj = GP.projb
    #    mu = GP.mu
    return aout = alg2(u0, fobj, fgrad, proj, pdata, R, tau, epsilon, maxit)
end
