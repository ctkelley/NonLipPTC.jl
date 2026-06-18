function alg2g(u0, fobj, fgrad, proj, pdata, R, tau, maxit; tol=1.e-6)
    uex = pdata.uex1d
    precond = pdata.precond
    E0 = norm(u0 - uex)
    E = E0
    reshist = Float64[]
    errhist = Float64[]
    v = copy(u0)
    vhalf = copy(v)
    R .= fgrad(v, pdata)
    N0 = norm(R)
    push!(reshist, 1.0)
    push!(errhist, E / E0)
    RX = copy(R)
    Rhalf = copy(R)
    u = copy(u0)
    ux = copy(u)
    fc = fobj(u, pdata)
    rrnrm = 1.0
    grat = 1.0
#    for ix in 1:maxit
    ix=1
#    while (ix <= maxit) && (grat < 2.0)
    while (ix <= maxit) && (rrnrm > tol)
        ix +=1
#
# Two stage method
#
        vhalf .= proj(v - tau * R)
        Rhalf .= fgrad(vhalf, pdata)
        v .= proj(v - tau * Rhalf)
        R .= fgrad(v,pdata)
        ft = fobj(v, pdata)
#        if (ft > fc)
        if false
            u .= ux
        else
            fc = ft
            u .= v
            ux .= u
        end
        rrnrmx = rrnrm
        rrnrm = norm(R) / N0
        grat=rrnrm/rrnrmx
        push!(reshist, rrnrm)
        E = norm(u - uex)
        push!(errhist, E / E0)
    end
    return (sol = u, reshist = reshist, errhist = errhist)
end

function alg2g(GP::GD_Prob, R, tau, maxit; tol=1.e-6)
    u0 = GP.u0
    fobj = GP.fobj
    fgrad = GP.fgrad
    pdata = GP.pdata
    proj = GP.projb
    aout = alg2g(u0, fobj, fgrad, proj, pdata, R, tau, maxit; tol=tol)
    return aout
end
