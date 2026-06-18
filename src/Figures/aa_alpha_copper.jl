function aa_alpha_copper(
        n; nu = 1.0, tau=.1, tol=1.e-12,
        maxit = 20000, pvec = [0.5, 0.4, 0.2, 0.1], m = 1,
        precond = true, GF=GF1!
    )
    exnum=ExNum(GF)
    replhead=manefesto(exnum, precond, tau, "tau", "alpha")
    println(replhead)
    nalpha = length(pvec)
    labelarray = String[]
    for ist in 1:nalpha
        push!(labelarray, string(pvec[ist]))
    end
    avals = Vector{Array}(undef, nalpha)
    rvals = Vector{Array}(undef, nalpha)
    for ialpha in 1:nalpha
        (rhist, errhist) = aa_nonlip(
            n; nu = nu, alpha = pvec[ialpha], tol=tol,
            maxit = maxit, tau = tau, m = m, precond=precond,
            GF=GF
        )
        avals[ialpha] = errhist
        rvals[ialpha] = rhist
    end
    ts = L"\tau_0"
    as = L"\alpha"
    ptex = "" * ts * " = "* string(tau) * "; vary $as"
    ctex = "Anderson($m)"
    Plot_NL(avals, rvals, as, ptex, ctex, labelarray,exnum)
    return
end
