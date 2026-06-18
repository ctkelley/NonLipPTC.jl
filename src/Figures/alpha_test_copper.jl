function alpha_test_copper(
        n; nu = 1.0, tau=.1, tol=1.e-12, precond=false, 
        maxit = 20000, pvec = [0.5, 0.4, 0.2, 0.1], 
        algfun=alg1e1
    )
    exnum=ExNum(algfun)
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
       aout = algfun(
            n; tau0 = tau, alpha=pvec[ialpha],
            precond = precond, maxit = maxit
        )
        avals[ialpha] = aout.errhist
        rvals[ialpha] = aout.reshist
    end
    ts = L"\tau_0"
    as = L"\alpha"
    ptex = "" * ts * " = "* string(tau) * "; vary $as"
    precond ? (cval = "on") : (cval = "off")
    ctex = "Preconditioning $cval"
    Plot_NL(avals, rvals, as, ptex, ctex, labelarray,exnum)
    return
end
