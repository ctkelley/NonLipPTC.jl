function aa_test_copper(
        n; nu = 1.0, alpha = 0.5, tol=1.e-12,
        maxit = 20000, tauvec = [0.2, 0.1, 0.05, 0.01], m = 1,
        precond = true, GF = GF1!
    )
    exnum=ExNum(GF)
    replhead=manefesto(exnum, precond, alpha, "alpha", "tau")
    println(replhead)
    ntau = length(tauvec)
    labelarray = String[]
    for ist in 1:ntau
        push!(labelarray, string(tauvec[ist]))
    end
    avals = Vector{Array}(undef, ntau)
    rvals = Vector{Array}(undef, ntau)
    for itau in 1:ntau
        (rhist, errhist) = aa_nonlip(
            n; nu = nu, alpha = alpha, tol=tol,
            maxit = maxit, tau = tauvec[itau], m = m,
            precond = precond, GF=GF
        )
        avals[itau] = errhist
        rvals[itau] = rhist
    end
    as = L"\alpha"
    ts = L"\tau"
#    ptex = "" * as * " = $alpha"
    ptex = "" * as * " = " *string(alpha) * "; vary $ts"
    ctex = "Anderson($m)"
    Plot_NL(avals, rvals, as, ptex, ctex, labelarray,exnum)
    return
end
