function aa_nonlip(
        n; nu = 1.0, alpha = 0.5,
        maxit = 20000, tau = 0.02, tol = 1.0e-6, m = 3
    )
    precond = true
    pdata = build_problem(n, uefun_ex12d; alpha = alpha, nu = nu, precond = precond)
    u0 = pdata.u0
    N = n * n
    Vstore = zeros(N, 3 * m + 4)
    ntau = 1
    tauvec = [tau]
    labelarray = String[]
    for ist in 1:ntau
        push!(labelarray, string(tauvec[ist]))
    end
    aaout = aasol(
        GF!, u0, m, Vstore; pdata = pdata, beta = tau,
        maxit = maxit, rtol = tol, atol = tol, keepsolhist = true
    )
    lh = length(aaout.history)
    evec = zeros(lh)
    as = L"\alpha"
    ptex = "" * as * " = $alpha"
    ctex = "Anderson($m)"
    sole = pdata.uex1d
    for iv in 1:lh
        evec[iv] = norm(aaout.solhist[:, iv] - sole, Inf)
    end
    return (rhist = aaout.history, errhist = evec)
end

function aa_test_copper(
        n; nu = 1.0, alpha = 0.1, tol=1.e-12,
        maxit = 20000, tauvec = [0.2, 0.1, 0.05, 0.01], m = 1
    )
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
            maxit = maxit, tau = tauvec[itau], m = m
        )
        avals[itau] = errhist
        rvals[itau] = rhist
    end
    as = L"\alpha"
    ptex = "" * as * " = $alpha"
    ctex = "Anderson($m)"
    subplot(1, 2, 1)
    title(ptex)
    plothist(avals, labelarray, "error")
    subplot(1, 2, 2)
    title(ctex)
    return plothist(rvals, labelarray, "residual")
    return
end
