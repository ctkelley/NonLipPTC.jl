function alg2ge1(
        n; nu = 0.5, alpha = 0.5, tau0 = 0.1, maxit = 20000,
        precond = false
    )
    pdata = build_problem(n, uefun_ex12d; alpha = alpha, nu = nu, precond = precond)
    u0 = pdata.u0
    R = FEX1(u0, pdata)
    N = length(u0)
    #    tau=tau0/N
    tau = tauev(tau0, N, pdata)
    mu = 2.0 * pi * pi
    GDE1 = GD_Prob(fobj, FEX1, proj0, pdata, u0, mu)
    algout = alg2g(GDE1, R, tau, maxit)
    reshist = algout.reshist
    errhist = algout.errhist
    sol = algout.sol
    return alg1out = (reshist = reshist, errhist = errhist, sol = sol)
end

function alg2ge2(
        n; nu = 0.5, alpha = 0.5, tau0 = 0.1, maxit = 20000,
        precond = false
    )
    pdata = build_problem2(n; alpha = alpha, precond = precond)
    u0 = pdata.u0
    R = FEX2(u0, pdata)
    N = length(u0)
    #    tau=tau0/N
    tau = tauev(tau0, N, pdata)
    mu = 2.0 * pi * pi
    GDE1 = GD_Prob(fobj2, FEX2, proj0, pdata, u0, mu)
    algout = alg2g(GDE1, R, tau, maxit)
    reshist = algout.reshist
    errhist = algout.errhist
    sol = algout.sol
    return alg1out = (reshist = reshist, errhist = errhist, sol = sol)
end
