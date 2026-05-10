function alg3e1(
        n; nu = 0.5, alpha = 0.5, tau0 = 0.1, maxit = 20000,
        epsilon = 1.0e-6
    )
    pdata = build_problem(n, uefun_ex12d; alpha = alpha, nu = nu)
    u0 = pdata.u0
    R = FEX1(u0, pdata)
    N = length(u0)
    mu = pdata.mu
    GDE2 = GD_Prob(fobj, FEX1, proj0, pdata, u0, mu)
    tau = tauev(tau0, N)
    algout = alg3(GDE2, R, tau, epsilon, maxit)
    reshist = algout.reshist
    errhist = algout.errhist
    return alg2out = (reshist = reshist, errhist = errhist)
end

function alg3e2(
        n; alpha = 0.5, p = 1.5, tau0 = 0.1, delta = 20.0, maxit = 20000,
        epsilon = 1.0e-6
    )
    pdata = build_problem2(n; p = p, alpha = alpha, delta = delta)
    u0 = pdata.u0
    R = FEX2(u0, pdata)
    N = length(u0)
    tau = tauev(tau0, N)
    mu = pdata.mu
    GDE2 = GD_Prob(fobj2, FEX2, proj2, pdata, u0, mu)
    algout = alg3(GDE2, R, tau, epsilon, maxit)
    reshist = algout.reshist
    errhist = algout.errhist
    sol = algout.sol
    return alg2out = (sol = sol, reshist = reshist, errhist = errhist)
end
