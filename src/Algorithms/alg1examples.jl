function alg1e1(
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
    algout = alg1(GDE1, R, tau, maxit)
    #    algout = alg1(u0, fobj, FEX1, proj0, pdata, R, tau, maxit)
    reshist = algout.reshist
    errhist = algout.errhist
    sol = algout.sol
    return alg1out = (reshist = reshist, errhist = errhist, sol = sol)
end

function alg1e2(
        n; alpha = 0.5, p = 1.5, tau0 = 0.1,
        delta = 20.0, maxit = 20000
    )
    pdata = build_problem2(n; p = p, alpha = alpha, delta = delta)
    u0 = pdata.u0
    R = FEX2(u0, pdata)
    N = length(u0)
    #    tau=tau0/N
    tau = tauev(tau0, N)
    mu = 2.0 * pi * pi
    GDE2 = GD_Prob(fobj2, FEX2, proj2, pdata, u0, mu)
    algout = alg1(GDE2, R, tau, maxit)
    #    algout = alg1(u0, fobj2, FEX2, proj2, pdata, R, tau, maxit)
    reshist = algout.reshist
    errhist = algout.errhist
    sol = algout.sol
    return alg2out = (sol = sol, reshist = reshist, errhist = errhist)
end


function tauev(tau0, N, pdata)
    precond = pdata.precond
    n = Int(sqrt(N))
    h = 1.0 / (n + 1.0)
    precond ? tau = tau0 : (tau = tau0 * h * h)
    return tau
end
