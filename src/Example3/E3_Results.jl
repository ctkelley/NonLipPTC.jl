function alg1e3(n; alpha=.5, tau0=.5, c=5.0, maxit=10000, 
            precond=false, k=50.0, tol=1.e-6)
    u0=ones(n);
    pdata=build_problem3(n; c=c, alpha=.5, precond=precond, k=k);
    R = FEX3(u0, pdata);
    tau=tauev1d(tau0, n, pdata)
    GDE3=GD_Prob(fobj3, FEX3, projp, pdata, u0, 1.0)
    algout = alg1(GDE3, R, tau, maxit; tol=tol)
    reshist = algout.reshist
    errhist = algout.errhist
    sol = algout.sol
    alg1out = (sol = sol, reshist = reshist, errhist = errhist)
    return alg1out
end

function alg2ge3( n; nu = 0.5, alpha = 0.5, tau0 = 0.01, 
              maxit = 10000, precond = false, c=5.0, k=50.0, tol=1.e-6
    )
    pdata=build_problem3(n; c=c, alpha=.5, precond=precond, k=k);
    u0=ones(n)
    R = FEX3(u0, pdata)
    N = length(u0)
    tau = tauev1d(tau0, N, pdata)
    GDE3 = GD_Prob(fobj3, FEX3, projp, pdata, u0, 1.0)
    algout = alg2g(GDE3, R, tau, maxit; tol=tol)
    reshist = algout.reshist
    errhist = algout.errhist
    sol = algout.sol
    return alg2out = (reshist = reshist, errhist = errhist, sol = sol)
end

function alg_base3(n; alg=alg1,
    alpha=.5, tau0=.01, c=5.0, maxit=10000, precond=false, k=50.0)
    u0=ones(n);
    pdata=build_problem3(n; c=c, alpha=.5, precond=precond, k=k);
    R = FEX3(u0, pdata);
    tau=tauev1d(tau0, n, pdata)
    GDE3=GD_Prob(fobj3, FEX3, projp, pdata, u0, 1.0)
    algout = alg(GDE3, R, tau, maxit)
    reshist = algout.reshist
    errhist = algout.errhist
    sol = algout.sol
    alg1out = (sol = sol, reshist = reshist, errhist = errhist)
    return alg1out
end 

