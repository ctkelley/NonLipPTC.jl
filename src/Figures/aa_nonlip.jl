function aa_nonlip(
        n; nu = 1.0, alpha = 0.5,
        maxit = 20000, tau = 0.02, tol = 1.0e-6, m = 1,
        precond=true, GF=GF1!, c=5.0, k=50.0
    )
    exnum=ExNum(GF)
    if exnum==1
    pdata = build_problem(n, uefun_ex12d; alpha = alpha, 
            nu = nu, precond = precond, tau0=tau)
    elseif exnum==2
    pdata = build_problem2(n; alpha = alpha, 
            precond = precond, tau0=tau)
    else
    pdata = build_problem3(n; alpha = alpha, 
            precond = precond, tau0=tau, c=c, k=k)
    end
    u0 = pdata.u0
    h = 1.0/(n+1.0); h2=h^2; 
    exnum == 3 ? N=n : N=n*n
    Vstore = zeros(N, 3 * m + 4)
#    beta = 1.0
#    precond ? (beta = tau) : (beta = tau*h2)
    aaout = aasol(
        GF, u0, m, Vstore; pdata = pdata, 
        maxit = maxit, rtol = tol, atol = tol*tol, keepsolhist = true
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
    return (rhist = aaout.history, errhist = evec, sol=aaout.solution)
end
