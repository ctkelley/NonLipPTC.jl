function testP(n; maxit = 10000, tau = 0.001, taunp = 0.001, algnum = 2, tol = 1.0e-3)
    algfuns = [alg1_test, alg3_test]
    aoutp = algfuns[algnum](n; pcond = true, maxit = maxit, tau = tau, epsilon = tol)
    aoutnp = algfuns[algnum](n; tau = taunp, maxit = maxit, epsilon = tol)
    #aoutp=alg3_test(n;pcond=true, maxit=maxit, tau=tau)
    #aoutnp=alg3_test(n;tau=tau, maxit=maxit)
    h = 1.0 / (n + 1.0)
    x = h:h:(1.0 - h)
    figure(1)
    semilogy(aoutp.reshist, "k-", aoutnp.reshist, "k--")
    legend(["precond", "no precond"])
    figure(3)
    subplot(1, 2, 1); plot(x, aoutp.sol); xticks(0:0.1:1.0); ylim(bottom = 0)
    title("precond")
    subplot(1, 2, 2); plot(x, aoutnp.sol); xticks(0:0.1:1.0); ylim(bottom = 0)
    title("no precond")
    println(norm(aoutp.sol - aoutnp.sol) / norm(aoutnp.sol))
    return println(minimum(aoutp.sol), "   ", minimum(aoutnp.sol))
end

function aa_test(
        n, u0 = nothing; maxit = 10000, rtol = 1.0e-1, atol = 1.0e-1, beta = 0.001,
        lambda = 200.0, p = 0.1, m = 1
    )
    pdata = setdata(n, lambda, p)
    if (u0 == nothing)
        u0 = ones(n)
    end
    Vstore = zeros(n, 3 * m + 4)
    aout = aasol(
        GF1d!, u0, m, Vstore; beta = beta, pdata = pdata, maxit = maxit,
        rtol = rtol, atol = atol
    )
    return aout
end

function HF1d!(G, u, pdata)
    beta = 0.001
    G .= GF1d!(G, u, pdata)
    G .= (1.0 - beta) * u + beta * G
    return G
end
