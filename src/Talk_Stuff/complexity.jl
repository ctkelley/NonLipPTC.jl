function complexity(h, alpha, epsilon)
    M = h^-2
    p1 = log(M^(1.0 + alpha) / (4.0 * epsilon))
    t1 = 2.0 * (1.0 - alpha) / (1.0 + alpha)
    p2 = M * epsilon^(-t1)
    tau = 1.0 / p2
    itc = Int(floor(p1 * p2))
    println("tau = $tau, iterations = $itc")
    return (itc, tau)
end
