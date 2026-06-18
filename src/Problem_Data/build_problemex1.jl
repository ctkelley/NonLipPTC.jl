function build_problem(n, uex; alpha = 0.5, nu = 0.5, precond = false, 
     tau0=.5)
    #
    # build rhs and boundary conditions
    #
    fdata = fishinit(n)
    h = 1.0 / (1.0 + n)
    #mu=1.0/(4.0*h*h)
    mu = 2 * pi * pi
    N = n * n
    #rhs=ones(N)
    #rhs2d=reshape(rhs,(n,n))
    D2 = Lap2d(n)
    X = h:h:(1.0 - h)
    uex2d = zeros(n, n)
    uex2d .= [uex(x, y) for x in X, y in X]
    uex1d = reshape(uex2d, (N, 1))
    #
    # Now fix the boundary
    #
    bvec = zeros(N)
    bvec .= fix_rhs!(bvec, uex)
    #
    # Does initial iterate satisfy the bc?
    #
    u0 = D2 \ bvec
    u02d = reshape(u0, (n, n))
    #figure(1)
    #mesh(u02d)
    #figure(2)
    #mesh(u02d-uex2d)
    nlterm = nu * uex1d .^ alpha
    rhs_eg1 = D2 * uex1d - bvec + nlterm
    lapex = -[lapeval_ex12d(x, y) for x in X, y in X]
    lapex1d = reshape(lapex, (N, 1))
    rhs_exact = lapex1d + nlterm
    return (
        bvec = bvec,
        u0 = u0,
        uex1d = uex1d,
        rhs_eg1 = rhs_eg1,
        rhs_exact,
        D2 = D2,
        alpha = alpha,
        nu = nu,
        mu = mu,
        fdata = fdata,
        precond = precond,
        tau0 = tau0
    )
end
