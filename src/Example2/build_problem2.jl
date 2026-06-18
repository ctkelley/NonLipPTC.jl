function build_problem2(n; p=.5, alpha=.5, delta=2.0, precond=false, tau0=.1)
deltamin=(p/alpha) + 1.e-5
delta = max(delta, deltamin)
fdata = fishinit(n)
N=n*n;
h=1.0/(n+1.0);
X=h:h:(1.0-h);
D2=Lap2d(n);
mu = 2.0*pi*pi
#
# We are using the function identically one to play
# the role of the analytic solution so we can set
# the boundary data like we did in example 1
#
uex2d=zeros(n,n);
#uex2d .= [uexe2(x, y) for x in X, y in X];
uex1d = reshape(uex2d, (N, 1));
#
# Now fix the boundary and make u0 satisfy the bc.
#
bvec = zeros(N)
bvec .= fix_rhs!(bvec, uexe2)
u0 = D2\bvec
rhs_eg2 = zeros(N)
rhs_exact = zeros(N)
return(
        bvec = bvec,
        u0 = u0,
        uex1d = uex1d,
        rhs_eg2 = rhs_eg2,
        rhs_exact,
        D2 = D2,
        p = p,
        alpha = alpha,
        mu = mu,
        delta = delta,
        tau0=tau0,
        fdata=fdata,
        precond = precond
)
end

function uexe2(x,y)
#uex=1.0
uex=.5-sin(x)*sin(y)
return uex
end

function proj2(x)
r1=min.(x,1.0)
r2=max.(r1,-1.0)
return r2
end

