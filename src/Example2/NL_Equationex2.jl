#
# Residual
#
function FEX2(u, pdata)
    #
    # Problem parameters
    #
    N = length(u); n = Int(sqrt(N)); h = 1.0 / (1.0 + n)
    precond=pdata.precond
    alpha=pdata.alpha
    p=pdata.p
    delta=pdata.delta
    #
    # boundary data
    #
    bvec=pdata.bvec
    #
    # negative Laplacian
    #
    D2 = pdata.D2
    #
    FG1 = D2*u + nlterm_ex2.(u, p, alpha, delta) - bvec
    v1 = reshape(FG1, (n,n))
    fg21 = fish2d(v1,pdata.fdata)
    FG2 = reshape(fg21, (n*n,1))
#    FG2 = D2\FG1
    precond ? (FG=FG2) : (FG=FG1)
    return FG
end

#
# Easy to understand Fixed point map
#
function GF1e2!(G, u, pdata)
    N = length(u); n = Int(sqrt(N));
    h=1.0/(n+1); h2=h*h
    tau0 = pdata.tau0
    precond=pdata.precond
    precond ? (tau = tau0) : (tau = tau0*h2)
    G .= proj2(u - tau*FEX2(u,pdata))
    return G
end
#   
# Fixed point map for two stage method
#   
function GF2e2!(G,u,pdata)
    N = length(u); n = Int(sqrt(N));
    h=1.0/(n+1); h2=h*h
    tau0 = pdata.tau0
    precond=pdata.precond
    precond ? (tau = tau0) : (tau = tau0*h2)
#    FC1 = FEX2(u, pdata)
    y = proj2(u - tau*FEX2(u,pdata))
    G .= proj2(u - tau*FEX2(y,pdata))
end


#
# This is the nonlinear term in the operator.
# the parameter delta depends on the example
#
function nlterm_ex2(u, p, alpha, delta)
    nlout1=(delta*abs.(u).^alpha)*sign.(u)
    nlout2=abs.(u).^(p-1) .* u
    nlout = nlout1 - nlout2
    return nlout
end

#
# objective function
#
function fobj2(u, pdata)
p=pdata.p
alpha=pdata.alpha
delta=pdata.delta
fcons=dot(u,pdata.bvec)
# Linear term
D2=pdata.D2
fl = 0.5*(u'*D2*u)
#
fnl1=(delta/(1.0+alpha)) * sum(abs.(u).^(alpha+1.0))
fnl2=(1.0/(p+1)) * sum(abs.(u).^(p+1))
fobj = fl-fcons + fnl1 - fnl2
end


