function build_problem3(n; c=1.0, k=1.0, alpha=.5, precond=false, 
        tau0=.5)
exsol(x) = sin(x)*exp(x)
u0=ones(n)
uc=copy(u0)
up=copy(u0)
L1=Lap1d(n)
D1=FD1(n)
h=1.0/(n+1); h=1.0/(n+1.0); x=h:h:1.0-h;
uex=exsol3.(x)
rhs = L1*uex + c*D1*uex + k*max.(uex,0).^alpha
pdata=(c=c, alpha=alpha, exsol=exsol, precond=precond,
       L1=L1, D1=D1, rhs=rhs, uex1d=uex, tau0=tau0, u0=u0, k=k)
end

function fobj3(u,pdata)
fobj3=1.0
end

function exsol3(x)
exs=1.0+2.0*sin(4.0*pi*x)*exp(x)
return exs
end

function FEX3(u,pdata)
c=pdata.c; alpha=pdata.alpha; precond=pdata.precond
L1=pdata.L1; D1=pdata.D1; rhs=pdata.rhs; k=pdata.k
n=length(u)
bcfix=sbc(n) + c* sbd(n)
#F3U = L1*u + c*D1*u + k*(max.(u,0).^alpha) - rhs + bcfix
F3U = L1*u + c*D1*u + k*(max.(u,0).^alpha) + bcfix 
F3P = L1\F3U
precond ? F3=F3P : F3=F3U
return F3
end

#
# Easy to understand Fixed point map
#
function GF1e3!(G, u, pdata)
    n = length(u);
    h=1.0/(n+1);
    tau0 = pdata.tau0
    precond=pdata.precond
    precond ? (tau = tau0) : (tau = tau0*h)
    G .= projp(u - tau*FEX3(u,pdata))
    return G
end

#
# Fixed point map for two stage method
#
function GF2e3!(G,u,pdata)
    n = length(u); 
    h=1.0/(n+1);
    tau0 = pdata.tau0
    precond=pdata.precond
    precond ? (tau = tau0) : (tau = tau0*h)
    y = projp(u - tau*FEX3(u,pdata))
    G .= projp(u - tau*FEX3(y,pdata))
    return G
end


function FD1(n)
h=1.0/(n+1.0)
ssdiag=ones(n-1)/(2.0*h)
updiag = Pair(1, ssdiag)
lowdiag = Pair(-1, -ssdiag)
D1=spdiagm(lowdiag, updiag)
return D1
end

function tauev1d(tau0,n,pdata)
    precond=pdata.precond
    h = 1.0 / (n + 1.0)
    precond ? tau = tau0 : (tau = tau0 * h * h)
    return tau
end 
