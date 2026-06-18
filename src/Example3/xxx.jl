function build_problem3(n; alpha=.5, tau0=.1, precond=true, u0=[])
h=1.0/n;
X=.5*h:h:(1.0-.5*h);
W=h
K=[sin( 2*pi*(x - y)) for x in X, y in X]
K .*= h
ux3 = uexact3.(X)
nl3=nl_term3.(ux3,alpha)
rhs3 = ux3 - K*nl3 
if length(u0)==0
u0=ones(n)
end

pdata=(X=X, W=W, K=K, rhs3=rhs3, alpha=alpha, 
       uex1d=ux3, ux3=ux3, u0=u0, precond=true)
end

function nl_term3(u, alpha)
   nl=max(u,0.0)^alpha
   return nl
end

function uexact3(x)
#   ux3 = x*(1.0-x)
   ux3 = x.*sqrt.(x).*cos.(10.0*x)
   return ux3
end

function FEX3(u, pdata)
alpha=pdata.alpha
rhs=pdata.rhs3
K=pdata.K
nl3=nl_term3.(u,alpha)
rhs3=pdata.rhs3
FG3 = u - K*nl3 - rhs3
return FG3
end 

function fobj3(u, pdata)
   return 1.0
end

function test3(pdata, tau)
ux3=pdata.ux3
u0=pdata.u0
alpha=pdata.alpha
n=length(u0)
u0=ones(n)
up=copy(u0)
uc=copy(u0)
R = FEX3(u0,pdata)
ernrm = norm(u0-ux3,Inf)
rnrm=norm(R,Inf)
println(ernrm,"  ",rnrm)
for i=1:100
    up .= uc - tau*R
    ernrm=norm(up - ux3,Inf)
    R=FEX3(up,pdata)
    rnrm=norm(R,Inf)
    uc .= up
    println(ernrm,"  ",rnrm)
end
end
