function OneD(u, pdata)
bvec=pdata.bvec
L=pdata.L
p=pdata.p
lambda=pdata.lambda
nlt=lambda*(proj.(u).^p)
lt=L*u
resid=nlt+lt-bvec
return resid
end

function POneD(u,pdata)
bvec=pdata.bvec
L=pdata.L
p=pdata.p
lambda=pdata.lambda
nlt=lambda*(proj.(u).^p)-bvec
nltr=L\nlt
lt= u
resid=nltr+lt
return resid
end


function fobj1D(u, pdata)
L=pdata.L
p=pdata.p
lambda=pdata.lambda
# Boundary term
fcons = dot(u,pdata.bvec)
# Linear term
fl = .5*(u'*L*u) 
# Nonlinear term
upp1=(proj.(u)) .^ (1.0+p)
fnlv=(lambda/(1.0+p)) * upp1
fnl=sum(fnlv)
fobj=fl-fcons + fnl
return fobj
end

function proj(u)
proj=max.(u,0.0)
return proj
end

function setdata(n, lambda, p)
L = Lap1d(n)
h=1.0/(n+1.0)
bvec=zeros(n)
bvec[1] += h^(-2)
bvec[n] += h^(-2)
uex1d=2.0*ones(n)
mu = pi*pi
pdata=(L=L, lambda=lambda, p=p, bvec=bvec, uex1d=uex1d, mu=mu)
return pdata
end
