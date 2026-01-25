function alg0_test(n; p=.1, lambda=100.0, maxit=1000, tau=.01, 
epsilon=.001, pcond=false)
pdata=setdata(n, lambda, p)
u0=ones(n)
R=ones(n)
if pcond
aout=alg0(u0, fobj1D, POneD, proj, pdata, R, tau, maxit)
else
aout=alg0(u0, fobj1D, OneD, proj, pdata, R, tau, maxit)
end
return aout
end

function alg1_test(n; p=.1, lambda=100.0, maxit=1000, tau=.01, 
epsilon=.001, pcond=false)
pdata=setdata(n, lambda, p)
u0=ones(n)
R=ones(n)
if pcond
aout=alg1(u0, fobj1D, POneD, proj, pdata, R, tau, maxit)
else
aout=alg1(u0, fobj1D, OneD, proj, pdata, R, tau, maxit)
end
return aout
end

function alg3_test(n; p=.1, lambda=100.0, maxit=1000, tau=.01, 
epsilon=.001, pcond=false)
pdata=setdata(n, lambda, p)
u0=ones(n)
R=ones(n)
if pcond
aout=alg3(u0, fobj1D, POneD, proj, pdata, R, tau, epsilon, maxit)
else
aout=alg3(u0, fobj1D, OneD, proj, pdata, R, tau, epsilon, maxit)
end
return aout
end

