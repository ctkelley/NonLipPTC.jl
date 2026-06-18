"""
rhs to impose bc for u''
"""
function sbc(n)
b2=zeros(n)
h=1.0/(n+1.0)
hm2=1.0/(h*h)
b2[1]=-hm2
b2[n]=-hm2
return b2
end

function sbd(n)
b1=zeros(n)
h=1.0/(n+1.0)
b1[1]=-.5/h
b1[n]=.5/h
return b1
end

function D1test(n, ft, ftp)
h=1.0/(n+1)
x=h:h:1.0-h;
fdt=ft.(x);
bcf=sbd(n)
D1=FD1(n)
fpe=ftp.(x)
fpd = D1*fdt + bcf
println(norm(fpe-fpd, Inf))
plot(fpe-fpd)
end

function D2test(n, ft, ftpp)
h=1.0/(n+1)
x=h:h:1.0-h;
fdt=ft.(x);
bcf=sbc(n)
D2 = Lap1d(n)
fppd=D2*fdt + bcf
fppe=ftpp.(x)
d2d=norm(fppd-fppe, Inf)
plot(fppd - fppe)
println(d2d)
end

