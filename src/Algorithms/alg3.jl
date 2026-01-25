function alg3(u0, fobj, fgrad, proj, pdata, R, tau, epsilon, maxit)
    uex=pdata.uex1d
    E0=norm(u0 - uex)
    E = E0
    mu=pdata.mu
    reshist=Float64[]
    errhist=Float64[]
    R .= fgrad(u0,pdata)
    N0=norm(R)
    push!(reshist, 1.0)
    push!(errhist, E/E0)
    RX=copy(R)
    N=length(u0)
#    nu = tau0/N
#    nu = tauev(tau0,N)
    nu = tau
    eta = nu/(1.0 + nu)
    u = copy(u0)    
    w = copy(u0)    
    v = copy(u0)    
    z = copy(u0)    
    ux = copy(u0)
    fc = fobj(u,pdata)
    for ix=1:maxit
        v .= (1.0 - eta)*u + eta*proj(w)
        z .= proj(proj(w) - (nu/mu)*fgrad(v,pdata))
        ux .= (1.0 - eta)*u + eta*z
        RX .= fgrad(ux, pdata)
        ft = fobj(ux,pdata)
#
# I am taking out the requirement that the gradient norm
# decrease for testing purposes. Lei explained why I should
# do that.
#
#        if (norm(R) > norm(RX))
           u .= ux
           R .= RX
           fc=ft
           w .= (1.0 - eta)*w + eta*v - (eta/mu)*fgrad(v,pdata)
#        else
#           RX .= R
#           ux .= u
#           ft = fc
#        end    
        rrnrm=norm(R) / N0
        push!(reshist, rrnrm)
        E=norm(u - uex)
        push!(errhist, E/E0)
    end
    return (sol = u, reshist = reshist, errhist = errhist)
end

function alg3(GP::GD_Prob, R, tau, epsilon, maxit)
    u0=GP.u0
    fobj=GP.fobj
    fgrad=GP.fgrad
    pdata=GP.pdata
    proj=GP.projb
#    mu = GP.mu
aout=alg3(u0, fobj, fgrad, proj, pdata, R, tau, epsilon, maxit)
end

