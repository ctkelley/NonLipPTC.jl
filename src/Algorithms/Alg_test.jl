function tau_test1(
        n; nu = 0.5, alpha = 0.5, maxit = 20000, algfun = alg1e1,
        tauvec = [0.2, 0.1, 0.05, 0.01]
    )
    ntau = length(tauvec)
    labelarray = String[]
    for ist in 1:ntau
        push!(labelarray, string(tauvec[ist]))
    end
    # ([string(tauvec[1]), string(tauvec[2]), string(tauvec[3]), string(tauvec[4])])
    #errmat=zeros(maxit+1,ntau)
    avals = Vector{Array}(undef, 4)
    for itau in 1:ntau
        aout = algfun(n; tau0 = tauvec[itau], maxit = maxit, alpha = alpha)
        avals[itau] = aout.errhist
        #    errmat[:,itau] .= aout.errhist
        #    semilogy(aout.errhist)
    end
    return plothist(avals, labelarray, "error")
end

function p_test1(
        n; nu = 0.5, tau0 = 0.1, algfun = alg1e1, maxit = 20000,
        pvec = [0.1, 0.2, 0.4, 0.5]
    )
    np = length(pvec)
    labelarray = String[]
    for ist in 1:np
        push!(labelarray, string(pvec[ist]))
    end
    #    labelarray = ([string(pvec[1]), string(pvec[2]), string(pvec[3]), string(pvec[4])])
    avals = Vector{Array}(undef, np)
    for ip in 1:np
        aout = algfun(n; tau0 = tau0, alpha = pvec[ip], maxit = maxit)
        avals[ip] = aout.errhist
    end
    return plothist(avals, labelarray, "error")
end

function p_test2(
        n; nu = 0.5, tau0 = 0.1, algfun = alg1e2, maxit = 20000,
        pvec = [0.1, 0.2, 0.3, 0.4]
    )
    labelarray = ([string(pvec[1]), string(pvec[2]), string(pvec[3]), string(pvec[4])])
    np = length(pvec)
    avals = Vector{Array}(undef, 4)
    for ip in 1:np
        aout = algfun(n; tau0 = tau0, alpha = pvec[ip], maxit = maxit)
        avals[ip] = aout.reshist
    end
    return plothist(avals, labelarray, "residual")
end


function res_test1(
        n; nu = 0.5, alpha = 0.5,
        tau0 = 0.1, algfun = alg3e1, maxit = 20000
    )
    aout = algfun(n; tau0 = tau0, alpha = alpha, maxit = maxit)
    #x=0:maxit-1
    semilogy(aout.errhist, "k-")
    semilogy(aout.reshist, "k--")
    legend(["error", "grad norm"])
    as = L"\alpha"
    return title(as * "=$alpha")
end

function Alg3_alpha_test(n = 15; tau0 = 0.1, pvec = [0.5, 0.6, 0.7, 0.8], maxit = 2000)
    labelarray = ([string(pvec[1]), string(pvec[2]), string(pvec[3]), string(pvec[4])])
    np = length(pvec)
    avals = Vector{Array}(undef, 4)
    for ip in 1:np
        aout = alg3e2(n; tau0 = tau0, alpha = pvec[ip], maxit = maxit)
        avals[ip] = aout.reshist
    end
    return plothist(avals, labelarray, "residual")
end
