function tau_test_copper(
        n; precond = false, nu = 1.0, alpha = 0.5,
        maxit = 20000, algfun = alg1e1,
        tauvec = [0.2, 0.1, 0.05, .01]
#        tauvec = [0.1, 0.05, 0.025, .01]
    )
    ntau = length(tauvec)
    labelarray = String[]
    for ist in 1:ntau
        push!(labelarray, string(tauvec[ist]))
    end
    # ([string(tauvec[1]), string(tauvec[2]), string(tauvec[3]), string(tauvec[4])])
    #errmat=zeros(maxit+1,ntau)
    avals = Vector{Array}(undef, ntau)
    rvals = Vector{Array}(undef, ntau)
    for itau in 1:ntau
        aout = algfun(
            n; tau0 = tauvec[itau],
            precond = precond, maxit = maxit, alpha = alpha
        )
        avals[itau] = aout.errhist
        rvals[itau] = aout.reshist
        #    errmat[:,itau] .= aout.errhist
        #    semilogy(aout.errhist)
    end
    as = L"\alpha"
    ptex = "" * as * " = $alpha"
    precond ? (cval = "on") : (cval = "off")
    ctex = "Preconditioning $cval"
    subplot(1, 2, 1)
    title(ptex)
    plothist(avals, labelarray, "error")
    subplot(1, 2, 2)
    title(ctex)
    return plothist(rvals, labelarray, "residual")
end
