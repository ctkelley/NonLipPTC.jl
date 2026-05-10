function fobj(u, pdata)
    #
    # objective function
    #
    chat = pdata.bvec + pdata.rhs_eg1
    p = pdata.alpha
    nu = pdata.nu
    fcons = dot(u, chat)
    upp1 = (proj0.(u)) .^ (p + 1.0)
    fnlv = (nu / (p + 1.0)) * upp1
    fnl = sum(fnlv)
    D2 = pdata.D2
    fl = 0.5 * (u' * D2 * u)
    fobj = fl + fnl - fcons
    return fobj
end
