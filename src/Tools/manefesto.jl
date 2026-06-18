function manefesto(exnum, precond, fval, fname, vname)
Estr="Ex "*string(exnum)*"; "
Pstr="Precond = "*string(precond)*"; "
Fstr="Fixed "*fname* " = "*string(fval)*"; vary "*vname
Manefesto=Estr * Pstr * Fstr
end

