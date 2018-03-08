program define rclogit, eclass
  version 7.0
  syntax varlist (min=2), GRoup(varname) [Robust CLuster(varname) OR *] 
  if "`cluster'"!="" {
    local robust "robust"
  }
  if "`robust'"!="" {
    if "`cluster'"=="" {
      local cluster "`group'"
    }
    if "`or'"!="" {
      local eform "eform(Odds Ratio)"
    }
    preserve
    sort `group' `cluster'
    tempvar nest
    by `group': gen byte `nest'= (`cluster'[1]==`cluster'[_N])
    quietly count if !`nest'
    if r(N)>0 {
      di in red "Groups are not nested within clusters"
      exit
    }
    drop `nest'
    quietly clogit `varlist', group(`group') `options'
    local nobs = e(N)
    local df = e(df_m)
    local vdep "`e(depvar)'"
    tempvar res
    tempname b V M R Wald npar
    matrix `b'=e(b)
    matrix `V'=e(V)
    quietly predict `res'
    tokenize `varlist'
    quietly replace `res' = `1' - `res'
    local indv : colnames `b'
    tokenize  `indv'
    local scores
    while "`1'"!="" {
      tempvar sc
      local scores "`scores' `sc'"
      gen double `sc' = `1'*`res'
      mac shift
    }
    sort `cluster'
    collapse (sum) `scores', by(`cluster')
    local nc = _N
    quietly mat accum `M' = `scores', noconstant
    mat `M' = `M' * `V'
    mat `R' = `V' * `M'
    mat `V' = syminv(`R')
    mat `Wald' = `b' * `V'
    mat `Wald' = `Wald' * `b''
    scalar `npar' = colsof(`b')
    local chi2 = `Wald'[1,1]
    estimates post `b' `R'
    estimates scalar chi2 = `chi2'
    estimates local chi2type "Wald"
    estimates scalar df_m = `df'
    estimates scalar Nclus= `nc'
    estimates local depvar "`vdep'"
    estimates local group "`group'"
    estimates local cluster "`cluster'"
    estimates local predict "clogit_p"
    estimates local cmd "rclogit"
    di 
    di "Conditional (fixed-effects) logistic regression" _continue
    di  _column(51) "Number of obs"  _column(67) "=" _column(70) %9.0g =`nobs'
    di "(robust standard errors)" _continue
    di  _column(51) "Number of clus"  _column(67) "=" _column(70) %9.0f =`nc'
    di  _column(51) "Wald chi2(`df')"  _column(67) "=" _column(70) %9.2f = /*
      */ `chi2'
    di  _column(51) "Prob > chi2"  _column(67) "=" _column(70) %9.4f = /*
      */ chi2tail(`df',`chi2')
    di 
    estimates display, `eform'
  }
  else {
    clogit `varlist', group(`group') `or' `options'
  }
  end
