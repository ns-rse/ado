*! version 1.4 Apr, 2003 DGC
program define gtab, rclass
  version 7.0
  syntax varlist(min=2 max=2) [if] [,Gen(string) REF(string) RECode]
  marksample touse
  tokenize `varlist'
  local a1 `1'
  local a2 `2'
  local alab : value label `a1'
  if "`ref'"!="" {
    if "`gen'"=="" {
      di in red "ref() option only relevant with gen() option"
      exit
    }
    if substr("`ref'",1,1)=="_" & "`ref'"!="_most" {
      if "`alab'"=="" {
        di in red "Reference by label requires value labels for alleles"
        exit
      }
      local lref=substr("`ref'", 2, .)
    }
  }
  preserve
  quietly keep if `touse'
  keep `a1' `a2'
  tempvar ln al
  gen long `ln'=_n
  quietly expand 2
  sort `ln'
  quietly by `ln': gen `al' = cond(_n==1,`1', `2')
  sort `al'
  quietly by `al': keep if _n==1
  local i = 0
  local codes
  local rule 
  while `i'<_N {
    local i = `i' + 1
    local cd = `al'[`i']
    if `cd'!=. {
      local codes "`codes' `cd'"
      local rule "`rule' `cd'=`i'"
    }
  }
  restore
  local k: word count `codes'
  local most = 0
  local alldum 
  local vref
  tokenize `codes'
  local i = 0
  di
  di "---------------------------------------------------------------------------------"
  di _skip(25) "    Total  " _continue
  di "     Homozygots   " _continue
  di "    Heterozygots   " _continue
  di "     Z"
  di  _skip(19) "Allele"  _skip(2) "frequency" _continue
  di "  Obsvd  ( Expctd)" _continue
  di "  Obsvd  ( Expctd)    (HWE)"
  di "---------------------------------------------------------------------------------"
  while "`1'"!="" {
    local i = `i'+1
    if "`gen'"=="" {
      tempvar dummy
    }
    else {
      local dummy "`gen'`1'"
    }
    local alldum "`alldum' `dummy'"
    quietly gen `dummy'= cond((`a1'!=. & `a2'!=.), /*
              */ (`a1'==`1') + (`a2'==`1'), .) if `touse'==1
    quietly count if `dummy'==0 
    local f0 = r(N)
    quietly count if `dummy'==1 
    local f1 = r(N)
    quietly count if `dummy'==2 
    local f2 = r(N)
    local f = 2*`f2' + `f1'
    local nc = `f0' + `f1' + `f2'
    local e0 = ((2*`f0' + `f1')^2)/(4*`nc')
    local e2 = ((2*`f2' + `f1')^2)/(4*`nc')
    local e1 = `nc' - `e0' - `e2'
    local z = (`f2'-`e2')*sqrt(1/`e0' + 1/`e2' + 4/`e1')
    if "`alab'"!="" {
      local lab : label `alab' `1'
      if "`gen'"!="" {
        label variable `dummy' "`lab'"
      }
    }
    else {
      local lab ""
    }
    di %25s = "`1': `lab'"  _continue
    di %9.0f = `f' _skip(2) _continue
    di %7.0f = `f2' "  (" %7.1f = `e2' ")" _continue
    di %7.0f = `f1' "  (" %7.1f = `e1' ")" _continue
    di %9.3f = `z'
    if "`ref'"!="" {
      if "`lref'"!="" {
        if "`lref'"=="`lab'" {
          local vref "`dummy'"
        }
      }
      else if "`ref'"!="_most"{
        if "`ref'" == "`1'" {
          local vref "`dummy'"
        }
      }
      else if `f'>`most' {
        local vref "`dummy'"
        local most = `f'
      }
    } 
    mac shift
  }
  quietly kappa `alldum' if `touse'
  local kappa =  r(kappa)
  local zval = r(z)
  local pval = chiprob(1, (`zval')^2)
  di "---------------------------------------------------------------------------------"
  di
  di "Global kappa statistic for Hardy-Weinberg equilibrium = "%6.3f = `kappa'
  di "                                             (Z-value = "%6.3f =`zval' ")"
  di "                                             (p-value = "%6.4f =`pval' ")"
  return scalar p_val = `pval'
  return scalar z = `zval'
  return scalar kappa = `kappa'
  if ("`gen'"!="") {
    di
    if "`ref'"!="" {
      if "`vref'" != "" {
        tokenize `alldum'
        local dum
        while "`1'"!="" {
          if "`1'"!="`vref'" {
            local dum "`dum' `1'"
          }
          else {
            drop `1'
          }
          mac shift
        }
        return local ivars `dum'
        di "Indicator variables generated: `dum'"
        di "(Variable `vref', corresponding to reference category, has been dropped)"
      }
      else {
        di "Reference category, `ref', not found"
        di "Indicator variables generated: `alldum'"
        return local ivars `alldum'
      }
    }
    else {
      di "Indicator variables generated: `alldum'"
      return local ivars `alldum'
    }
  }
  else {
    drop `alldum'
  }
  if "`recode'"!="" {
    quietly recode `a1' `rule'
    quietly recode `a2' `rule'
    di "Alleles have been recoded to consecutive integers"
  }
end

