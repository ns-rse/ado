*! version 1.5, Jan 15, 2002

program define htsubsets, rclass
  syntax [varlist] [fw aw iw pw] [, MIn(integer 1) MAx(integer 5) /*
     */ Dots(integer 10) Saving(string) Replace CHap]
  preserve
  tempvar wt cwt  hap use
  quietly 
  if "`weight'"!="" {
    gen `wt' `exp'
  }
  else {
    gen `wt' = 1
  }
  gen byte `use'=1
  markout `use' `varlist' `wt'
  quietly keep if `use'
  local loci
  foreach locus of local varlist {
    quietly count if `use' & `locus'!=1 & `locus'!=2
    if r(N) > 0 {
      di "`locus' is not coded 1 or 2 - this locus ignored"
    }
    else {
      local loci "`loci' `locus'"
    }
  }
  if "`loci'"=="" {
    di in red "No diallelic loci found"
    exit
  }
  /* Collapse into haplotype file */
  quietly {
    gsort `loci', gen(`hap')
    by `hap': replace `wt' = sum(`wt')
    by `hap': keep if _n==_N
    gen `cwt' = sum(`wt')
    replace `cwt' = `cwt'/`cwt'[_N]
  }
  /* Sort haplotypes in decreasing order of frequency */
  gsort -`wt'
  /* Do search */
  local nloc : word count `loci'
  if `max'>`nloc' {
    local max = `nloc'
  }
  if "`chap'"!="" {
    di "Searching subsets of size up to `max' of `nloc' loci"
    di "Only subsets which linearly identify the most common haplotypes "/*
       */ "will be considered as eligible"
  }
  else {
    di "Searching all subsets of size up to `max' of `nloc' loci"
  }
  if "`saving'"!="" {
    local pvars "Size"
    forvalues size = 1/`max' {
      local pvars "`pvars' Locus_`size'"
    }
    local pvars "`pvars' PDE_tot PDE_min RMD_tot RMD_max"
    tempname pfile
    tempfile tfile
    postfile `pfile' `pvars' using `tfile'
  }
  tempname index
  forvalues size = `min'/`max' {
    local maxt = 0
    local maxm = 0
    local mint = `nloc'
    local minm = 1
    local sz1 = `size'+1
    matrix one = J(`sz1', 1, 1)
    di
    di "Subsets of size `size': "
    local elig = 0
    local idot = 0 
    matrix  `index' = J(1, `size', 0)
    forvalues i = 1/`size' {
      matrix `index'[1,`i'] = `i'
    }
    local i = `size'
    while `i'>0 {
      if `i'==`size' {
        local subset
        local pvars "(`size')"
        forvalues j = 1/`size' {
          local ij = `index'[1,`j']
          local add : word `ij' of `loci'
          local subset "`subset' `add'"
          local pvars "`pvars' (`ij')"
        }
        /* Check whether this tags most common haplotypes */
        if "`chap'"!="" {
          mkmat `subset' in 1/`sz1', matrix(x)
          matrix define X = one,x
          cap matrix define x = inv(X)
          local cont = (_rc==0)
        }
        else {
          local cont 1
        }
        if `cont' {
          local elig = `elig' + 1
          quietly hapdiv `loci' [pw=`wt'], ht(`subset')
          local rm = r(rmd_max)
          local rt = r(rmd_tot)
          local wstr "`r(rmd_worst)'"
          local pm = r(p_min)
          local pt = r(p_tot)
          local wstp "`r(p_worst)'"
          if "`saving'" != "" {
            if `size'<`max' {
              local s1 = `size' + 1
              forvalues j = `s1'/`max' {
                local pvars "`pvars' (.)"
              }
            }
            local pvars "`pvars' (`pt') (`pm') (`rt') (`rm')"
            post `pfile' `pvars'
          }
          if `rt'<`mint' {
            local mint = `rt'
            local ssrt "`subset'"
          }
          if `rm'<`minm' {
            local minm = `rm'
            local ssrm "`subset'"
            local lminm "`wstr'"
          }
          if `pt'>`maxt' {
            local maxt = `pt'
            local sspt "`subset'"
          }
          if `pm'>`maxm' {
            local maxm = `pm'
            local sspm "`subset'"
            local lmaxm "`wstp'"
          }
        }
        if `dots'>0 {
          local idot = `idot'+1
          if `idot'==`dots' {
            di "." _continue
            local idot = 0
          }
        }
      }
      local imax = `nloc' - `size' + `i'
      if `index'[1,`i'] < `imax' {
        matrix `index'[1,`i'] = `index'[1,`i'] + 1
        local i1 = `i'+1
        forvalues j = `i1'/`size' {
          matrix `index'[1,`j'] = `index'[1,`i'] + `j'-`i'
        }
        local i = `size'
      }
      else {
        local i = `i'-1
      }
    }
    di
    if `elig'>0 {
      di "Number of eligible subsets = " `elig'
      di
      di %11s = "Criterion" %10s = "Value" _skip(12) "Subset"
      di %11s = "---------" %10s = "-----" _skip(12) "------"
      di %11s = "RMD (Total)" %10.3f = `mint' _skip(12) "`ssrt'"
      di %11s = "(Worst)" %10.3f = `minm' %-12s = " (`wstr')" "`ssrm'"
      di %11s = "PDE (Total)" %10.3f = `maxt' _skip(12) "`sspt'"
      di %11s = "(Worst)" %10.3f = `maxm' %-12s = " (`wstp')" "`sspm'"
    }
    else {
      di "No eligible subset could be found"
    }
  }
  if "`saving'"!="" {
    postclose `pfile'
    preserve
    use `tfile', clear
    gsort -Size -PDE_tot
    tempname vars
    label define `vars', add
    forvalues i = 1/`nloc' {
      local var : word `i' of `loci'
      label define `vars' `i' `"`var'"', add
    }
    forvalues i = 1/`max' {
      label values Locus_`i' `vars'
    }
    di
    save `saving', replace
    restore
  }
  return scalar p_total  = `maxt'
  return scalar p_worst  = `maxm'
  return scalar r_total  = `mint'
  return scalar r_worst  = `minm'
  return local  worst_p "`wstp'"
  return local  worst_r "`wstr'"
  return local pt_subset "`sspt'"
  return local pw_subset "`sspm'"
  return local rt_subset "`ssrt'"
  return local rw_subset "`ssrm'"
end
    
