*! version 1.3, DGC, 24 Oct 2003

program define htsearch, rclass
  syntax [varlist] [fw aw iw pw] [, MIn(integer 1) MAx(integer 0) /*
     */ CRIteria(string) INClude(varlist) EXClude(varlist) Dots(integer 10) /*
     */  UNtil(string) RA(real 0.0) Saving(string) Replace CHap FULL]
  preserve
  if "`full'"!="" & "`until'"=="" {
    di in red "full option requires until option"
    exit
  }
  if `ra'<0 | `ra'>1 {
    di in red "Illegal ra() option"
    exit
  }
  local pqmin = `ra'*(1-`ra')
  local posscrit "pde r2 r2c r2a"
  if "`criteria'" == "" {
    local criteria "`posscrit'"
  }
  else {
    foreach crit of local criteria {
      if index("`posscrit'","`crit'")==0 {
        di in red "Criterion {bf:`crit'} not recognized"
        exit
      }
    }
  }
  if "`until'"!="" {
    local n: word count `until'
    if `n'!=3 {
      di in red "Illegal {bf:until(criterion min|mean value)} option"
      exit
    }
    local until_c: word 1 of `until'
    if index("`criteria'","`until_c'")==0 {
      di in red "Illegal {bf:until(criterion ... ...)} option"
      exit
    }
    local until_m: word 2 of `until'
    if "`until_m'"!="mean" & "`until_m'"!="min" {
      di in red "Illegal {bf:until(... min|mean ...)} option"
      exit
    }
    local until_v: word 3 of `until'
    cap confirm number `until_v'
    if _rc!=0 {
      di in red "Illegal {bf:until(... ... value)} option"
      exit
    }
  }
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
      di "{bf:`locus'} is not coded 1 or 2 - this locus ignored"
    }
    else {
      quietly {
        summ `locus' [aw=`wt'], meanonly
        local pq = (r(mean)-1)*(2-r(mean))
      }
      if `pq'<`pqmin' {
        di "Minor allele of {bf:`locus'} is too rare - ignored"
        local exclude : subinstr local exclude "`locus'" "", word
        local exclude : subinstr local exclude "  " " "
      }
      else {
        local loci "`loci' `locus'"
      }
    } 
  }
  local wuse
  local wign
  foreach locus of local include {
    local loci : subinstr local loci "`locus'" "", word count(local there)
    local loci : subinstr local loci "  " " "
    if `there' {
      local wuse "`wuse' `locus'"
    }
    else {
      local wign "`wign' `locus'"
    }
  }
  local include "`wuse'"
  if "`wign'"!="" {
    di "Loci specified in include() option are not in varlist --- ignored: " /*
       */ "{bf: `wign'}"
  }
  local wuse
  local wign
  foreach locus of local exclude {
    local loci : subinstr local loci "`locus'" "", word count(local there)
    local loci : subinstr local loci "  " " "
    if `there' {
      local wuse "`wuse' `locus'"
    }
    else {
      local wign "`wign' `locus'"
    }
  }  
  local exclude "`wuse'"
  if "`wign'"!="" {
    di "Loci specified in exclude() option are not in varlist --- ignored: " /*
       */ "{bf: `wign'}"
  }
  if "`loci'"=="" {
    di in red "No additional diallelic loci found"
    exit
  }
  /* Collapse into haplotype file */
  quietly {
    gsort `include' `loci' `exclude', gen(`hap')
    by `hap': replace `wt' = sum(`wt')
    by `hap': keep if _n==_N
    gen `cwt' = sum(`wt')
    replace `cwt' = `cwt'/`cwt'[_N]
  }
  /* Sort haplotypes in decreasing order of frequency */
  gsort -`wt'
  /* Do search */
  local nloc : word count `loci'
  if `max'<`min' | `max'>`nloc' {
    local max = `nloc'
  }
  di "Searching subsets of sizes from `min' up to `max' of `nloc' loci:"
  di "{bf:`loci'}"
  if "`include'"!="" {
    di "(added to {bf:`include'})"
  }
  local nexcl : word count `exclude'
  if `nexcl'>0 {
    di "The following `nexcl' loci are to be predicted but cannot be htSNPs: "
    di "{bf:`exclude'}"
  }
  di "Criteria used: {bf:`criteria'}"
  if "`until'"!="" {
    di "Search will be stopped when {bf:`until_m'} of " _continue
    di "{bf:`until_c'} exceeds {bf:`until_v'}"
  }
  local ninc : word count `include'
  if "`saving'"!="" {
    local pvars "Size"
    forvalues size = 1/`ninc'{
      local pvars "`pvars' Inc_`size'"
    }
    forvalues size = 1/`max' {
      local pvars "`pvars' Add_`size'"
    }
    foreach crit of local criteria {                  
      local pvars "`pvars' `crit'_mean `crit'_min"
    }
    tempname pfile
    tempfile tfile
    postfile `pfile' `pvars' using "`tfile'"
  }
  tempname index
  forvalues size = `min'/`max' {
/*    local maxt = 0
    local maxm = 0
    local mint = `nloc'
    local minm = 1  */
    local sz1 = `size'+1
    foreach crit of local criteria {
      local `crit'_mean 0
      local `crit'_min 0
      local `crit'_worst
    }
    matrix one = J(`sz1', 1, 1)
    di
    di "Subsets of size `size': " _continue
    local idot = 0 
    matrix  `index' = J(1, `size', 0)
    forvalues i = 1/`size' {
      matrix `index'[1,`i'] = `i'
    }
    local i = `size'
    while `i'>0 {
      if `i'==`size' {
        local subset
        if "`saving'"!="" {
          local pvars "(`size')"
          forvalues j = 1/`ninc' {
            local pvars "`pvars' (-`j')"
          }
        }
        forvalues j = 1/`size' {
          local ij = `index'[1,`j']
          local add : word `ij' of `loci'
          local subset "`subset' `add'"
          if "`saving'"!="" {
            local pvars "`pvars' (`ij')"
          }
        }
        quietly {
          haptag `loci' `exclude' [pw=`wt'], ht(`include' `subset') ra(`ra')
        }
        if "`saving'" != "" {
          if `size'<`max' {
            local s1 = `size' + 1
            forvalues j = `s1'/`max' {
              local pvars "`pvars' (.)"
            }
          }
        }
        foreach crit of local criteria { 
          local cmin = r(`crit'_min)
          local cmean = r(`crit'_mean)
          local worst = r(`crit'_worst)
          if "`saving'" != "" {
            local pvars "`pvars' (`cmean') (`cmin')"
          }
          if `cmin' > ``crit'_min' {
            local `crit'_min  `cmin'
            local `crit'_worst `worst'
            local `crit'_min_ss `include' `subset'
          }
          if `cmean' > ``crit'_mean' {
            local `crit'_mean  `cmean'
            local `crit'_mean_ss `include' `subset'
          }
        }
        if "`saving'"!="" {
          post `pfile' `pvars'
        }
      }
      if `dots'>0 {
        local idot = `idot'+1
        if `idot'==`dots' {
          di "." _continue
          local idot = 0
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
    di
    di %13s = "Criterion" %10s = "Value" _skip(12) "Subset"
    di %13s = "---------" %10s = "-----" _skip(12) "------"
    foreach crit of local criteria {
      di %13s = "`crit'(Mean)" %10.3f = ``crit'_mean' _continue
      di _skip(12) "{bf:``crit'_mean_ss'}"
      di %13s = "`crit'(Min)" %10.3f = ``crit'_min' _continue
      di _skip(12) "{bf:``crit'_min_ss'}"
      di %13s = "(Worst locus)"  _skip(5) "({bf:``crit'_worst'})"
    }
    if "`full'"!="" {
      di
      di "Full listing of predictive efficiency"
      haptag `loci' `exclude' [pw=`wt'], ht(``until_c'_`until_m'_ss') ra(`ra')
    }
    if "`until'"!="" {
      if ``until_c'_`until_m'' >= `until_v' {
         continue , break
      }
    }
  }
  if "`saving'"!="" {
    postclose `pfile'
    use "`tfile'", clear
    gsort -Size
    tempname vars
    label define `vars', add
    forvalues i = 1/`ninc' {
      local var : word `i' of `include'
      label define `vars' -`i' `"`var'"', add
    }
    forvalues i = 1/`nloc' {
      local var : word `i' of `loci'
      label define `vars' `i' `"`var'"', add
    }
    forvalues i = 1/`ninc' {
      label values Inc_`i' `vars'
    }
    forvalues i = 1/`max' {
      label values Add_`i' `vars'
    }
    di
    save `saving', `replace'
    restore
  }
  foreach crit of local criteria {
    return scalar `crit'_mean = ``crit'_mean'
    return local `crit'_mean_ss ``crit'_mean_ss'
    return scalar `crit'_min = ``crit'_min'
    return local `crit'_worst ``crit'_worst'
    return local `crit'_min_ss ``crit'_min_ss'
  }
end
    





