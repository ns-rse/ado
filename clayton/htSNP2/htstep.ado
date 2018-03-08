*! version 1.3, DGC, 24 Oct 2003

program define htstep, rclass
  syntax [varlist] [fw aw iw pw] [, MIn(integer 1) MAx(integer 0) /*
     */ CRIterion(string) INClude(varlist) EXClude(varlist) Until(real 0) /*
     */ RA(real 0.0) UP DOWN]
  preserve
  local dir "down"
  if "`up'"!="" {
    local dir "up"
  }
  else {
    local dir "down"
  }
  if `ra'<0 | `ra'>1 {
    di in red "Illegal ra() option"
    exit
  }
  local pqmin = `ra'*(1-`ra')
  local posscrit "pde r2 r2c r2a"
  if "`criterion'" == "" {
    local critg "r2"
    local mcrit "min"
  }
  else {
    local ncr : word count `criterion'
    if `ncr'>2 {
      di in red "Illegal criterion() option"
      exit
    }
    else {
      local critg : word 1 of `criterion' 
      if index("`posscrit'","`critg'")==0 {
        di in red "Criterion {bf:`criterion'} not recognized"
        exit
      }
      if `ncr'==2 {
        local mcrit : word 2 of `criterion'
        if "`mcrit'"!="mean" & "`mcrit'"!="min" {
          di in red "Criterion modifier {bf:`mcrit'} not recognized"
          exit
        }
      }
      else {
        local mcrit min
      }
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
  local nloc : word count `loci'
  if `nloc'==0 {
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
  if `max'==0 | `max'>=`nloc' {
    local max = `nloc'-1
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
  di "Criterion used: {bf:`mcrit'} of {bf:`critg'} " _continue
  di "in a step-`dir' algorithm"
  if `until'>0 {
    di "Search will be stopped when {bf:`mcrit'} of {bf:`critg'} " _continue
    di "exceeds {bf:`until'}"
  }
  local ninc : word count `include'
  if "`dir'"=="up" {
    local size 0
    local subset
  }
  else {
    local size `nloc'
    local subset "`loci'"
  }
  while 1 {
    if "`dir'"=="up" {
      local crit = -1
      foreach locus of local loci{
        local ss "`subset' `locus'"
        quietly {
          haptag `loci' `exclude' [pw=`wt'], ht(`include' `ss') ra(`ra')
        }
        if r(`critg'_`mcrit')>`crit' {
          local crit = r(`critg'_`mcrit')
          local add `locus'
        }
      }
      local size = `size' + 1
      local loci : subinstr local loci "`add'" "", word
      local subset "`subset' `add'"
      di "  Adding {bf:`add'}, `mcrit' of `critg' = " %6.4f = `crit'
      if `size' == `max' | (`until'>0 & `crit'>`until') {
        continue , break
      }
    }
    else {
      local crit = -1
      foreach locus of local subset{
        local ss : subinstr local subset "`locus'" "", word
        quietly haptag `loci' [pw=`wt'], ht(`include' `ss') ra(`ra')
        if r(`critg'_`mcrit')>`crit' {
          local crit = r(`critg'_`mcrit')
          local omit `locus'
        }
      }
      if `until'>0 & `crit'<`until' {
        continue, break
      }
      else {
        local size = `size' - 1
        local subset : subinstr local subset "`omit'" "", word
        di "  Omitting {bf:`omit'}, `mcrit' of `critg' = " %6.4f = `crit'
        if `size'==`min' {
          continue, break
        }
      }
    }
  }
  local ss
  foreach locus of local varlist {
    local subset : subinstr local subset "`locus'" "", word count(local ifin)
    if `ifin'==1 {
      local ss "`ss' `locus'"
    }
  }
  haptag `loci' `exclude' [pw=`wt'], ht(`include' `ss') ra(`ra')
  return scalar `critg'_`mcrit' = r(`critg'_`mcrit')
  return scalar size = `size' + `ninc'
  return local htsnps "`include' `ss'"
end






