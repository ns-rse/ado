*! version 2.6, 9 Oct 2003, DGC
program define tdt, rclass
  version 7.0
  syntax varlist(min=2 max=2) [iweight/] [, /*
               */ Ped(varname) Id(varname) Mother(varname) Father(varname) /*
               */ Sex(varname) AFfect(varname) ACode(numlist) EMin(real 5) /*
               */ MOrigin POrigin PACode(numlist) MAle FEmale noWARN /*
               */ ROBust CLuster(string) ]
  preserve
  tempfile idord
  tempvar o1 o2 m1 m2 f1 f2 tid mid fid w sw mg fg o12 o21 inf use wt psel last
  tempvar chwt anum score var parent homoz
  if "`acode'"=="" {
    local acode 2
  }
  local origin "`porigin'`morigin'"
  if "`cluster'"!="" {
    local robust "robust"
  }
  if "`robust'"!="" {
    if "`cluster'"=="" {
      local cl "n"
    }
    else {
      local cl = substr("`cluster'", 1, 1)
    }
    if "`cl'"!="n" & "`cl'"!="p" & "`cl'"!="t" {
      di in red "Invalid cluster() option"
    }
  }
  if ("`porigin'"!="") & ("`morigin'"!="") {
    di in red "morigin and porigin options cannot be used together"
    exit
  }
  if ("`male'"!="") & ("`female'"!="") {
    di "male and female options cannot be used together"
    exit
  }
  di
  di "Analysis of"  _continue
  if "`morigin'"!="" {
    di " maternal"  _continue
  }
  if "`porigin'"!="" {
    di " paternal" _continue
  }
  di " transmissions to" _continue
  if "`male'"!="" {
    di " male" _continue
  }
  if "`female'"!="" {
    di " female" _continue
  }
  di " offspring with disease status codes `acode'"
  if "`pacode'"!="" {
    di "(restricted to transmissions from parents" _continue
    di " with disease status codes `pacode')"
  }
  if "`weight'" != "" {
    di "Weights will be used: `exp'"
    gen `wt' = `exp'
  }
  else {
    gen `wt' = 1
  }
    
  tokenize  `varlist'
  local v1 "`1'"
  local v2 "`2'"
  unab vars : _all
  tokenize `vars'
  if "`ped'"=="" {
    local ped "`1'"
  }
  if "`ped'"=="" {
    local ped "`1'"
  }
  if "`id'"=="" {
    local id "`2'"
  }
  if "`mother'"=="" {
    local mother "`4'"
  }
  if "`father'"=="" {
    local father "`3'"
  }
  if "`sex'"=="" {
    local sex "`5'"
  }
  if "`affect'"=="" {
    local affect "`6'"
  }

  keep `ped' `id' `mother' `father' `sex' `affect' `varlist' `wt'
  sort `ped' `id'
  quietly {
    if "`pacode'"!="" {
      gen `psel' = `affect'
      recode `psel' `pacode' = 1 * = 0
    }
    else {
      gen byte `psel' = 1
    }
    recode `affect' `acode' = 1 * = 0
    save "`idord'"
    gen byte `use' = (`affect' & `mother'!=. & `father'!=. /*
              */ & `v1'!=. & `v2'!=.)
    if "`male'"!="" {
      replace `use' = 0 if `sex'!=1
    }
    if "`female'"!="" {
      replace `use' = 0 if `sex'!=2
    }
    keep if `use'
    drop `affect' `psel' `use'
    rename `wt'  `chwt'
  }
  local naff = _N
  di
  di "Total number of genotyped affected offspring = `naff'"
  if "`weight'"!="" {
    quietly count if `chwt'==0
    local zerowt = r(N)
    if `zerowt'>0 {
      di "`zerowt' have zero weights"
      /* quietly drop if `chwt'==0 */
    }
  }
  quietly {
    rename `v1' `o1'
    rename `v2' `o2'
    rename `id' `tid'
    rename `mother' `mid'
    rename `father' `fid'

    rename `mid' `id'
    sort `ped' `id'
    merge `ped' `id' using "`idord'"
    drop if _merge==2
    drop _merge `mother' `father'
    rename `id' `mid'
    rename `v1' `m1'
    rename `v2' `m2'
    if "`porigin'"=="" {
      gen byte _PS_1 = `psel'     
    }
    else {
      gen byte _PS_1 = 0
    }
    drop `psel'
    
    rename `fid' `id'
    sort `ped' `id'
    merge `ped' `id' using "`idord'", update
    drop if _merge==2
    drop _merge `mother' `father'
    rename `v1' `f1'
    rename `v2' `f2'
    if "`morigin'"=="" {
      gen byte _PS_3 = `psel'     
    }
    else {
      gen byte _PS_3 = 0
    }
    drop `psel'

    rename `id' `father'
    rename `mid' `mother'
    rename `tid' `id'
    
    gen byte `mg' = (`m1'!=.) & (`m2'!=.)
    gen byte `fg' = (`f1'!=.) & (`f2'!=.)
    gen byte `o12' = ((!`mg')|`o1'==`m1'|`o1'==`m2') & /*
      */ ((!`fg')|`o2'==`f1'|`o2'==`f2')
    gen byte `o21' = ((!`mg')|`o2'==`m1'|`o2'==`m2') & /*
      */ ((!`fg')|`o1'==`f1'|`o1'==`f2')
    count if (!`o12')&(!`o21')
  }
  if r(N)>0 {
    if "`warn'"=="" {
      di
      di "Warning " r(N) " offspring with genotyping errors or ex-paternity :"
      list `ped' `id' if (!`o12')&(!`o21'), noobs
      di
    }
    quietly drop if (!`o12')&(!`o21')
  }
  if "`origin'"=="" {
    gen byte `inf' = ((`o12'&!`o21')|(`o21'&!`o12')|(`mg'&`fg'))
  }
  else {
    gen byte `inf' = ((`o12'&!`o21')|(`o21'&!`o12')|(`mg'&`fg'&(`o1'==`o2')))
  }
  quietly keep if `inf'
  gen byte `sw' = `o21'&!`o12'
  gen `w' = `o1'
  quietly replace `o1' = `o2' if `sw'
  quietly replace `o2' = `w' if `sw'
  quietly drop `w' `sw' `o12' `o21' `mg' `fg' `inf' 
  /* PA, PT transmission of maternal Codes (1&2) or paternal Codes(3&4) */
  gen byte _PT_1 = (`o1'==`m1')
  rename `m1' Allele1
  gen byte _PT_2 = (`o1'==`m2')
  rename `m2' Allele2
  gen byte _PT_3 = (`o2'==`f1')
  rename `f1' Allele3
  gen byte _PT_4 = (`o2'==`f2')
  rename `f2' Allele4
  quietly {
    drop `o1' `o2'
    replace _PS_1 = 0 if (Allele1==Allele2)
    replace _PS_3 = 0 if (Allele3==Allele4)
    gen byte _PS_2 = _PS_1
    gen byte _PS_4 = _PS_3
    reshape long Allele _PT_ _PS_, i(`ped' `id') j(_PC_)
    gen byte `parent' = cond(_PC_ < 3, 1, 2)
    drop if Allele==. | !_PS_ 
    drop _PS_
  }
  if _N==0 {
    di "No informative transmissions"
    exit
  }

  /* Find and save allele value label and numeric codes,
     Merge infrequent alleles */

  quietly {
    summarize Allele, meanonly 
    local allmiss = r(max) + 1
    di "`allmiss'"
    sort Allele
    by Allele: gen byte `last' = (_n==_N)
    local allab : value label Allele
    if "`allab'"=="" {
      local allab  _Allele
      label define `allab', modify
    }
    label define `allab' `allmiss' "Others", add 
    tempfile labsave
    label save `allab' using "`labsave'"
    by Allele: gen E =  _N/2
  }
  quietly count if `last' & E<`emin'
  if r(N)>0 {
    if "`warn'"=="" {
      di
      di "There are alleles with expected transmission frequencies < " `emin'
    }
    local ezero 0
    local thresh 0
    while `ezero'<`emin' {
      quietly summarize E if `last' & E>`thresh'
      if r(N)==0 {
        di in red "Minimum number of expected transmissions cannot be achieved"
        exit
      }
      local thresh = r(min)
      local ezero = `ezero' + `thresh'
    }
    if "`warn'"=="" {
     di "The following alleles have been grouped together:"
     list  Allele if `last' & E<=`thresh', noobs
    }
    quietly {
      replace Allele = `allmiss' if E<=`thresh'
      sort `ped' `id' `parent'
      by `ped' `id' `parent': gen byte `homoz' = /*
        */ (Allele[1]==`allmiss'  & Allele[2]==`allmiss') 
      count if `homoz'
      local nxh = r(N)/2
      drop if `homoz'
    }
    if `nxh'>0 & "`warn'"=="" {
      di "(`nxh' parents became effectively homozygous as a result)"
    }
  }
  di "Number of informative transmissions = "  _N/2
  if _N==0 {
    exit
  }

  quietly {
    sort Allele
    by Allele: gen byte `anum' = (_n==1)
    mkmat Allele if `anum', matrix(all)
    replace `anum' = sum(`anum')
    local numa = `anum'[_N]
    by Allele: replace `last' = (_n==_N)
    by Allele: gen O = sum(_PT_)
    by Allele: replace E =  _N/2
    gen _OmE = _PT_-0.5
    gen _Exp = 0.5
    by Allele: gen `score' = sum(`chwt'*_OmE)
    by Allele: gen `var' = sum((`chwt'*_Exp)^2)
    gen Chisq = `score'^2/`var' if `last' 
    gen P_val = chiprob(1, Chisq) if `last'
  }
  di
  di "Observed and expected transmission counts, " _continue
  di "and 1 df chisquared tests"
  quietly do "`labsave'"
  label value Allele `allab'
  list Allele O E Chisq P_val if `last', noobs
  drop Chisq P_val
    
  /* Global test */

  quietly {
    keep `ped' `id' `mother' `father' `parent' `anum' `chwt' _OmE _Exp
    sort `ped' `id' `parent' `anum'
    reshape wide _OmE _Exp, i(`ped' `id' `parent') j(`anum')
    mvencode _OmE* _Exp*, mv(0) override
    matrix vecaccum u = `chwt' _OmE*, nocons
    matrix accum V = _Exp* [iw=`chwt'^2], nocons
  }
  /* Calculate score test */
  matrix def d = vecdiag(V)
  matrix def V = 2*diag(d) - V
  matrix def A = syminv(V)
  local df = 0
  local i = 0
  while `i'<`numa' {
    local i = `i'+1
    if A[`i',`i']>0 {
      local df = `df'+1
    }
  }
  matrix c = u*A*u'
  local chisq = c[1,1]
  di
  di "Global chi-squared (" `df' " df) = " `chisq' /*
     */ ", P-value = " chiprob(`df', `chisq')

  /* Robust tests */

  if "`robust'"!="" {
    quietly drop _Exp*
    di
    di "Robust tests allowing for clustering at the " _continue
    if "`cl'"=="n" {
      local by "`ped' `mother' `father'"
      di "nuclear family " _continue
    }
    else if "`cl'"=="p" {
      local by "`ped'"
      di "pedigree " _continue
    }
    else { 
      local by "`ped' `id'"
      di "triad " _continue
    }
    di "level"
    sort `by'
    quietly {
      by `by': gen byte `last' = (_n==_N)
      local i 0
      while `i'<`numa' {
        local i = `i'+1
        by `by': replace _OmE`i' = sum(`chwt'*_OmE`i')
      }
      keep if `last'
    }
    di "(number of clusters = "  _N ")"
    quietly matrix accum V =  _OmE*, nocons /* dev */
    matrix V = V  /* * (_N/(_N-1)) */
    matrix def A = syminv(V)
    local i = 0
    quietly do "`labsave'"
    di
    di "1 df tests for each allele:"
    di
    di _column(5) "Allele" _column(15) "Chi-squared" _column(34) "P-value"
    local df = 0
    while `i'<`numa' {
      local i = `i'+1
      if A[`i',`i']>0 {
        local df = `df'+1
      }
      if V[`i',`i']>0 {
        local ai = all[`i',1]
        local anam : label `allab' `ai'
        local chisq = u[1,`i']^2/V[`i', `i']
        di %10s = "`anam'" _column(16) %10.4f = `chisq' /*
           */ _column(31) %10.3g =  chiprob(1, `chisq') 
      }
    }
    matrix c = u*A*u'
    local chisq = c[1,1]
    di
    di "Global chi-squared (" `df' " df) = " `chisq' /*
     */ ", P-value = " chiprob(`df', `chisq')
  }
  return local chi2type "Score"
  return scalar p_val = chiprob(`df', `chisq')
  return scalar df = `df'
  return scalar chi2 = `chisq'
  end


