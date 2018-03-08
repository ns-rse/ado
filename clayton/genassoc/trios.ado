*! version 1.5, 8 Oct 2003, DGC
program define trios 
  version 7.0
  syntax varlist(min=2 max=2) [, /*
               */ Ped(varname) Id(varname) Mother(varname) Father(varname) /*
               */ SEX(varname) AFfect(varname) ACode(numlist) FIrst /*
               */ SAVing(string) REPlace POrigin MISS]
  if "`saving'"!="" | "`replace'"=="" {
    preserve
  }
  if "`miss'"!="" & "`porigin'"!="" {
    di in red "porigin and miss options cannot be used together"
    exit
  }
  tokenize  `varlist'
  local v1 "`1'"
  local v2 "`2'"
  local allab : value label `v1'
  tempfile labsav
  quietly label save `allab' using "`labsav'"
  /* Put genotype into order --- smallest allele first */
  tempvar swap work error last sv1 sv2 
  quietly {
    gen byte `swap' = `v1'>`v2'
    gen `work' = `v1'
    replace `v1' = `v2' if `swap'
    replace `v2' = `work' if `swap'
    drop `swap' `work'
  }
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
  if "`acode'"=="" {
    local acode 2
  }
  /* New variable names */
  local ma1 "mother_1"
  local ma2 "mother_2"
  local pa1 "father_1"
  local pa2 "father_2"
  local ch1 "child_1"
  local ch2 "child_2"
  local ma "Mother"
  local pa "Father"
  local ch "Child"
  local me "Mendelian"
  local mating "Mating"
  keep `ped' `id' `mother' `father' `affect' `varlist' 
  tempfile idord
  sort `ped' `id'
  quietly save "`idord'"
  di "Triads selected from offspring with variable `affect' == `acode'"
  quietly recode `affect' `acode' = 2 * = 1 if `affect'!=. 
  quietly keep if `affect'==2 & `v1'!=. & `v2'!=. & `mother' != . & `father'!=.
  local naff = _N
  di "Total number of genotyped affected offspring = `naff'" 
  tempvar o1 o2 m1 m2 f1 f2 mid fid tid gtname res
  rename `v1' `o1'
  rename `v2' `o2'
  rename `id' `tid'
  rename `mother' `mid'
  rename `father' `fid'

  rename `mid' `id'
  sort `ped' `id'
  quietly merge `ped' `id' using "`idord'"
  quietly drop if _merge==2
  quietly drop _merge `mother' `father'
  rename `v1' `m1'
  rename `v2' `m2'
  rename `id'  `mid'

  rename `fid' `id'
  sort `ped' `id'
  quietly merge `ped' `id' using "`idord'", update
  quietly drop if _merge==2
  quietly drop _merge `mother' `father'
  rename `v1' `f1'
  rename `v2' `f2'
  rename `id'  `fid'
  rename `tid' `id'
  rename `o1' `v1'
  rename `o2' `v2'
  drop `affect'
  rename `mid' `mother'
  rename `fid' `father'
  if "`miss'"=="" {
    quietly drop if `m1'==. | `m2'==. | `f1'==. | `f2'==. | `v1'==. | `v2'==.
    di "Number of triads with complete data = " _N
  }
  if "`first'"=="" {
    tempvar if1
    sort `ped' `id'
    by `ped': gen byte `if1' = (_n==1)
    quietly count if !`if1'
    local ndr = r(N)
    if `ndr'>0  & "`warn'"=="" {
      noi di " `ndr' triads dropped - first case restriction"
    }
    quietly drop if !`if1'
    drop `if1'
  }
  rename `m1' `ma1'
  rename `m2' `ma2'
  rename `f1' `pa1'
  rename `f2' `pa2'
  rename `v1' `ch1'
  rename `v2' `ch2'
  sort `ma1' `ma2' `pa1' `pa2' `ch1' `ch2'
  gen Count = 1
  collapse (count) Count, by(`ma1' `ma2' `pa1' `pa2' `ch1' `ch2')
  quietly {
    egen `ma' = gtype(`ma1' `ma2')
    egen `pa' = gtype(`pa1' `pa2')
    egen `ch' = gtype(`ch1' `ch2')
    gen `me' = /*
      */ ((`ch1'==`ma1'|`ch1'==`ma2'|`ma'==.)& /*
      */  (`ch2'==`pa1'|`ch2'==`pa2'|`pa'==.)) | /*
      */ ((`ch1'==`pa1'|`ch1'==`pa2'|`pa'==.)& /*
      */  (`ch2'==`ma1'|`ch2'==`ma2'|`ma'==.))
    gen `mating'= cond(`ma'<`pa',`ma'+(`pa'*(`pa'-1))/2,`pa'+(`ma'*(`ma'-1))/2)
    label define `me' 0 "No" 1 "Yes"
    label values `me' `me'
  }
  if "`porigin'"=="" {
    di "List of all triads observed:"
    gsort - `me' + `mating'
    list `mating' `ma' `pa' `ch' `me' Count, noobs
  }
  else {
    quietly drop if `ma'==`pa'
    quietly count if !`me'
    if r(N)>0 {
      di "The following non-Mendelian triads were observed:"
      list `mating' `ma' `pa' `ch' Count if !`me', noobs
      quietly keep if `me'
    }
    di
    quietly summ Count 
    di "Number of triads informative for parent/parent-of-origin effects = " /*
        */ r(sum)
    quietly summ Count if `ch1'!=`ch2'
    di "Number of these with heterozygous offspring = " r(sum)
    sort `mating' `ch'  `ma' `pa'
    di
    di "Sorted list of all informative triads observed:"
    list `mating' `ma' `pa' `ch' Count, noobs
  }
  if "`saving'"!="" {
    save `saving', `replace'
  }
  end





