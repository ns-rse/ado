*! version 1.4, 25 April 2003, DGC
program define origin
  version 7.0
  syntax varlist(min=2 max=2) [, /*
               */ Ped(varname) Id(varname) Mother(varname) Father(varname) /*
               */ SEX(varname) AFfect(varname) ACode(numlist) /*
               */ noFIrst MAT GType ROBust CLuster(string) EST REF(string)]
  preserve
  if "`robust'"!="" & "`cluster'"=="" {
    local cl "n"
  }
  if "`cluster'"!="" {
    local cl = substr("`cluster'",1,1)
    if "`cl'"!="n" & "`cl'"!="p" {
      di in red "Illegal cluster option --- must be nuclear or pedigree"
      exit
    }
  }
  if "`first'"!="" {
    di "Warning: you have overridden the FIRST option"
    di "This analysis is not valid when families are ascertained for multiple cases"
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
  /* alleles */
  if "`allab'"=="" {
    _alleles `v1' `v2'
  }
  else {
    _alleles `v1' `v2', lab(`allab')
  }
  local alleles `s(alleles)'
  local anames `s(labels)'
  local amin `s(first)'
  local amax `s(last)'
  local nall : word count `alleles'
  /* Reference category */
  if "`ref'"=="" {
    local ref `s(major)' 
  }
  else {
    local i 0
    while `i'<`nall' {
      local i = `i'+1
      local code : word `i' of `alleles'
      local name : word `i' of `anames'
      if "`ref'"=="`code'"/"`ref'"=="`name'" {
        local ref `code'
        local i `nall'
      }
    }
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
  tempvar mg fg cg g1 g2 ones

  keep `ped' `id' `mother' `father' `affect' `varlist' 
  tempfile idord
  sort `ped' `id'
  quietly save "`idord'"
  di
  di "Triads selected have {bf: `affect' == `acode'} for offspring"
  quietly recode `affect' `acode' = 2 * = 1 if `affect'!=. 
  quietly keep if `affect'==2 & `mother'!=. & `father'!=. /*
              */ & `v1'!=. & `v2'!=.
  local naff = _N
  di
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
  quietly drop if `m1'==. | `m2'==. | `f1'==. | `f2'==. | `v1'==. | `v2'==.
  di "Number of triads with complete data = " _N
  if "`cl'"=="`n'" {
    tempvar nuc
    sort `ped' `mother' `father'
    by `ped' `mother' `father': gen `nuc' = (_n==1)
    quietly  replace `nuc'=sum(`nuc')
  }
  keep `ped' `id' `nuc'  `m1' `m2' `f1' `f2' `v1' `v2'
  gen byte `ones' = 1

  if "`allab'"!="" {
    quietly do "`labsav'"
    label values `m1' `allab'
    label values `m2' `allab'
    label values `f1' `allab'
    label values `f2' `allab'
    label values `v1' `allab'
    label values `v2' `allab'
    
    decode `m1', gen(`g1')
    decode `m2', gen(`g2')
    quietly replace `g1' = `g1'+"/"+`g2'
    encode `g1', gen(`mg')
    drop `g2'

    decode `f1', gen(`g1')
    decode `f2', gen(`g2')
    quietly replace `g1' = `g1'+"/"+`g2'
    encode `g1', gen(`fg')
    drop `g2'

    decode `v1', gen(`g1')
    decode `v2', gen(`g2')
    quietly replace `g1' = `g1'+"/"+`g2'
    encode `g1', gen(`cg')
    drop `g2'
  }
  else {
    gen str12 `g1' = string(`m1')+"/"+string(`m2')
    encode `g1', gen(`mg')
    drop `g1'
    
    gen str12 `g1' = string(`f1')+"/"+string(`f2')
    encode `g1', gen(`fg')
    drop `g1'
    
    gen str12 `g1' = string(`v1')+"/"+string(`v2')
    encode `g1', gen(`cg')
    drop `g1'
  }
  quietly drop if `mg'==`fg'
  di "Number of triads with parents who are dissimilar at the locus = " _N

  /* Calculate order of child genotype */

  tempvar ord1 ord2
  gen byte `ord1' = ((`v1'==`m1')|(`v1'==`m2')) & ((`v2'==`f1')|(`v2'==`f2'))
  gen byte `ord2' = ((`v1'==`f1')|(`v1'==`f2')) & ((`v2'==`m1')|(`v2'==`m2'))
  quietly count if !(`ord1'|`ord2')
  if r(N)>0 {
    di "The following are misinheritances and will be dropped"
    list `ped' `id' if  !(`ord1'|`ord2')
    quietly drop if !(`ord1'|`ord2')
  }
  sort `ped' `id'
  tempvar if1
  by `ped': gen byte `if1' = (_n==1)
  quietly count if !`if1'
  local ndr = r(N)
  if `ndr'>0 {
    di
    di "Warning: there are multi-case families: this analysis is suspect"
    if "`first'"=="" {
      di "Using only first case in such families: " _continue
      di "`ndr' triads dropped"
      quietly drop if !`if1'
    }
  } 
     
  /* Indicator variables for maternal origin */

  local povars
  if "`gtype'"=="" {
    foreach a1 of local alleles {
      if "`a1'"!="`ref'" {
        tempvar x1
        gen byte `x1' = cond(`v1'==`v2', 0, /*
          */ cond((`ord1'&(`v1'==`a1'))|(`ord2'&(`v2'==`a1')), +1, /*
          */ cond((`ord1'&(`v2'==`a1'))|(`ord2'&(`v1'==`a1')), -1, /*
          */ 0)))
        local povars "`povars' `x1'"
      }
    }
  }
  else {
    local done
    foreach a2 of local alleles {
      foreach a1 of local done {
        tempvar x1
        gen byte `x1' = cond((`v1'==`a1')&(`v2'==`a2'), cond(`ord1',+1,-1), 0)
        local povars "`povars' `x1'"
      }
      local done "`done' `a2'"
    }
  }
    
  /* Indicator variables for direct maternal effects */

  if "`mat'"!="" {
    local done
    foreach a2 of local alleles {
      local done "`done' `a2'"
      foreach a1 of local done {
        tempvar x1
        gen byte `x1'=((`m1'==`a1')&(`m2'==`a2')) - ((`f1'==`a1')&(`f2'==`a2'))
        local mgvars "`mgvars' `x1'"
      }
    }
  }

  /* Hopefully I will be able to put the following line in eventually ...

  logit `ones' `povars' `mgvars', nocons

  */
  di
  if "`cl'"=="" {
    di "Model with " _continue
    if "`mat'"!="" {
      di "maternal genotype but " _continue
    }
    di "no parental origin effects, "
    if "`mat'"!="" {
      quietly glm `ones' `mgvars', nocons fam(bin)
      local nll = e(ll)
      local ndf = e(df)
    }
    else {
      local nll = -_N*log(2)
      local ndf =  _N
    }
    di _column(8) "Log likelihood (" `ndf' " df) = " `nll'
  }
  else {
    di "Robust variance estimates requested, " _continue
    di "taking the units of clustering as " _continue
    if "`cl'"=="`n'" {
      di "nuclear families"
      local cl "cluster(`nuc')"
    }
    else {
      di "pedigrees"
      local cl "cluster(`ped')"
    }
    di
  }
  di "Model with " _continue
  di "parental origin " _continue
  if "`mat'"!="" {
    di "and maternal genotype " _continue
  }
  di "effects, "
  quietly glm `ones' `povars' `mgvars', nocons fam(bin) `cl'
  if "`cl'"=="" {
    local logl = e(ll)
    local df = e(df)
    di _column(8) "Log likelihood (" `df' " df) = " `logl'
    local chi = 2*(`logl'-`nll')
    local df = `ndf' - `df'
    local pval = chi2tail(`df', `chi')
    di "Likelihood ratio chi-squared test = " %7.3f = `chi' /*
      */ " (" `df' " df), p=" %6.2e = `pval'
  }
  quietly testparm `povars'
  local chi = r(chi2)
  local df = r(df)
  local pval = r(p)
  di "Wald test" _skip(24) " = " /*
    */ %7.3f = `chi' " (" `df' " df), p=" %6.2e = `pval'

  /* Estimates */
  if "`est'"!="" {
    matrix b = e(b)
    local fitted : colnames b
    local i = 0
    foreach pov of local povars {
      local i = `i' + 1
      local found_`i' 0
      foreach fit of local fitted {
        if "`pov'"=="`fit'"{
          local found_`i' 1
        }
      }
    }
    matrix V = e(V)
    di
    di "Estimated maternal origin effects"
    mat def R = J(`nall', `nall', 1)
    mat rowname R = `anames'
    mat roweq R = Maternal
    mat colname R = `anames'
    mat coleq R = Paternal
    
    local j  0
    di
    if "`gtype'"!="" {
      local ij = 0
      di %10s = "Genotype" _continue
    }
    else {
      local ja = 0
      di %10s = "Allele" _continue
    }
    di _column(15) "Effect" _column(26) "Ratio" _continue
    di _column(37) "S.E." _column(44) "z-value"
    while `j'<`nall' {
      local j = `j'+1
      local nj: word `j' of `anames'
      if "`gtype'"=="" {
        di %10s = "`nj'" _continue
        local aj: word `j' of `alleles'
        if "`aj'"!="`ref'" {
          local ja = `ja' + 1
          if `found_`ja'' {
            local var : word `ja' of `povars'
            matrix bj = - b[1,"`ones':`var'"]
            matrix Vj = V["`ones':`var'", "`ones':`var'"] 
            local beta =  - bj[1,1]
            local seb = sqrt(V[1,1])
            di %10.3g = `beta' %10.3g = exp(`beta') _continue
            di %10.3g = `seb' %10.3f = `beta'/`seb'
          }
          else {
            matrix bj = J(1,1,0)
            di _column(14) "Dropped" _column(26) "1.0"
          }
        }
        else {
          matrix bj = J(1,1,0)
          di _column(18) "Ref" _column(26) "1.0"

        }
        local ia = 0
      }
      local i  1
      while `i'<`j' {
        if "`gtype'"=="" {
          local ai: word `i' of `alleles'
          if "`ai'"!="`ref'" {
            local ia = `ia' + 1
            local var : word `ia' of `povars'
            if `found_`ia'' {
              matrix bij= bj+b[1,"`ones':`var'"]
            }
            else {
              matrix bij = bj
            }
          }
          else {
            matrix bij = bj
          }
        }
        else {
          local ni : word `i' of `anames'
          di %10s = "`ni'|`nj'" _continue
          local ij = `ij'+1
          if `found_`ij'' {
            local var: word `ij' of `povars'
            matrix bij = b[1,"`ones':`var'"]
            matrix Vij = V["`ones':`var'", "`ones':`var'"]
            local beta = bij[1,1]
            local seb = sqrt(Vij[1,1])
            di %10.3g = `beta' %10.3g = exp(`beta') _continue
            di %10.3g = `seb' %10.3f = `beta'/`seb'
          }
          else {
            matrix bij = J(1,1,0)
            di _column(14) "Dropped" _column(26) "1.0"
          }
        }
        local coef = bij[1,1]
        matrix R[`i',`j'] = exp(`coef')
        matrix R[`j',`i'] = 1/ R[`i',`j']
        if `coef'!=0 & "`egi'"=="" {
          local egi `i'
          local egj `j'
         }
        local i = `i'+1
     }
    }
    di
    di "Fitted ratio of risks for genotype i|j vs genotype j|i"
    matrix list R, noheader format(%7.3g)
    di
    di "Eg. genotype `egi'|`egj' (Maternal|Paternal) carries " _continue
    di %7.3g = R[`egi',`egj'] " times the risk of `egj'|`egi'"
  }   
  end

program define  _alleles, sclass
  version 7.0
  syntax varlist [, LABel(string)]
  local amin 9999
  local amax -9999
  local most 0
  local acode
  local aname
  local afreq
  tokenize `varlist'
  while "`1'"!="" {
    quietly summarize `1'
    if r(min)<`amin' {
      local amin = r(min)
    }
    if r(max)>`amax' {
      local amax = r(max)
    }
    mac shift
  }
  local code `amin'
  while `code'<=`amax' {
    tokenize `varlist'
    local count 0
    while "`1'"!="" {
      quietly count if `1'==`code'
      local count = `count'+r(N)
      mac shift
    }
    if `count'>0 {
      local acode "`acode' `code'"
      local afreq "`afreq' `count'"
      if "`label'"!="" {
        local name : label `label' `code'
        local aname "`aname' `name'"
      }
      else {
        local aname "`aname' `code'"
    }
    if `count'>`most' {
      local most `count'
      local maj `code'
    }
    local code = `code'+1
  }
  sreturn local alleles `acode'
  sreturn local labels `aname'
  sreturn local counts `afreq'
  sreturn local major `maj'
  sreturn local first `amin'
  sreturn local last `amax'
  end

